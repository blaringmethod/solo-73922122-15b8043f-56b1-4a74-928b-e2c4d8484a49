"""Read-only inspection of a published lane's quarantine bucket.

Scans the published paired quarantine FASTQ files, validates records and
R1/R2 read-name pairing, builds the 16-base index from the first 8 bases of
each mate, and exactly counts index combinations per read pair.

Memory is bounded independently of the data: each 16-base index is appended
(raw, 16 bytes) to one of 256 hash-spill partition files chosen by CRC32,
and partitions are counted one at a time. Even a full 512 MiB lane of
distinct indices therefore fits in a 128 MiB container. Line reads are
windowed, so a single malformed giant line can never exhaust memory either.
Hot loops stay in C: buffered gzip reads, translate/min/max validation and
CRC32 partitioning.

The scan never writes to the published tree. Scratch space lives in a
locked ``work/inspect-<uuid>`` directory (the flock makes startup recovery
on the sibling instance skip it) and is always removed afterwards. Any
failure aborts the whole check: no partial statistics are produced.
"""
from __future__ import annotations

import fcntl
import gzip
import heapq
import io
import shutil
import struct
import uuid
import zlib
from collections import Counter
from pathlib import Path
from typing import BinaryIO, Iterator

from . import config

TOP_N = 20
MAX_SEQ_LINE = 1000   # bytes of sequence content per mate
MAX_QUAL_LINE = 1000  # bytes of quality content per mate
MAX_NAME = 256        # bytes of read name (header without the '@')
_MAX_PLUS_LINE = 1024  # memory-safety bound; the spec limits seq/qual/name only
_PARTITIONS = 256      # hash-spill partitions; must be a power of two
_FLUSH_BYTES = 32768   # partition buffer flush threshold (2048 indices)
_SIZE_BUDGET = config.MAX_TOTAL_UNCOMPRESSED  # 512 MiB decompressed per lane
_INDEX_LEN = 16        # 8 bases from R1 + 8 bases from R2

_DELETE_VALID = b"ACGTN"  # seq.translate(None, this) is empty iff all-ACGTN


class InspectError(Exception):
    """The check failed; the caller must get an error, never partial stats."""


class MissingDataError(InspectError):
    """Published quarantine files are missing (incomplete published data)."""


# ---------------- bounded FASTQ parsing ----------------

def _line(rl, limit: int, mate: str, index: int,
          errmsg: str) -> bytes | None:
    """Read one FASTQ line (content, newline stripped) with a hard cap.

    ``rl`` is a buffered readline; the read window is limit+3 bytes, so
    memory is capped regardless of the input. Any content longer than
    ``limit`` fails the whole check.
    """
    raw = rl(limit + 3)
    if raw == b"":
        return None
    if raw.endswith(b"\n"):
        content = raw[:-1]
        if content.endswith(b"\r"):
            content = content[:-1]
    else:
        content = raw  # last line without newline, or an over-limit fragment
    if len(content) > limit:
        raise InspectError(f"{mate}: record {index}: {errmsg}")
    return content


def _read_record(rl, mate: str,
                 index: int) -> tuple[bytes, bytes, bytes, bytes] | None:
    """Read and validate one record. Returns (name, seq, plus, qual); None at EOF.

    ``name`` is the header content without '@'. Limits: name <= 256 bytes,
    sequence and quality lines <= 1000 bytes each.
    """
    name_line = _line(rl, MAX_NAME + 1, mate, index,
                      f"read name exceeds {MAX_NAME} bytes")
    if name_line is None:
        return None
    seq = _line(rl, MAX_SEQ_LINE, mate, index,
                f"sequence line exceeds {MAX_SEQ_LINE} bytes")
    plus = _line(rl, _MAX_PLUS_LINE, mate, index,
                 f"'+' line exceeds {_MAX_PLUS_LINE} bytes")
    qual = _line(rl, MAX_QUAL_LINE, mate, index,
                 f"quality line exceeds {MAX_QUAL_LINE} bytes")
    if seq is None or plus is None or qual is None:
        raise InspectError(f"{mate}: record {index} is truncated")
    if not name_line.startswith(b"@"):
        raise InspectError(
            f"{mate}: record {index}: header does not start with '@'")
    name = name_line[1:]
    if not name.split():
        raise InspectError(f"{mate}: record {index} has an empty read name")
    if not plus.startswith(b"+"):
        raise InspectError(
            f"{mate}: record {index}: expected '+' separator line")
    if len(seq) == 0:
        raise InspectError(f"{mate}: record {index} has an empty sequence")
    if seq.translate(None, _DELETE_VALID):
        raise InspectError(
            f"{mate}: record {index} contains characters outside A/C/G/T/N")
    if len(seq) != len(qual):
        raise InspectError(
            f"{mate}: record {index} sequence/quality length mismatch")
    if min(qual) < 33 or max(qual) > 126:
        raise InspectError(
            f"{mate}: record {index} has quality bytes outside Phred+33 range")
    return name, seq, plus, qual


def _paired(r1_path: Path, r2_path: Path,
            ) -> Iterator[tuple[tuple[bytes, bytes, bytes, bytes],
                                tuple[bytes, bytes, bytes, bytes]]]:
    """Yield validated R1/R2 pairs; enforce name/order/count equality and the
    512 MiB decompressed budget. Any violation aborts the whole check."""
    used = 0
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        try:
            # BufferedReader makes readline a C-speed buffered operation.
            g1 = io.BufferedReader(gzip.GzipFile(fileobj=f1),
                                   buffer_size=1 << 20)
            g2 = io.BufferedReader(gzip.GzipFile(fileobj=f2),
                                   buffer_size=1 << 20)
            try:
                rl1, rl2 = g1.readline, g2.readline
                index = 0
                while True:
                    try:
                        r1 = _read_record(rl1, "R1", index)
                        r2 = _read_record(rl2, "R2", index)
                    except (OSError, EOFError, zlib.error) as exc:
                        raise InspectError(
                            f"gzip stream corrupted: {exc}") from exc
                    if r1 is None and r2 is None:
                        return
                    if r1 is None or r2 is None:
                        raise InspectError(
                            "R1 and R2 contain different numbers of records")
                    n1 = r1[0].split()[0]
                    n2 = r2[0].split()[0]
                    if n1 != n2:
                        raise InspectError(
                            f"record {index}: R1/R2 read names differ "
                            f"({n1.decode('ascii', 'replace')} vs "
                            f"{n2.decode('ascii', 'replace')})")
                    used += (len(r1[0]) + len(r1[1]) + len(r1[2]) + len(r1[3])
                             + len(r2[0]) + len(r2[1]) + len(r2[2])
                             + len(r2[3]) + 16)
                    if used > _SIZE_BUDGET:
                        raise InspectError(
                            f"decompressed total exceeds {_SIZE_BUDGET} bytes")
                    yield r1, r2
                    index += 1
            finally:
                g1.close()
                g2.close()
        except InspectError:
            raise
        except (OSError, EOFError, zlib.error) as exc:
            raise InspectError(f"gzip stream corrupted: {exc}") from exc


# ---------------- exact counting with bounded memory ----------------

class _SpillCounter:
    """Exact index counter that spills raw 16-byte indices to partition files.

    Each index lands in exactly one partition (CRC32 of the index), so
    counting partitions one at a time yields exact global counts while only
    one partition's worth of data is ever held in memory.
    """

    def __init__(self, scratch: Path):
        self._dir = scratch
        self._bufs = [bytearray() for _ in range(_PARTITIONS)]
        self._handles: dict[int, BinaryIO] = {}

    def add(self, idx16: bytes) -> None:
        part = zlib.crc32(idx16) & (_PARTITIONS - 1)
        buf = self._bufs[part]
        buf += idx16
        if len(buf) >= _FLUSH_BYTES:
            self._flush(part)

    def _flush(self, part: int) -> None:
        buf = self._bufs[part]
        if not buf:
            return
        fh = self._handles.get(part)
        if fh is None:
            fh = open(self._dir / f"part-{part:03d}.bin", "wb")
            self._handles[part] = fh
        fh.write(buf)
        self._bufs[part] = bytearray()

    def finish(self) -> None:
        for part in range(_PARTITIONS):
            self._flush(part)
        for fh in self._handles.values():
            fh.close()
        self._handles.clear()

    def top(self, n: int) -> list[tuple[bytes, int]]:
        """Global top-n as [(index, count)] ordered by count desc, index asc.

        A global top-n index is always inside its own partition's top-n, so
        merging per-partition top-n candidates is exact.
        """
        key = lambda kv: (-kv[1], kv[0])  # noqa: E731
        best: list[tuple[bytes, int]] = []
        fmt = f"{_INDEX_LEN}s"
        for path in sorted(self._dir.glob("part-*.bin")):
            with open(path, "rb") as fh:
                data = fh.read()
            if not data:
                continue
            counts = Counter(struct.iter_unpack(fmt, data))
            candidates = [(k[0], v) for k, v in
                          heapq.nsmallest(n, counts.items(), key=key)]
            best = sorted(best + candidates, key=key)[:n]
        return best


# ---------------- public entry point ----------------

def inspect_lane_quarantine(lane: str, r1_path: Path, r2_path: Path,
                            samples: list[dict]) -> dict:
    """Scan the published quarantine pair files and return exact statistics.

    ``samples`` is this lane's barcode table (dicts with sample_id,
    barcode_r1, barcode_r2). Raises MissingDataError when published files are
    absent and InspectError on any data problem; no partial result is ever
    returned. Published files are opened read-only and never modified.
    """
    root = config.PUBLISHED_DIR.resolve()
    paths: list[Path] = []
    for p in (r1_path, r2_path):
        rp = Path(p).resolve()
        if root not in rp.parents:
            raise InspectError("quarantine path escapes the published root")
        if not rp.is_file():
            raise MissingDataError(
                f"published quarantine file is missing: {rp.name}")
        paths.append(rp)

    scratch = config.WORK_DIR / f"inspect-{uuid.uuid4().hex}"
    scratch.mkdir(parents=True)
    lock = open(scratch / "lock", "w")
    try:
        # The flock makes startup recovery on the sibling instance treat this
        # scratch directory as live work and leave it alone.
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock.write(f"{uuid.getnode():x}\n")
        lock.flush()

        counter = _SpillCounter(scratch)
        total = 0
        ranked = 0
        try:
            for r1, r2 in _paired(paths[0], paths[1]):
                total += 1
                s1 = r1[1][:8]
                s2 = r2[1][:8]
                if len(s1) < 8 or len(s2) < 8 or b"N" in s1 or b"N" in s2:
                    continue  # short index or N in first 8: counted separately
                ranked += 1
                counter.add(s1 + s2)
        finally:
            counter.finish()
        top = counter.top(TOP_N)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        finally:
            lock.close()

    # Exact unique barcode lookups per end, for cross-sample flagging.
    r1_map: dict[str, list[str]] = {}
    r2_map: dict[str, list[str]] = {}
    for s in samples:
        r1_map.setdefault(s["barcode_r1"], []).append(s["sample_id"])
        r2_map.setdefault(s["barcode_r2"], []).append(s["sample_id"])

    top_out: list[dict] = []
    flagged: set[str] = set()
    for rank, (idx_bytes, count) in enumerate(top, 1):
        index = idx_bytes.decode()
        idx_r1, idx_r2 = index[:8], index[8:]
        m1 = r1_map.get(idx_r1, [])
        m2 = r2_map.get(idx_r2, [])
        s1 = m1[0] if len(m1) == 1 else None
        s2 = m2[0] if len(m2) == 1 else None
        cross = s1 is not None and s2 is not None and s1 != s2
        if cross:
            flagged.update((s1, s2))
        top_out.append({
            "rank": rank,
            "index": index,
            "index_r1": idx_r1,
            "index_r2": idx_r2,
            "pairs": count,
            "r1_sample": s1,
            "r2_sample": s2,
            "suspected_cross_sample": cross,
        })

    return {
        "lane": lane,
        "total_pairs": total,
        "ranked_pairs": ranked,
        "unranked_pairs": total - ranked,
        "top": top_out,
        "suspected_cross_samples": sorted(flagged),
    }
