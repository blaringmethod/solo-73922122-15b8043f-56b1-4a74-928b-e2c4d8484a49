"""Tests for POST /lanes/{lane}/quarantine/inspect."""
import gzip
import hashlib
import json
import resource
import time
from pathlib import Path

import pytest

from app import processing
from app import quarantine_inspect as qi

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
    {"sample_id": "s2", "barcode_r1": "GGGGGGGG", "barcode_r2": "TTTTTTTT"},
]


def _post(client, lane, barcodes, r1_path, r2_path):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": json.dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )


def _wait(client, job_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(f"/jobs/{job_id}").json()
        if row["status"] in ("published", "failed"):
            return row
        time.sleep(0.1)
    raise AssertionError("timeout waiting for job")


def _rec(name, idx1, idx2, tail1=b"ACGT", tail2=b"TGCA"):
    """One quarantine pair: first 8 bases of each mate form the index."""
    if isinstance(name, str):
        name = name.encode()
    seq1 = idx1 + tail1
    seq2 = idx2 + tail2
    return ([(name, seq1, b"I" * len(seq1))],
            [(name, seq2, b"J" * len(seq2))])


def _publish_directly(data_dir, lane, samples, r1_records, r2_records):
    """Create a published lane without running the worker: real quarantine
    files on disk plus the DB rows a published lane would have."""
    from app import config as cfg, db as dbmod

    lane_dir = cfg.PUBLISHED_DIR / "lanes" / lane
    lane_dir.mkdir(parents=True, exist_ok=True)
    counts = []
    for s in samples:
        sid = s["sample_id"]
        for mate in ("r1", "r2"):
            with gzip.open(lane_dir / f"{sid}.{mate}.fastq.gz", "wb") as g:
                g.write(b"")
        counts.append({"sample_id": sid,
                       "r1_file": str(lane_dir / f"{sid}.r1.fastq.gz"),
                       "r2_file": str(lane_dir / f"{sid}.r2.fastq.gz"),
                       "read_pairs": 0})
    q1 = lane_dir / f"{processing.QUARANTINE}.r1.fastq.gz"
    q2 = lane_dir / f"{processing.QUARANTINE}.r2.fastq.gz"
    with gzip.open(q1, "wb") as g:
        for name, seq, qual in r1_records:
            g.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")
    with gzip.open(q2, "wb") as g:
        for name, seq, qual in r2_records:
            g.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")
    counts.append({"sample_id": processing.QUARANTINE,
                   "r1_file": str(q1), "r2_file": str(q2),
                   "read_pairs": len(r1_records)})
    job_id = f"job-{lane}"
    dbmod.create_job(
        {"job_id": job_id, "lane": lane, "input_hash": "x",
         "barcode_table": json.dumps(samples),
         "r1_incoming": "-", "r2_incoming": "-",
         "created_at": dbmod.utcnow()},
        samples)
    dbmod.publish_job(job_id, lane, counts)
    return q1, q2


def _inspect(client, lane):
    return client.post(f"/lanes/{lane}/quarantine/inspect")


# ---------------- happy path through the real submit pipeline ----------------

def test_inspect_happy_path(client, helpers, tmp_path):
    r1, r2 = [], []
    def add(name, idx1, idx2):
        a, b = _rec(name, idx1, idx2)
        r1.extend(a)
        r2.extend(b)

    for i in range(5):
        add(f"c1_{i}", b"CCCCCCCC", b"AAAAAAAA")   # ranked, count 5
    for i in range(3):
        add(f"c2_{i}", b"TTTTTTTT", b"GGGGGGGG")   # ranked, count 3
    for i in range(2):
        add(f"c3_{i}", b"AAAAAAAA", b"TTTTTTTT")   # cross s1.r1 x s2.r2
    # short R1 (7 bases) -> quarantined by demux, unranked for inspect
    r1.append((b"short1", b"AAAAAAA", b"I" * 7))
    r2.append((b"short1", b"CCCCCCCCTGCA", b"J" * 12))
    # N inside the first 8 bases of R1 -> unranked
    add(b"nread", b"NNNNNNNN", b"ACGTACGT")

    p1, p2 = helpers.make_pair(r1, r2, tmp_path)
    resp = _post(client, "inspLane", BARCODES, p1, p2)
    assert resp.status_code == 201, resp.text
    done = _wait(client, resp.json()["job"]["job_id"])
    assert done["status"] == "published", done
    assert {c["sample_id"]: c["read_pairs"]
            for c in client.get("/lanes/inspLane").json()["counts"]
            }["__quarantine__"] == 12

    resp = _inspect(client, "inspLane")
    assert resp.status_code == 200, resp.text
    info = resp.json()
    assert info["lane"] == "inspLane"
    assert info["total_pairs"] == 12
    assert info["ranked_pairs"] == 10
    assert info["unranked_pairs"] == 2
    assert info["published_pairs"] == 12

    top = info["top"]
    assert [(t["index"], t["pairs"]) for t in top] == [
        ("CCCCCCCC" + "AAAAAAAA", 5),
        ("TTTTTTTT" + "GGGGGGGG", 3),
        ("AAAAAAAA" + "TTTTTTTT", 2),
    ]
    assert [t["rank"] for t in top] == [1, 2, 3]
    cross = top[2]
    assert cross["index_r1"] == "AAAAAAAA"
    assert cross["index_r2"] == "TTTTTTTT"
    assert cross["r1_sample"] == "s1"
    assert cross["r2_sample"] == "s2"
    assert cross["suspected_cross_sample"] is True
    assert top[0]["suspected_cross_sample"] is False
    assert top[0]["r1_sample"] is None and top[0]["r2_sample"] is None
    assert info["suspected_cross_samples"] == ["s1", "s2"]


# ---------------- ranking rules (direct published-lane setup) ----------------

def test_inspect_top20_limit_and_tiebreak(client, data_dir):
    # 25 distinct combos: counts 1..25 plus three tied combos at count 30.
    r1, r2 = [], []
    expected = {}
    bases = [b"CCCCCCCC", b"CCCCCCCA", b"CCCCCCCG", b"CCCCCCCT", b"CCCCCCAA"]
    n = 0
    for i, b1 in enumerate(bases):
        for j, b2 in enumerate(bases):
            count = i * 5 + j + 1  # 1..25, all distinct
            for k in range(count):
                a, b = _rec(f"r{n}", b1, b2)
                r1.extend(a)
                r2.extend(b)
                n += 1
            expected[(b1 + b2).decode()] = count
    # Three combos tied at 30: lexicographic order decides their ranks.
    for b1, b2 in ((b"AAAAACCC", b"GGGGTTTT"), (b"TTTTAAAA", b"CCCCGGGG"),
                   (b"GGGGCCCC", b"AAAATTTT")):
        for k in range(30):
            a, b = _rec(f"t{n}", b1, b2)
            r1.extend(a)
            r2.extend(b)
            n += 1
        expected[(b1 + b2).decode()] = 30

    _publish_directly(data_dir, "topLane", BARCODES, r1, r2)
    resp = _inspect(client, "topLane")
    assert resp.status_code == 200, resp.text
    info = resp.json()
    assert info["total_pairs"] == sum(expected.values())
    assert info["ranked_pairs"] == info["total_pairs"]
    assert info["unranked_pairs"] == 0
    assert len(info["top"]) == 20  # capped at 20 of the 28 combos

    want = sorted(expected.items(), key=lambda kv: (-kv[1], kv[0]))[:20]
    assert [(t["index"], t["pairs"]) for t in info["top"]] == want
    assert [t["rank"] for t in info["top"]] == list(range(1, 21))


def test_inspect_cross_sample_flag_rules(client, data_dir):
    samples = BARCODES + [
        {"sample_id": "s3", "barcode_r1": "AAAAAAAA",
         "barcode_r2": "GGGGGGGG"},  # shares its R1 barcode with s1
    ]
    combos = [
        (b"AAAAAAAA", b"TTTTTTTT"),  # r1 hits s1 AND s3 -> not unique
        (b"GGGGGGGG", b"CCCCCCCC"),  # r1=s2 only, r2=s1 only -> cross
        (b"CCCCCCCC", b"GGGGGGGG"),  # r1 no hit; r2=s3 only -> no flag
        (b"TTTTTTTT", b"AAAAAAAA"),  # no hits at all
        (b"GGGGGGGG", b"TTTTTTTT"),  # exact s2 combo: same sample, no flag
    ]
    r1, r2 = [], []
    for i, (b1, b2) in enumerate(combos):
        a, b = _rec(f"x{i}", b1, b2)
        r1.extend(a)
        r2.extend(b)
    _publish_directly(data_dir, "crossLane", samples, r1, r2)

    info = _inspect(client, "crossLane").json()
    assert info["total_pairs"] == len(combos)
    by_index = {t["index"]: t for t in info["top"]}
    assert by_index["AAAAAAAATTTTTTTT"]["r1_sample"] is None  # s1 and s3
    assert by_index["AAAAAAAATTTTTTTT"]["r2_sample"] == "s2"
    assert by_index["AAAAAAAATTTTTTTT"]["suspected_cross_sample"] is False

    flagged = by_index["GGGGGGGGCCCCCCCC"]
    assert flagged["r1_sample"] == "s2"
    assert flagged["r2_sample"] == "s1"
    assert flagged["suspected_cross_sample"] is True

    assert by_index["CCCCCCCCGGGGGGGG"]["suspected_cross_sample"] is False
    assert by_index["CCCCCCCCGGGGGGGG"]["r2_sample"] == "s3"
    assert by_index["TTTTTTTTAAAAAAAA"]["suspected_cross_sample"] is False
    assert by_index["GGGGGGGGTTTTTTTT"]["suspected_cross_sample"] is False
    assert info["suspected_cross_samples"] == ["s1", "s2"]


def test_inspect_unranked_short_and_n_bases(client, data_dir):
    r1 = [
        (b"a", b"AAAAAAA", b"I" * 7),                    # R1 too short
        (b"b", b"CCCCCCCCACGT", b"I" * 12),
        (b"c", b"NNNNNNNNACGT", b"I" * 12),              # N in first 8 of R1
        (b"d", b"CCCCCCCCACGT", b"I" * 12),
        (b"e", b"CCCCCCCCNACG", b"I" * 12),              # N at pos 9: ranked
    ]
    r2 = [
        (b"a", b"CCCCCCCCTGCA", b"J" * 12),
        (b"b", b"GGGGGG", b"J" * 6),                     # R2 too short
        (b"c", b"CCCCCCCCTGCA", b"J" * 12),
        (b"d", b"CCCCNCCCTGCA", b"J" * 12),              # N in first 8 of R2
        (b"e", b"TTTTTTTTTGCA", b"J" * 12),
    ]
    _publish_directly(data_dir, "nLane", BARCODES, r1, r2)
    info = _inspect(client, "nLane").json()
    assert info["total_pairs"] == 5
    assert info["unranked_pairs"] == 4
    assert info["ranked_pairs"] == 1
    assert [t["index"] for t in info["top"]] == ["CCCCCCCCTTTTTTTT"]


# ---------------- input limits ----------------

def test_inspect_line_limits(client, data_dir):
    def lane_with(name1, seq1, qual1, name2=None, seq2=None, qual2=None,
                  lane="limLane"):
        r1 = [(name1, seq1, qual1)]
        r2 = [(name2 or name1, seq2 or b"CCCCCCCC", qual2 or b"J" * 8)]
        _publish_directly(data_dir, lane, BARCODES, r1, r2)
        return _inspect(client, lane)

    # Boundary values are accepted: 256-byte name, 1000-byte seq/qual.
    ok = lane_with(b"n" * 256, b"CCCCCCCC" + b"A" * 992, b"I" * 1000,
                   seq2=b"CCCCCCCC", qual2=b"J" * 8, lane="okLane")
    assert ok.status_code == 200, ok.text
    assert ok.json()["total_pairs"] == 1

    resp = lane_with(b"a", b"C" * 1001, b"I" * 1001, lane="bigSeq")
    assert resp.status_code == 422
    assert "sequence line exceeds 1000 bytes" in resp.json()["detail"]

    resp = lane_with(b"a", b"C" * 1000, b"I" * 1001, lane="bigQual")
    assert resp.status_code == 422
    assert "quality line exceeds 1000 bytes" in resp.json()["detail"]

    resp = lane_with(b"n" * 257, b"CCCCCCCC", b"I" * 8, lane="bigName")
    assert resp.status_code == 422
    assert "read name exceeds 256 bytes" in resp.json()["detail"]


# ---------------- failure modes: all-or-nothing ----------------

def test_inspect_missing_quarantine_file(client, data_dir):
    q1, q2 = _publish_directly(
        data_dir, "missLane", BARCODES,
        [(b"a", b"CCCCCCCCACGT", b"I" * 12)],
        [(b"a", b"AAAAAAAATGCA", b"J" * 12)])
    q2.unlink()
    resp = _inspect(client, "missLane")
    assert resp.status_code == 409
    assert "missing" in resp.json()["detail"].lower()


def test_inspect_corrupt_gzip(client, data_dir):
    q1, q2 = _publish_directly(
        data_dir, "gzLane", BARCODES,
        [(b"a", b"CCCCCCCCACGT", b"I" * 12)],
        [(b"a", b"AAAAAAAATGCA", b"J" * 12)])
    q1.write_bytes(q1.read_bytes()[:-6])  # truncate the gzip stream
    resp = _inspect(client, "gzLane")
    assert resp.status_code == 422
    assert "gzip" in resp.json()["detail"].lower()


def test_inspect_unpaired_mates(client, data_dir):
    # Different record counts.
    _publish_directly(
        data_dir, "cntLane", BARCODES,
        [(b"a", b"CCCCCCCCACGT", b"I" * 12),
         (b"b", b"CCCCCCCCACGT", b"I" * 12)],
        [(b"a", b"AAAAAAAATGCA", b"J" * 12)])
    resp = _inspect(client, "cntLane")
    assert resp.status_code == 422
    assert "different numbers" in resp.json()["detail"]

    # Same count, mismatched read names.
    _publish_directly(
        data_dir, "nameLane", BARCODES,
        [(b"a", b"CCCCCCCCACGT", b"I" * 12)],
        [(b"b", b"AAAAAAAATGCA", b"J" * 12)])
    resp = _inspect(client, "nameLane")
    assert resp.status_code == 422
    assert "read names differ" in resp.json()["detail"]


def test_inspect_malformed_fastq(client, data_dir):
    _publish_directly(
        data_dir, "fmtLane", BARCODES,
        [(b"a", b"CCCCCCCCACGT", b"I" * 11)],   # seq/qual length mismatch
        [(b"a", b"AAAAAAAATGCA", b"J" * 12)])
    resp = _inspect(client, "fmtLane")
    assert resp.status_code == 422
    assert "mismatch" in resp.json()["detail"]


def test_inspect_unknown_or_unpublished_lane(client, helpers, tmp_path):
    assert _inspect(client, "ghostLane").status_code == 404

    # A submitted-but-failed lane is not published: inspect is a 404.
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    r1.write_bytes(b"not a gzip")
    helpers.write_fastq_gz(r2, [(b"a", b"CCCCCCCC", b"I" * 8)])
    resp = _post(client, "failedLane", BARCODES, r1, r2)
    job = _wait(client, resp.json()["job"]["job_id"])
    assert job["status"] == "failed"
    assert _inspect(client, "failedLane").status_code == 404


def test_inspect_is_read_only_and_cleans_scratch(client, data_dir):
    from app import config as cfg
    q1, q2 = _publish_directly(
        data_dir, "roLane", BARCODES,
        [(b"a", b"CCCCCCCCACGT", b"I" * 12)],
        [(b"a", b"AAAAAAAATGCA", b"J" * 12)])
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    before = (digest(q1), digest(q2), q1.stat().st_mtime_ns,
              q2.stat().st_mtime_ns)

    resp = _inspect(client, "roLane")
    assert resp.status_code == 200
    after = (digest(q1), digest(q2), q1.stat().st_mtime_ns,
             q2.stat().st_mtime_ns)
    assert before == after  # published files untouched
    assert list(cfg.WORK_DIR.glob("inspect-*")) == []  # scratch removed

    # Same on the failure path: corrupt the file, check still cleans up.
    q1.write_bytes(q1.read_bytes()[:-6])
    resp = _inspect(client, "roLane")
    assert resp.status_code == 422
    assert list(cfg.WORK_DIR.glob("inspect-*")) == []


# ---------------- exact counting under spilling ----------------

def test_inspect_exact_counts_with_spill(client, data_dir, monkeypatch):
    # Tiny buffers/partitions force real spill+merge behaviour.
    monkeypatch.setattr(qi, "_FLUSH_BYTES", 16 * 16)  # 16 indices per flush
    monkeypatch.setattr(qi, "_PARTITIONS", 8)

    def mer(i):
        return "".join("ACGT"[(i >> (2 * k)) & 3] for k in range(8)).encode()

    pool = [mer(i) for i in range(64)]  # 64 distinct 8-mers

    expected = {}
    r1, r2 = [], []
    n = 0
    for i, b1 in enumerate(pool):
        for j in range(3):
            b2 = pool[(i * 5 + j * 13) % len(pool)]
            count = (i * 37 + j * 11) % 23 + 1
            key = (b1 + b2).decode()
            expected[key] = expected.get(key, 0) + count
            for _ in range(count):
                a, b = _rec(f"s{n}", b1, b2)
                r1.extend(a)
                r2.extend(b)
                n += 1

    _publish_directly(data_dir, "spillLane", BARCODES, r1, r2)
    resp = _inspect(client, "spillLane")
    assert resp.status_code == 200, resp.text
    info = resp.json()
    assert info["total_pairs"] == sum(expected.values())
    want = sorted(expected.items(), key=lambda kv: (-kv[1], kv[0]))[:20]
    assert [(t["index"], t["pairs"]) for t in info["top"]] == want


def test_inspect_streaming_memory(client, data_dir):
    """Distinct-index-heavy quarantine stays in small RSS (spill counting)."""
    pairs = 300_000
    q1 = data_dir / "big.r1.fastq.gz"
    q2 = data_dir / "big.r2.fastq.gz"
    digits = "ACGT"

    def idx16(i):
        return "".join(digits[(i >> s) & 3] for s in range(0, 32, 2))

    with gzip.open(q1, "wb") as g1, gzip.open(q2, "wb") as g2:
        for i in range(pairs):
            idx = idx16(i).encode()
            g1.write(b"@r%d\n%sACGT\n+\n%s\n" % (i, idx, b"I" * 20))
            g2.write(b"@r%d\n%sTGCA\n+\n%s\n" % (i, idx[8:] + idx[:8],
                                                 b"J" * 20))
    # Note: R2 first-8 is the second half of idx16(i); every pair is a
    # distinct combination by construction of R1 alone.

    from app import config as cfg
    lane_dir = cfg.PUBLISHED_DIR / "lanes" / "bigLane"
    lane_dir.mkdir(parents=True, exist_ok=True)
    q1.replace(lane_dir / f"{processing.QUARANTINE}.r1.fastq.gz")
    q2.replace(lane_dir / f"{processing.QUARANTINE}.r2.fastq.gz")
    from app import db as dbmod
    dbmod.create_job(
        {"job_id": "job-bigLane", "lane": "bigLane", "input_hash": "x",
         "barcode_table": json.dumps(BARCODES),
         "r1_incoming": "-", "r2_incoming": "-",
         "created_at": dbmod.utcnow()}, BARCODES)
    dbmod.publish_job("job-bigLane", "bigLane", [
        {"sample_id": processing.QUARANTINE,
         "r1_file": str(lane_dir / f"{processing.QUARANTINE}.r1.fastq.gz"),
         "r2_file": str(lane_dir / f"{processing.QUARANTINE}.r2.fastq.gz"),
         "read_pairs": pairs}])

    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    resp = _inspect(client, "bigLane")
    assert resp.status_code == 200, resp.text
    info = resp.json()
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    growth_mib = (rss_after - rss_before) / 1024
    print(f"\ninspected {pairs} distinct pairs, RSS growth {growth_mib:.1f} MiB")
    assert growth_mib < 64, f"RSS grew {growth_mib:.1f} MiB during inspect"
    assert info["total_pairs"] == pairs
    assert info["ranked_pairs"] == pairs
    # All counts are 1, so the top 20 must be the 20 smallest indices.
    want = sorted((idx16(i) for i in range(pairs)))[:20]
    assert [t["index"] for t in info["top"]] == want
    assert all(t["pairs"] == 1 for t in info["top"])
