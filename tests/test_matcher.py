import json

from affine_match.aligner import Scoring
from affine_match.fasta import FastaRecord
from affine_match.matcher import (
    export_csv,
    export_json,
    match_queries,
    reverse_complement,
)

SCORING = Scoring()


def rec(ident, seq):
    return FastaRecord(ident, ident, seq)


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACC") == "GGTT"


def test_reverse_strand_hit_coordinates():
    # reference contains GGTT; query AACC only matches on its reverse strand
    refs = [rec("ref1", "TTTGGTTAAA")]
    queries = [rec("q1", "AACC")]
    (res,) = match_queries(queries, refs, SCORING)
    assert res.matched
    assert res.strand == "-"
    assert res.score == 8
    assert (res.ref_start, res.ref_end) == (3, 7)
    assert (res.query_start, res.query_end) == (0, 4)
    assert res.coverage == 1.0


def test_no_match_reported_without_interval():
    refs = [rec("ref1", "CCCC")]
    queries = [rec("q1", "AAAA")]
    (res,) = match_queries(queries, refs, SCORING)
    assert not res.matched
    assert res.score == 0
    assert res.reference_id == ""


def test_tie_breaks_by_reference_id_then_strand():
    refs = [rec("refB", "ACGT"), rec("refA", "ACGT")]
    queries = [rec("q1", "ACGT")]
    (res,) = match_queries(queries, refs, SCORING)
    assert res.reference_id == "refA"
    assert res.strand == "+"


def test_query_order_preserved_and_exports(tmp_path):
    refs = [rec("ref1", "AAAAAAAA")]
    queries = [rec("q2", "AACC"), rec("q1", "CCCC")]
    results = match_queries(queries, refs, SCORING)
    assert [r.query_id for r in results] == ["q2", "q1"]

    json_path = tmp_path / "out.json"
    csv_path = tmp_path / "out.csv"
    export_json(results, json_path)
    export_csv(results, csv_path)
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert [p["query_id"] for p in payload] == ["q2", "q1"]
    assert payload[0]["matched"] is True
    assert payload[1]["matched"] is False
    assert payload[1]["ref_interval"] is None
    lines = csv_path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 3
    assert lines[1].startswith("q2,True,ref1")
