import json

from affine_align.align import Scoring, reverse_complement
from affine_align.fasta import FastaRecord
from affine_align.matcher import best_match, match_all
from affine_align.report import write_csv, write_json

SCORING = Scoring(match=2, mismatch=1, gap_open=3, gap_extend=1)


def rec(identifier, sequence, description=""):
    return FastaRecord(identifier, description, sequence)


def test_reverse_strand_hit_coordinates():
    # Query only matches the reference on its reverse complement.
    query = rec("q1", "TTACGGAA")
    ref = rec("r1", "GG" + reverse_complement("TTACGGAA") + "GG")
    result = best_match(query, [ref], SCORING)
    assert result.matched
    assert result.strand == "-"
    assert result.alignment.score == 16
    # Full-length reverse hit maps back to the whole original query.
    assert (result.query_start, result.query_end) == (0, 8)
    assert (result.ref_start, result.ref_end) == (2, 10)
    assert result.coverage == 1.0


def test_no_match_reported_for_zero_score():
    result = best_match(rec("q", "AAAA"), [rec("r", "CCCC")], SCORING)
    assert not result.matched
    assert result.reference is None
    assert result.query_start is None
    assert result.coverage == 0.0


def test_tie_breaks_by_reference_id_then_strand():
    query = rec("q", "ACGTACGT")
    ref_b = rec("b", "ACGTACGT")
    ref_a = rec("a", "ACGTACGT")
    result = best_match(query, [ref_b, ref_a], SCORING)
    assert result.reference.identifier == "a"
    assert result.strand == "+"


def test_match_all_preserves_query_order(tmp_path):
    queries = [rec("q2", "AAAA"), rec("q1", "CCCC")]
    refs = [rec("r", "AAAA")]
    results = match_all(queries, refs, SCORING)
    assert [r.query.identifier for r in results] == ["q2", "q1"]

    json_path = tmp_path / "out.json"
    csv_path = tmp_path / "out.csv"
    write_json(results, SCORING, json_path)
    write_csv(results, csv_path)

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert [e["query_id"] for e in payload["results"]] == ["q2", "q1"]
    assert payload["results"][0]["matched"] is True
    assert payload["results"][1]["matched"] is False
    assert "query_start" not in payload["results"][1]

    lines = csv_path.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].startswith("query_id,")
    assert lines[1].startswith("q2,yes,r,+,")
    assert lines[2].startswith("q1,no,,,,")
