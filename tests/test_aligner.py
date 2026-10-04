from affine_match.aligner import Scoring, align

SCORING = Scoring(match=2, mismatch=1, gap_open=3, gap_extend=1)


def rescore(aln, scoring=SCORING):
    """Recompute the score from the alignment strings."""
    total = 0
    gap_len = 0
    for op, q, r in zip(aln.operations, aln.aligned_query, aln.aligned_ref):
        if op in "ID":
            gap_len += 1
            continue
        if gap_len:
            total -= scoring.gap_open + (gap_len - 1) * scoring.gap_extend
            gap_len = 0
        total += scoring.match if q == r else -scoring.mismatch
    if gap_len:
        total -= scoring.gap_open + (gap_len - 1) * scoring.gap_extend
    return total


def test_exact_match():
    aln = align("ACGT", "TTACGTAA", SCORING)
    assert aln.score == 8
    assert (aln.ref_start, aln.ref_end) == (2, 6)
    assert (aln.query_start, aln.query_end) == (0, 4)
    assert aln.operations == "===="
    assert rescore(aln) == aln.score


def test_contiguous_gap_affine_cost():
    # query has a 3-base insertion relative to the reference; the flanking
    # matches make the gapped alignment beat any ungapped local one
    aln = align("ACGTGGGACGT", "ACGTACGT", SCORING)
    assert aln.aligned_query == "ACGTGGGACGT"
    assert aln.aligned_ref == "ACGT---ACGT"
    assert aln.operations == "====III===="
    # 8 matches (16) minus gap of length 3 (3 + 2*1) = 11
    assert aln.score == 11
    assert rescore(aln) == aln.score


def test_gap_in_query():
    aln = align("ACGTACGT", "ACGTGGGACGT", SCORING)
    assert aln.aligned_query == "ACGT---ACGT"
    assert aln.aligned_ref == "ACGTGGGACGT"
    assert aln.operations == "====DDD===="
    assert aln.score == 11
    assert rescore(aln) == aln.score


def test_no_match_returns_zero():
    aln = align("AAAA", "CCCC", SCORING)
    assert aln.score == 0
    assert aln.operations == ""


def test_mismatch_penalized():
    aln = align("ACGT", "ACCT", SCORING)
    assert aln.score == 3 * 2 - 1
    assert "X" in aln.operations
    assert rescore(aln) == aln.score
