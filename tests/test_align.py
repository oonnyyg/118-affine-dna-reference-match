from affine_align.align import Scoring, reverse_complement, smith_waterman

SCORING = Scoring(match=2, mismatch=1, gap_open=3, gap_extend=1)


def check(alignment):
    assert alignment.recomputed_score(SCORING) == alignment.score
    assert len(alignment.aligned_query) == len(alignment.aligned_ref)
    assert len(alignment.operations) == len(alignment.aligned_query)


def test_exact_match():
    aln = smith_waterman("ACGT", "TTACGTAA", SCORING)
    assert aln.score == 8
    assert (aln.query_start, aln.query_end) == (0, 4)
    assert (aln.ref_start, aln.ref_end) == (2, 6)
    assert aln.operations == "===="
    check(aln)


def test_contiguous_gap_single_penalty():
    # Reference has a 3-base insertion relative to the query: one gap of
    # length 3 must cost gap_open + 2 * gap_extend, not three opens.
    aln = smith_waterman("ACGTACGT", "ACGTGGGACGT", SCORING)
    assert aln.score == 8 * 2 - (3 + 2 * 1)
    assert "D" in aln.operations
    assert aln.operations.count("D") == 3
    check(aln)


def test_gap_in_query():
    aln = smith_waterman("ACGTGGGACGT", "ACGTACGT", SCORING)
    assert aln.score == 16 - 5
    assert aln.operations.count("I") == 3
    check(aln)


def test_all_mismatch_scores_zero():
    aln = smith_waterman("AAAA", "CCCC", SCORING)
    assert aln.score == 0
    assert aln.operations == ""


def test_local_alignment_clamps_at_zero():
    aln = smith_waterman("CCCCACGTCCCC", "ACGT", SCORING)
    assert aln.score == 8
    assert aln.aligned_query == "ACGT"
    check(aln)


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGT") == "ACTT"
