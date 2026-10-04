"""Smith-Waterman local alignment with affine gap penalties.

Scoring: match/mismatch via `match`/`mismatch` (mismatch is subtracted),
a gap of length k costs `gap_open + (k - 1) * gap_extend`. All scores are
clamped at zero (local alignment).

Tie-breaking (fixed, deterministic):
- When several predecessor states give the same score, the order of
  preference is M (diagonal), then Ix (gap in the reference), then Iy
  (gap in the query).
- When opening a gap, the source state preference is M, then the other
  gap state, then continuing the same gap state is preferred over opening
  a new one when scores tie (extension checked first).
- The best-scoring cell is the first cell (scanning rows top-to-bottom,
  left-to-right) attaining the maximum score.
"""

from __future__ import annotations

from dataclasses import dataclass

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Scoring:
    match: int = 2
    mismatch: int = 1
    gap_open: int = 3
    gap_extend: int = 1

    def __post_init__(self) -> None:
        for name in ("match", "mismatch", "gap_open", "gap_extend"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")


@dataclass(frozen=True)
class Alignment:
    score: int
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    aligned_query: str
    aligned_ref: str
    operations: str  # '=', 'X', 'I' (insertion vs ref), 'D' (deletion vs ref)

    def recomputed_score(self, scoring: Scoring) -> int:
        total = 0
        index = 0
        while index < len(self.operations):
            op = self.operations[index]
            if op in "=X":
                total += scoring.match if op == "=" else -scoring.mismatch
                index += 1
            else:
                length = 0
                while index < len(self.operations) and self.operations[index] == op:
                    length += 1
                    index += 1
                total -= scoring.gap_open + (length - 1) * scoring.gap_extend
        return total


# States
M, IX, IY = 0, 1, 2


def smith_waterman(query: str, ref: str, scoring: Scoring) -> Alignment:
    """Align `query` to `ref` locally; returns the best alignment."""
    n, m = len(query), len(ref)
    neg = float("-inf")
    # Score matrices per state; H is the clamped best score per cell.
    m_score = [[0] * (m + 1) for _ in range(n + 1)]
    ix_score = [[0] * (m + 1) for _ in range(n + 1)]
    iy_score = [[0] * (m + 1) for _ in range(n + 1)]
    # Traceback pointers: for each state, (prev_state, prev_i, prev_j).
    trace_m = [[None] * (m + 1) for _ in range(n + 1)]
    trace_ix = [[None] * (m + 1) for _ in range(n + 1)]
    trace_iy = [[None] * (m + 1) for _ in range(n + 1)]

    best_score = 0
    best_cell = (0, 0)

    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            # Ix: gap in reference, consume query (vertical move).
            cand_ix = [
                (ix_score[i - 1][j] - scoring.gap_extend, IX, 0),
                (m_score[i - 1][j] - scoring.gap_open, M, 1),
                (iy_score[i - 1][j] - scoring.gap_open, IY, 2),
            ]
            value, state, _ = max(cand_ix, key=lambda c: (c[0], -c[2]))
            ix_score[i][j] = value
            trace_ix[i][j] = (state, i - 1, j)
            # Iy: gap in query, consume reference (horizontal move).
            cand_iy = [
                (iy_score[i][j - 1] - scoring.gap_extend, IY, 0),
                (m_score[i][j - 1] - scoring.gap_open, M, 1),
                (ix_score[i][j - 1] - scoring.gap_open, IX, 2),
            ]
            value, state, _ = max(cand_iy, key=lambda c: (c[0], -c[2]))
            iy_score[i][j] = value
            trace_iy[i][j] = (state, i, j - 1)
            # M: diagonal.
            sub = scoring.match if qi == ref[j - 1] else -scoring.mismatch
            cand_m = [
                (m_score[i - 1][j - 1], M, 0),
                (ix_score[i - 1][j - 1], IX, 1),
                (iy_score[i - 1][j - 1], IY, 2),
            ]
            value, state, _ = max(cand_m, key=lambda c: (c[0], -c[2]))
            m_score[i][j] = value + sub
            trace_m[i][j] = (state, i - 1, j - 1)

            cell_best = max(m_score[i][j], ix_score[i][j], iy_score[i][j], 0)
            if cell_best > best_score:
                best_score = cell_best
                best_cell = (i, j)

    return _traceback(
        query, ref, best_score, best_cell, m_score, ix_score, iy_score,
        trace_m, trace_ix, trace_iy,
    )


def _traceback(query, ref, best_score, best_cell, m_score, ix_score, iy_score,
               trace_m, trace_ix, trace_iy):
    if best_score == 0:
        return Alignment(0, 0, 0, 0, 0, "", "", "")
    i, j = best_cell
    # Pick starting state with the fixed preference M, Ix, Iy.
    scores = [(m_score[i][j], M), (ix_score[i][j], IX), (iy_score[i][j], IY)]
    value, state = max(scores, key=lambda c: (c[0], -c[1]))
    aligned_q: list[str] = []
    aligned_r: list[str] = []
    ops: list[str] = []
    while True:
        if state == M:
            pointer = trace_m[i][j]
            aligned_q.append(query[i - 1])
            aligned_r.append(ref[j - 1])
            ops.append("=" if query[i - 1] == ref[j - 1] else "X")
        elif state == IX:
            pointer = trace_ix[i][j]
            aligned_q.append(query[i - 1])
            aligned_r.append("-")
            ops.append("I")
        else:
            pointer = trace_iy[i][j]
            aligned_q.append("-")
            aligned_r.append(ref[j - 1])
            ops.append("D")
        prev_state, i, j = pointer
        prev_scores = (m_score[i][j], ix_score[i][j], iy_score[i][j])
        if prev_scores[prev_state] <= 0 or (i == 0 and j == 0):
            break
        state = prev_state
    aligned_q.reverse()
    aligned_r.reverse()
    ops.reverse()
    return Alignment(
        score=best_score,
        query_start=i,
        query_end=best_cell[0],
        ref_start=j,
        ref_end=best_cell[1],
        aligned_query="".join(aligned_q),
        aligned_ref="".join(aligned_r),
        operations="".join(ops),
    )
