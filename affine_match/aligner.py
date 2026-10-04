"""Smith-Waterman local alignment with affine gap penalties.

A length-k gap costs ``gap_open + (k - 1) * gap_extend``. Match and
mismatch values are positive; scores are floored at zero (local).

Fixed tie-breaking rules (deterministic):
  * H cell: diagonal pair > gap in reference (Ix) > gap in query (Iy) > 0;
  * Ix / Iy cell: opening a new gap from H > extending an existing gap;
  * best-cell scan: earliest query end, then earliest reference end.
"""

from __future__ import annotations

from dataclasses import dataclass


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
    """Local alignment of query against reference."""

    score: int
    ref_start: int  # 0-based, half-open
    ref_end: int
    query_start: int
    query_end: int
    aligned_ref: str
    aligned_query: str
    operations: str  # '=' match, 'X' mismatch, 'I' gap in ref, 'D' gap in query


_STOP, _DIAG, _IX, _IY = 0, 1, 2, 3
_OPEN, _EXT = 0, 1


def align(query: str, reference: str, scoring: Scoring) -> Alignment:
    """Return the best local alignment of query against reference."""
    n, m = len(query), len(reference)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    Ix = [[0] * (m + 1) for _ in range(n + 1)]  # gap in reference
    Iy = [[0] * (m + 1) for _ in range(n + 1)]  # gap in query
    tr_h = [[_STOP] * (m + 1) for _ in range(n + 1)]
    tr_x = [[_OPEN] * (m + 1) for _ in range(n + 1)]
    tr_y = [[_OPEN] * (m + 1) for _ in range(n + 1)]

    best_score, best_i, best_j = 0, 0, 0
    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            open_x = H[i - 1][j] - scoring.gap_open
            ext_x = Ix[i - 1][j] - scoring.gap_extend
            if open_x >= ext_x:  # tie: prefer opening
                Ix[i][j], tr_x[i][j] = open_x, _OPEN
            else:
                Ix[i][j], tr_x[i][j] = ext_x, _EXT
            open_y = H[i][j - 1] - scoring.gap_open
            ext_y = Iy[i][j - 1] - scoring.gap_extend
            if open_y >= ext_y:
                Iy[i][j], tr_y[i][j] = open_y, _OPEN
            else:
                Iy[i][j], tr_y[i][j] = ext_y, _EXT
            diag = H[i - 1][j - 1] + (
                scoring.match if qi == reference[j - 1] else -scoring.mismatch
            )
            val, ptr = 0, _STOP
            if Iy[i][j] > val:
                val, ptr = Iy[i][j], _IY
            if Ix[i][j] > val:
                val, ptr = Ix[i][j], _IX
            if diag > val:
                val, ptr = diag, _DIAG
            H[i][j], tr_h[i][j] = val, ptr
            if val > best_score:
                best_score, best_i, best_j = val, i, j

    if best_score == 0:
        return Alignment(0, 0, 0, 0, 0, "", "", "")

    i, j = best_i, best_j
    query_end, ref_end = i, j
    aligned_q: list[str] = []
    aligned_r: list[str] = []
    ops: list[str] = []
    ptr = tr_h[i][j]
    while ptr != _STOP:
        if ptr == _DIAG:
            aligned_q.append(query[i - 1])
            aligned_r.append(reference[j - 1])
            ops.append("=" if query[i - 1] == reference[j - 1] else "X")
            i -= 1
            j -= 1
            ptr = tr_h[i][j]
        elif ptr == _IX:
            aligned_q.append(query[i - 1])
            aligned_r.append("-")
            ops.append("I")
            prev = tr_x[i][j]
            i -= 1
            ptr = _IX if prev == _EXT else tr_h[i][j]
        else:  # _IY
            aligned_q.append("-")
            aligned_r.append(reference[j - 1])
            ops.append("D")
            prev = tr_y[i][j]
            j -= 1
            ptr = _IY if prev == _EXT else tr_h[i][j]

    aligned_q.reverse()
    aligned_r.reverse()
    ops.reverse()
    return Alignment(
        score=best_score,
        ref_start=j,
        ref_end=ref_end,
        query_start=i,
        query_end=query_end,
        aligned_ref="".join(aligned_r),
        aligned_query="".join(aligned_q),
        operations="".join(ops),
    )
