"""Best-reference selection and report building."""

from __future__ import annotations

from dataclasses import dataclass

from .align import Alignment, Scoring, reverse_complement, smith_waterman
from .fasta import FastaRecord


@dataclass(frozen=True)
class MatchResult:
    query: FastaRecord
    reference: FastaRecord | None
    strand: str | None  # '+' or '-'
    alignment: Alignment | None  # coordinates on the oriented query
    query_start: int | None  # 0-based half-open, original query coordinates
    query_end: int | None
    ref_start: int | None
    ref_end: int | None
    coverage: float

    @property
    def matched(self) -> bool:
        return self.alignment is not None and self.alignment.score > 0


def best_match(
    query: FastaRecord, references: list[FastaRecord], scoring: Scoring
) -> MatchResult:
    """Find the highest-scoring reference for one query.

    Both the forward and reverse-complement orientations of the query are
    aligned against every reference. Ties are broken by reference identifier
    (lexicographic), then by strand (forward first).
    """
    best: tuple[int, str, str, Alignment] | None = None
    for ref in references:
        for strand, oriented in (
            ("+", query.sequence),
            ("-", reverse_complement(query.sequence)),
        ):
            alignment = smith_waterman(oriented, ref.sequence, scoring)
            if best is None or alignment.score > best[0]:
                best = (alignment.score, ref.identifier, strand, alignment)
            elif (
                alignment.score == best[0]
                and (ref.identifier, strand) < (best[1], best[2])
            ):
                best = (alignment.score, ref.identifier, strand, alignment)
    assert best is not None
    score, ref_id, strand, alignment = best
    ref = next(r for r in references if r.identifier == ref_id)
    if score == 0:
        return MatchResult(query, None, None, None, None, None, None, None, 0.0)
    if strand == "+":
        q_start, q_end = alignment.query_start, alignment.query_end
    else:
        q_start = len(query.sequence) - alignment.query_end
        q_end = len(query.sequence) - alignment.query_start
    aligned_bases = sum(1 for op in alignment.operations if op in "=XI")
    coverage = aligned_bases / len(query.sequence)
    return MatchResult(
        query=query,
        reference=ref,
        strand=strand,
        alignment=alignment,
        query_start=q_start,
        query_end=q_end,
        ref_start=alignment.ref_start,
        ref_end=alignment.ref_end,
        coverage=coverage,
    )


def match_all(
    queries: list[FastaRecord], references: list[FastaRecord], scoring: Scoring
) -> list[MatchResult]:
    return [best_match(query, references, scoring) for query in queries]
