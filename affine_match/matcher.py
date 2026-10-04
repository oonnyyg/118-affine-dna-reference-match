"""Batch matching of query reads against reference sequences."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from .aligner import Alignment, Scoring, align
from .fasta import FastaRecord

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class MatchResult:
    query_id: str
    query_description: str
    query_length: int
    matched: bool
    reference_id: str = ""
    reference_description: str = ""
    strand: str = ""  # '+' or '-'
    score: int = 0
    ref_start: int = 0
    ref_end: int = 0
    query_start: int = 0  # always in original query coordinates
    query_end: int = 0
    coverage: float = 0.0
    aligned_ref: str = ""
    aligned_query: str = ""
    operations: str = ""


def _pick_best(
    candidates: list[tuple[FastaRecord, str, Alignment]],
) -> tuple[FastaRecord, str, Alignment]:
    """Highest score; ties: reference id lexicographic, then '+' strand."""
    return max(
        candidates,
        key=lambda c: (
            c[2].score,
            tuple(-ord(ch) for ch in c[0].ident),
            1 if c[1] == "+" else 0,
        ),
    )


def match_queries(
    queries: list[FastaRecord],
    references: list[FastaRecord],
    scoring: Scoring,
) -> list[MatchResult]:
    """Align every query (both strands) to every reference, keep the best."""
    results: list[MatchResult] = []
    for query in queries:
        rc = reverse_complement(query.sequence)
        candidates: list[tuple[FastaRecord, str, Alignment]] = []
        for ref in references:
            candidates.append(
                (ref, "+", align(query.sequence, ref.sequence, scoring))
            )
            candidates.append((ref, "-", align(rc, ref.sequence, scoring)))
        ref, strand, best = _pick_best(candidates)
        if best.score == 0:
            results.append(
                MatchResult(
                    query_id=query.ident,
                    query_description=query.description,
                    query_length=len(query.sequence),
                    matched=False,
                )
            )
            continue
        if strand == "+":
            q_start, q_end = best.query_start, best.query_end
            aligned_query = best.aligned_query
        else:
            # map reverse-complement coordinates back to the original query
            qlen = len(query.sequence)
            q_start, q_end = qlen - best.query_end, qlen - best.query_start
            aligned_query = reverse_complement(best.aligned_query)
        results.append(
            MatchResult(
                query_id=query.ident,
                query_description=query.description,
                query_length=len(query.sequence),
                matched=True,
                reference_id=ref.ident,
                reference_description=ref.description,
                strand=strand,
                score=best.score,
                ref_start=best.ref_start,
                ref_end=best.ref_end,
                query_start=q_start,
                query_end=q_end,
                coverage=(q_end - q_start) / len(query.sequence),
                aligned_ref=best.aligned_ref,
                aligned_query=aligned_query,
                operations=best.operations,
            )
        )
    return results


def export_json(results: list[MatchResult], path: str | Path) -> None:
    payload = [
        {
            "query_id": r.query_id,
            "query_description": r.query_description,
            "query_length": r.query_length,
            "matched": r.matched,
            "reference_id": r.reference_id or None,
            "reference_description": r.reference_description or None,
            "strand": r.strand or None,
            "score": r.score,
            "ref_interval": [r.ref_start, r.ref_end] if r.matched else None,
            "query_interval": [r.query_start, r.query_end] if r.matched else None,
            "coverage": round(r.coverage, 6),
            "aligned_ref": r.aligned_ref or None,
            "aligned_query": r.aligned_query or None,
            "operations": r.operations or None,
        }
        for r in results
    ]
    Path(path).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def export_csv(results: list[MatchResult], path: str | Path) -> None:
    fields = [
        "query_id",
        "matched",
        "reference_id",
        "strand",
        "score",
        "ref_start",
        "ref_end",
        "query_start",
        "query_end",
        "coverage",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in results:
            writer.writerow(
                [
                    r.query_id,
                    r.matched,
                    r.reference_id,
                    r.strand,
                    r.score,
                    r.ref_start if r.matched else "",
                    r.ref_end if r.matched else "",
                    r.query_start if r.matched else "",
                    r.query_end if r.matched else "",
                    f"{r.coverage:.6f}" if r.matched else "",
                ]
            )
