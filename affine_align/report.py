"""JSON and CSV export of match results, preserving query input order."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .align import Scoring
from .matcher import MatchResult


def result_to_dict(result: MatchResult) -> dict:
    entry: dict = {
        "query_id": result.query.identifier,
        "query_description": result.query.description,
        "query_length": len(result.query.sequence),
        "matched": result.matched,
    }
    if not result.matched:
        return entry
    alignment = result.alignment
    assert alignment is not None and result.reference is not None
    entry.update(
        {
            "reference_id": result.reference.identifier,
            "reference_description": result.reference.description,
            "strand": result.strand,
            "score": alignment.score,
            "query_start": result.query_start,
            "query_end": result.query_end,
            "ref_start": result.ref_start,
            "ref_end": result.ref_end,
            "coverage": round(result.coverage, 6),
            "aligned_query": alignment.aligned_query,
            "aligned_ref": alignment.aligned_ref,
            "operations": alignment.operations,
        }
    )
    return entry


def write_json(
    results: list[MatchResult], scoring: Scoring, path: str | Path
) -> None:
    payload = {
        "scoring": {
            "match": scoring.match,
            "mismatch": scoring.mismatch,
            "gap_open": scoring.gap_open,
            "gap_extend": scoring.gap_extend,
        },
        "coordinate_system": "0-based half-open [start, end)",
        "results": [result_to_dict(r) for r in results],
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


CSV_FIELDS = [
    "query_id",
    "matched",
    "reference_id",
    "strand",
    "score",
    "query_start",
    "query_end",
    "ref_start",
    "ref_end",
    "coverage",
]


def write_csv(results: list[MatchResult], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for result in results:
            row = {field: "" for field in CSV_FIELDS}
            row["query_id"] = result.query.identifier
            row["matched"] = "yes" if result.matched else "no"
            if result.matched:
                assert result.reference is not None and result.alignment is not None
                row.update(
                    {
                        "reference_id": result.reference.identifier,
                        "strand": result.strand,
                        "score": result.alignment.score,
                        "query_start": result.query_start,
                        "query_end": result.query_end,
                        "ref_start": result.ref_start,
                        "ref_end": result.ref_end,
                        "coverage": f"{result.coverage:.6f}",
                    }
                )
            writer.writerow(row)
