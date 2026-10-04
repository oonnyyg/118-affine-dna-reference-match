"""Command line interface for the batch local alignment tool."""

from __future__ import annotations

import argparse
import sys

from .aligner import Scoring
from .fasta import FastaError, read_fasta
from .matcher import export_csv, export_json, match_queries


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="affine-match",
        description=(
            "Batch Smith-Waterman local alignment of short DNA queries "
            "against reference sequences (affine gaps, both strands)."
        ),
    )
    parser.add_argument("--reference", required=True, help="reference FASTA")
    parser.add_argument("--query", required=True, help="query FASTA")
    parser.add_argument("--json", required=True, help="JSON summary output path")
    parser.add_argument("--csv", required=True, help="CSV assignment output path")
    parser.add_argument("--match", type=int, default=2, help="match score (>0)")
    parser.add_argument(
        "--mismatch", type=int, default=1, help="mismatch penalty (>0)"
    )
    parser.add_argument(
        "--gap-open", type=int, default=3, help="gap open penalty (>0)"
    )
    parser.add_argument(
        "--gap-extend", type=int, default=1, help="gap extend penalty (>0)"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        scoring = Scoring(
            match=args.match,
            mismatch=args.mismatch,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
        )
        references = read_fasta(args.reference)
        queries = read_fasta(args.query)
    except (FastaError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    results = match_queries(queries, references, scoring)
    export_json(results, args.json)
    export_csv(results, args.csv)
    matched = sum(r.matched for r in results)
    print(f"{matched}/{len(results)} queries matched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
