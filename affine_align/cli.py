"""Command-line entry point for batch local alignment."""

from __future__ import annotations

import argparse
import sys

from .align import Scoring
from .fasta import parse_fasta
from .matcher import match_all
from .report import write_csv, write_json


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="affine-align",
        description=(
            "Batch Smith-Waterman local alignment of query DNA fragments "
            "against reference sequences with affine gap penalties."
        ),
    )
    parser.add_argument("--reference", required=True, help="reference FASTA file")
    parser.add_argument("--query", required=True, help="query FASTA file")
    parser.add_argument("--json", dest="json_out", required=True, help="JSON summary output")
    parser.add_argument("--csv", dest="csv_out", required=True, help="CSV assignment output")
    parser.add_argument("--match", type=int, default=2, help="match score (positive)")
    parser.add_argument("--mismatch", type=int, default=1, help="mismatch penalty (positive)")
    parser.add_argument("--gap-open", type=int, default=3, help="gap open penalty (positive)")
    parser.add_argument("--gap-extend", type=int, default=1, help="gap extend penalty (positive)")
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
        references = parse_fasta(args.reference)
        queries = parse_fasta(args.query)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    results = match_all(queries, references, scoring)
    write_json(results, scoring, args.json_out)
    write_csv(results, args.csv_out)
    matched = sum(1 for r in results if r.matched)
    print(f"{matched}/{len(results)} queries matched; wrote {args.json_out} and {args.csv_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
