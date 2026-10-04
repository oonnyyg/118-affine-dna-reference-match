# affine-dna-reference-match

Batch local alignment of short DNA query fragments against reference
sequences, implemented with the Python standard library only (no external
aligners).

## Usage

```bash
python -m affine_match --reference refs.fa --query reads.fa \
    --json out.json --csv out.csv \
    --match 2 --mismatch 1 --gap-open 3 --gap-extend 1
```

All four scoring parameters must be positive integers.

## Input

FASTA files with one or more records. Multi-line sequences are joined and
uppercased; only `A/C/G/T` bases are accepted. Empty sequences, invalid
bases, and duplicate identifiers within one file are rejected. The record
identifier is the first whitespace-delimited token of the header; the full
header text is kept as the description for reports.

## Scoring

Smith-Waterman local alignment with affine gap penalties. A contiguous gap
of length `k` costs `gap_open + (k - 1) * gap_extend`. Matches add
`match`, mismatches subtract `mismatch`, and local scores are floored at
zero.

Each query is aligned against every reference in both orientations (the
query as given, `+`, and its reverse complement, `-`). The best assignment
is the highest score; ties are broken by reference identifier
(lexicographic), then by preferring the `+` strand.

### Tie-breaking inside one alignment

Deterministic rules, applied in this order:

- cell score: diagonal (aligned pair) > gap in reference > gap in query >
  restart at zero;
- gap matrix: opening a new gap > extending an existing gap;
- best-cell scan: earliest query end, then earliest reference end.

## Coordinates and output

All intervals are 0-based, half-open `[start, end)`. Reverse-strand hits
are converted back to the original query coordinates, so `query_interval`
always refers to the query sequence as given in the input file. `coverage`
is the fraction of the query covered by the aligned part.

The JSON summary preserves query input order and includes, per query: best
reference, strand, score, reference/query intervals, coverage, and the
traceback result (gapped alignment strings plus an operations string:
`=` match, `X` mismatch, `I` gap in reference, `D` gap in query).
Re-scoring the reported alignment with the configured parameters
reproduces the reported score.

Queries whose best score is zero are reported as `matched: false` with
null intervals — no empty hits are fabricated. The CSV file is a flat
assignment table in the same query order.

## Development

```bash
python -m pytest tests
```
