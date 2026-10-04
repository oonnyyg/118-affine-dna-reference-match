# affine-dna-reference-match

Batch local alignment of short query DNA fragments against reference
sequences, implemented with the Python standard library only. For each
query it reports the best-matching reference, the alignment interval, and
a replayable gapped alignment, exported as JSON and CSV.

## Usage

```
python -m affine_align.cli \
    --reference refs.fa --query queries.fa \
    --json results.json --csv assignments.csv \
    --match 2 --mismatch 1 --gap-open 3 --gap-extend 1
```

All four scoring parameters must be positive integers; the values above
are the defaults.

## FASTA input rules

- Multi-line sequences are concatenated and converted to uppercase.
- Only `A`, `C`, `G`, `T` bases are accepted; anything else is rejected.
- Empty sequences and duplicate identifiers within one file are rejected.
- The identifier is the header text up to the first whitespace; the
  remaining header text is kept as the description for reporting.

## Scoring

Smith-Waterman local alignment with affine gap penalties:

- match: `+match`
- mismatch: `-mismatch`
- a contiguous gap of length `k`: `-(gap_open + (k - 1) * gap_extend)`

Local alignment scores are clamped at zero, so unrelated regions never
contribute negative score. A query whose best score is zero is reported
as unmatched (`matched: no`), with no fabricated interval.

Each query is aligned in both orientations (forward and reverse
complement) against every reference. The highest score wins; ties are
broken by reference identifier (lexicographic), then by strand (forward
before reverse).

### Tie-breaking inside one alignment

Deterministic rules, applied whenever candidates have equal score:

- Predecessor state preference: diagonal (match/mismatch), then gap in
  the reference, then gap in the query.
- Gap extension is preferred over opening a new gap; when opening, the
  diagonal state is preferred over the other gap state.
- The best-scoring cell is the first cell attaining the maximum score,
  scanning rows top-to-bottom and left-to-right.

## Coordinates

All intervals are 0-based half-open `[start, end)`:

- `query_start` / `query_end`: on the original query sequence, even for
  reverse-strand hits (converted back from the reverse complement).
- `ref_start` / `ref_end`: on the reference sequence.
- `coverage`: fraction of query bases covered by the aligned part
  (matched, mismatched, or inserted bases relative to the reference).

## Output

- JSON summary: scoring parameters, coordinate system note, and per-query
  entries in input order. Matched entries include `reference_id`,
  `strand`, `score`, intervals, `coverage`, and the replayable alignment
  (`aligned_query`, `aligned_ref`, `operations`). Operations are `=`
  (match), `X` (mismatch), `I` (insertion relative to the reference),
  `D` (deletion relative to the reference). Re-scoring the gapped
  alignment with the same parameters reproduces `score`.
- CSV assignment table: one row per query, in input order, with empty
  fields for unmatched queries.

## Development

```
python -m pytest tests
```
