"""Minimal strict FASTA reader (standard library only)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

VALID_BASES = frozenset("ACGT")


class FastaError(ValueError):
    """Raised when a FASTA file is malformed."""


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA entry.

    ``ident`` is the first whitespace-delimited token of the header;
    ``description`` keeps the full header text for reporting.
    """

    ident: str
    description: str
    sequence: str


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Parse a FASTA file, enforcing ACGT-only, non-empty, unique ids."""
    records: list[FastaRecord] = []
    seen: set[str] = set()
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal header, chunks
        if header is None:
            return
        ident = header.split()[0]
        sequence = "".join(chunks).upper()
        if not sequence:
            raise FastaError(f"empty sequence for record {ident!r} in {path}")
        bad = set(sequence) - VALID_BASES
        if bad:
            raise FastaError(
                f"invalid bases {sorted(bad)} in record {ident!r} in {path}"
            )
        if ident in seen:
            raise FastaError(f"duplicate identifier {ident!r} in {path}")
        seen.add(ident)
        records.append(FastaRecord(ident, header, sequence))
        header, chunks = None, []

    with open(path, "r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                if not header:
                    raise FastaError(f"empty header at {path}:{lineno}")
            else:
                if header is None:
                    raise FastaError(
                        f"sequence data before first header at {path}:{lineno}"
                    )
                chunks.append(line)
    flush()
    if not records:
        raise FastaError(f"no FASTA records found in {path}")
    return records
