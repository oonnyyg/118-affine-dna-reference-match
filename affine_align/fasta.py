"""FASTA parsing with strict validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

VALID_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


def parse_fasta(path: str | Path) -> list[FastaRecord]:
    """Parse a FASTA file.

    - Multi-line sequences are concatenated and uppercased.
    - Only A/C/G/T bases are accepted.
    - Empty sequences and duplicate identifiers are rejected.
    - The identifier is the first whitespace-delimited token of the header;
      the remaining header text is kept as the description.
    """
    records: list[FastaRecord] = []
    seen: set[str] = set()
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal header, chunks
        if header is None:
            return
        sequence = "".join(chunks).upper()
        if not sequence:
            raise ValueError(f"empty sequence for record {header!r} in {path}")
        invalid = set(sequence) - VALID_BASES
        if invalid:
            raise ValueError(
                f"invalid bases {sorted(invalid)} in record {header!r} of {path}"
            )
        identifier = header.split()[0]
        description = header[len(identifier):].strip()
        if identifier in seen:
            raise ValueError(f"duplicate identifier {identifier!r} in {path}")
        seen.add(identifier)
        records.append(FastaRecord(identifier, description, sequence))
        header, chunks = None, []

    with open(path, "r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty header at line {lineno} of {path}")
            else:
                if header is None:
                    raise ValueError(
                        f"sequence data before first header at line {lineno} of {path}"
                    )
                chunks.append(line)
    flush()
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    return records
