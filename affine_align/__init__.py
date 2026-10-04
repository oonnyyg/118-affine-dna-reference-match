"""Affine-gap Smith-Waterman batch DNA reference matching."""

from .align import Alignment, Scoring, reverse_complement, smith_waterman
from .fasta import FastaRecord, parse_fasta
from .matcher import MatchResult, best_match, match_all

__all__ = [
    "Alignment",
    "FastaRecord",
    "MatchResult",
    "Scoring",
    "best_match",
    "match_all",
    "parse_fasta",
    "reverse_complement",
    "smith_waterman",
]
