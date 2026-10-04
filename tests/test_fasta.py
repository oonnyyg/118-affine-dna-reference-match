import pytest

from affine_match.fasta import FastaError, read_fasta


def write(tmp_path, text):
    path = tmp_path / "in.fa"
    path.write_text(text, encoding="utf-8")
    return path


def test_multiline_uppercase_and_description(tmp_path):
    path = write(tmp_path, ">seq1 some description\nac\ngt\n")
    (rec,) = read_fasta(path)
    assert rec.ident == "seq1"
    assert rec.description == "seq1 some description"
    assert rec.sequence == "ACGT"


def test_rejects_invalid_base(tmp_path):
    with pytest.raises(FastaError):
        read_fasta(write(tmp_path, ">a\nACGN\n"))


def test_rejects_empty_sequence(tmp_path):
    with pytest.raises(FastaError):
        read_fasta(write(tmp_path, ">a\n>b\nACGT\n"))


def test_rejects_duplicate_ids(tmp_path):
    with pytest.raises(FastaError):
        read_fasta(write(tmp_path, ">a one\nACGT\n>a two\nACGT\n"))


def test_rejects_data_before_header(tmp_path):
    with pytest.raises(FastaError):
        read_fasta(write(tmp_path, "ACGT\n>a\nACGT\n"))
