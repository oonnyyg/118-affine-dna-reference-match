import pytest

from affine_align.fasta import parse_fasta


def write(tmp_path, text):
    path = tmp_path / "in.fa"
    path.write_text(text, encoding="utf-8")
    return path


def test_multiline_and_uppercase(tmp_path):
    path = write(tmp_path, ">seq1 some description\nac\ngt\n")
    (record,) = parse_fasta(path)
    assert record.identifier == "seq1"
    assert record.description == "some description"
    assert record.sequence == "ACGT"


def test_rejects_invalid_base(tmp_path):
    path = write(tmp_path, ">a\nACGN\n")
    with pytest.raises(ValueError, match="invalid bases"):
        parse_fasta(path)


def test_rejects_empty_sequence(tmp_path):
    path = write(tmp_path, ">a\n>b\nACGT\n")
    with pytest.raises(ValueError, match="empty sequence"):
        parse_fasta(path)


def test_rejects_duplicate_identifier(tmp_path):
    path = write(tmp_path, ">a first\nACGT\n>a second\nACGT\n")
    with pytest.raises(ValueError, match="duplicate identifier"):
        parse_fasta(path)


def test_rejects_missing_records(tmp_path):
    path = write(tmp_path, "ACGT\n")
    with pytest.raises(ValueError):
        parse_fasta(path)
