import pytest

from source.fa import *


FASTA_FILE_PATH = 'tests/sequence'
SEQUENCE = 'AAAAACCCCCTTTTTGGGGG'  # 5 x A, 5 x C, 5 x T, 5 x G
WRONG_SEQUENCE = 'VERY_WRONG_SEQUENCE_WITH_INVALID_CHARACTERS'


# read_fasta tests
def test_function_does_not_read_first_line():

    with open(FASTA_FILE_PATH, 'r') as seq:
        read_file = seq.read().splitlines()

    read_seq = read_fasta(FASTA_FILE_PATH).splitlines()

    assert read_file[0].startswith('>')
    assert not read_seq[0].startswith('>')
    assert read_seq[0].startswith(read_file[1].lower())


def test_sequence_raises_exception():

    with pytest.raises(Exception):
        read_fasta(SEQUENCE)


def test_wrong_path_raises_file_not_found():

    wrong_path = 'wrong/path'
    with pytest.raises(FileNotFoundError):
        read_fasta(wrong_path)


# read_fasta and read_dna tests
@pytest.mark.parametrize('function, correct_input', [
    (read_fasta, FASTA_FILE_PATH),
    (read_dna, SEQUENCE)
])
def test_sequence_is_string(function, correct_input):
    seq = function(correct_input)

    assert isinstance(seq, str)


@pytest.mark.parametrize('function, correct_input', [
    (read_fasta, FASTA_FILE_PATH),
    (read_dna, SEQUENCE)
])
def test_sequence_is_lowercase(function, correct_input):
    seq = function(correct_input)

    assert seq.islower()


# read_dna tests
def test_wrong_sequence_raises_exception():

    with pytest.raises(Exception):
        read_dna(WRONG_SEQUENCE)
