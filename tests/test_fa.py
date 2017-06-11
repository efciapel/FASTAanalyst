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


# lenght tests
def test_lenght_is_correct():

    seq = read_fasta(FASTA_FILE_PATH)
    correct_lenght = len(seq)
    assert length(seq) == correct_lenght


# nucleotide_occurence tests
@pytest.mark.parametrize('nucleotide, occurrence, sequence', [
    ('a', 1, 'AAAAAAAAAAAA'),
    ('c', 1, 'CCCCCCCCCCCC'),
    ('g', 1, 'GGGGGGGGGGGG'),
    ('t', 1, 'TTTTTTTTTTTT'),
    ('at', 1, 'AAAAAATTTTTT'),
    ('a', 0.5, 'AAAAAATTTTTT'),
    ('t', 0.5, 'AAAAAATTTTTT'),
    ('c', 0, 'AAAAAATTTTTT'),
    ('g', 0, 'AAAAAATTTTTT'),
])
def test_nucleotide_occurrence_is_correct(nucleotide, occurrence, sequence):

    seq = read_dna(sequence)
    ret = nucleotide_occurrence(seq, nucleotide)
    assert ret == {nucleotide.upper(): occurrence}


def test_nucleotide_occurence_can_take_multiple_args():
    seq = read_dna(SEQUENCE)
    ret = nucleotide_occurrence(seq, 'a', 'g', 'ct')
    assert ret['A'] == 0.25
    assert ret['G'] == 0.25
    assert ret['CT'] == 0.5


def test_letter_case_is_not_relevant():
    seq = read_dna(SEQUENCE)
    ret = nucleotide_occurrence(seq, 'A', 'g', 'Ct')
    assert ret['A'] == 0.25
    assert ret['G'] == 0.25
    assert ret['CT'] == 0.5


# molecular_mass tests
@pytest.mark.parametrize('test_seq, right_value', [
    ('AAAA', 251.25 * 4),
    ('CCCC', 227.22 * 4),
    ('GGGG', 267.25 * 4),
    ('TTTT', 242.24 * 4),
    ('UUUU', 228.14 * 4),
    ('AATT', 251.25 * 2 + 242.24 * 2),
])
def test_molecular_mass_is_correct(test_seq, right_value):
    seq = read_dna(test_seq)
    assert molecular_mass(seq) == right_value
