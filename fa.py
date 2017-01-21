#coding: utf-8
import re
import numpy as np

PATTERN = re.compile("[actgn]")


def read_fasta(fasta_file):
    """Odczytuje plik FASTA.
    Args:
        fasta_file (str): ścieżka pliku FASTA.
    """

    fasta_file = fasta_file.lower()
    try:
        with open(fasta_file, 'r') as f:
            # sekwencją jest wszystko oprócz pierwszej linii
            seq = f.read().splitlines()[1:]
    except FileNotFoundError:
        if len(PATTERN.findall(fasta_file)) == len(fasta_file):
            raise Exception("Nie znaleziono podanego pliku, a podana ścieżka"
                            "wygląda  jak sekwencja. "
                            "Spróbuj użyć funkcji 'read_dna'.")
        raise

    return ''.join(seq)


def read_dna(dna):
    """Odczytuje sekwencje DNA.
    Args:
        dna (str): sekwencja nukleotydowa.
    """
    if len(PATTERN.findall(dna)) != len(dna):
        raise Exception("Podana sekwencja posiada niepoprawne znaki. "
                        "Sekwencje mogą składać się z liter: "
                        "A, C, G, T, U i N. "
                        "Jeżeli chciałeś/aś wczytać sekwencję z pliku, "
                        "spróbuj użyć metody 'read_fasta'.")

    return dna.lower()


def sequence_type(seq):
    seq_type = 'DNA'

    if seq.find('t') != -1:
        seq_type = 'RNA'

    return seq_type


def length(seq):
    # Anna Osina
    """Zwraca długość podanej sekwencji.
    Args:
        seq (str): sekwencja nukleotydowa.
    """


def nucleotide_occurence(*args):
    # Anna Osina
    """Podaje procentowe występowanie danego nukleotydu w sekwencji lub grupy
    nukleotydów."""


def complement(seq):
    """Zwraca nić komplementarną sekwencji.
    Args:
        seq (str): sekwencja nukleotydowa.
    """
    nucleotides = {'c': 'g', 'g': 'c'}
    ret = ''
    seq_type = sequence_type(seq)

    if seq_type == 'RNA':
        nucleotides.update({'a': 'u', 'u': 'a'})
    else:
        nucleotides.update({'a': 't', 't': 'a'})

    for base in seq:
        ret += nucleotides.get(base)

    return ret


def reverse_complement(seq):
    """Zwraca odwróconą nić komplementarną sekwencji.
    Args:
        seq (str): sekwencja nukleotydowa.
    """
    return complement(seq)[::-1]


def molecular_mass(seq):
    # Anna Osina
    """Zwraca nić masę cząsteczkową łańcucha.
    Args:
        seq (str): sekwencja nukleotydowa.
    """


def melting_point(seq):
    """Zwraca temperaturę topnienia sekwencji.
    Args:
        seq (str): sekwencja nukleotydowa.
    """
    seq_type = sequence_type(seq)
    count_a = seq.count('a')
    count_c = seq.count('c')
    count_g = seq.count('g')
    count_ut = seq.count('t') if seq_type == 'DNA' else seq.count('u')

    return 4 * (count_g + count_c) + 2 * (count_a + count_ut)


def match(base_1, base_2):
    """Sprawdza czy dane nukleotydy są takie same i zwraca odpowiednią wartość
    Args:
        base_1 (str): pierwszy nukleotyd.
        base_2 (str): drugi nukleotyd.
    """

    if base_1 == base_2:
        return 1
    return -1


def alignment(dna_1, dna_2, gap):
    """Zwraca nić score przyrównania dwóch sekwencji.
    Wykorzystuje algorytm Needlemana-Wunscha.
    Args:
        dna_1 (str): pierwsza sekwencja nukleotydowa.
        dna_2 (str): druga sekwencja nukleotydowa.
        gap (int): kara za przerwę
    """

    # tworzenie macierzy F
    matrix = np.zeros((length(dna_1) + 1, length(dna_2) + 1))

    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            if row == 0 or col == 0:
                matrix[row][col] = max(row, col) * gap
            else:
                # sprawdzenie czy nukleotydy są takie same
                match_score = match(dna_1[row - 1], dna_2[col - 1])

                diagonal = match_score + matrix[row - 1][col - 1]
                top = matrix[row - 1][col] + gap
                left = matrix[row][col - 1] + gap

                matrix[row][col] = max(diagonal, top, left)
    alignment_1 = ''
    alignment_2 = ''
    i = length(dna_1)
    j = length(dna_2)

    while i > 0 and j > 0:
        match_score = match(dna_1[i - 1], dna_2[j - 1])

        diagonal = match_score + matrix[i - 1][j - 1]
        top = matrix[i - 1][j] + gap
        left = matrix[i][j - 1] + gap

        choices = {'diagonal': {'score': diagonal,
                                'i': -1,
                                'j': -1,
                                'al_1': dna_1[i - 1],
                                'al_2': dna_2[j - 1]},
                   'top': {'score': top,
                           'i': -1,
                           'j': 0,
                           'al_1': '-',
                           'al_2': dna_2[j - 1]},
                   'left': {'score': left,
                            'i': 0,
                            'j': -1,
                            'al_1': dna_1[i - 1],
                            'al_2': '-'}}

        for choice in choices:
            if matrix[i, j] == choices[choice]['score']:
                i += choices[choice]['i']
                j += choices[choice]['j']
                alignment_1 = choices[choice]['al_1'] + alignment_1
                alignment_2 = choices[choice]['al_2'] + alignment_2
