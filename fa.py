#coding: utf-8
import re

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


def alignment(seq_1, seq_2):
    # Ewa Król
    """Zwraca nić score przyrównania dwóch sekwencji.
    Wykorzystuje algorytm Needlemana-Wunscha.
    Args:
        dna (str): sekwencja nukleotydowa.
    """
