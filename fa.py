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


def length(dna):
    """Zwraca długość podanej sekwencji.
    Args:
        dna (str): sekwencja nukleotydowa.
    """
    return len(dna)


def nucleotide_occurrence(dna, *args):
    """Podaje procentowe występowanie danego nukleotydu w sekwencji lub grupy
    nukleotydów.
    Args:
        dna (str): sekwencja nukleotydowa.
    """
    ret = dict()
    for arg in args:
        occurrence = 0
        for i in range(len(arg)):
            occurrence += dna.count(arg[i])
        occurrence /= length(dna)
        ret.update({arg.upper(): round(occurrence, 2)})
    return ret


def complement(dna):
    # Martyna Kępska
    """Zwraca nić komplementarną sekwencji.
    Args:
        dna (str): sekwencja nukleotydowa.
    """


def reverse_complement(dna):
    # Martyna Kępska
    """Zwraca odwróconą nić komplementarną sekwencji.
    Args:
        dna (str): sekwencja nukleotydowa.
    """


def molecular_mass(dna):
    """Zwraca nić masę cząsteczkową łańcucha.
    Args:
        dna (str): sekwencja nukleotydowa.
    """
    weights = {'a': 251.25,
               'c': 227.22,
               'g': 267.25,
               't': 242.24,
               'u': 228.14}

    mass = 0
    for key in weights:
        mass += dna.count(key) * weights[key]

    return mass


def melting_point(dna):
    # Martyna Kępska
    """Zwraca temperaturę topnienia sekwencji.
    Args:
        dna (str): sekwencja nukleotydowa.
    """


def alignment(dna_1, dna_2):
    # Ewa Król
    """Zwraca nić score przyrównania dwóch sekwencji.
    Wykorzystuje algorytm Needlemana-Wunscha.
    Args:
        dna (str): sekwencja nukleotydowa.
    """
