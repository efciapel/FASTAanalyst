#coding: utf-8
import re

def read_fasta(fasta_file):
    """Odczytuje plik FASTA.
    Args:
        fasta_file (file): plik FASTA.
    """

    with open(fasta_file, 'r') as seq:
        seq = 1
        #TODO dokonczyc funkcje


def read_dna(dna):
    # Ewa Król
    """Odczytuje sekwencje DNA.
    Args:
        dna (str): sekwencja nukleotydowa.
    """


def length(dna):
    # Anna Osina
    """Zwraca długość podanej sekwencji.
    Args:
        dna (str): sekwencja nukleotydowa.
    """


def nucleotide_occurence(*args):
    # Anna Osina
    """Podaje procentowe występowanie danego nukleotydu w sekwencji lub grupy
    nukleotydów."""


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
    # Anna Osina
    """Zwraca nić masę cząsteczkową łańcucha.
    Args:
        dna (str): sekwencja nukleotydowa.
    """


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
