"""
fasta_iter.py

This script written by Brent Pedersen (UUtah), and taken from Biostars
https://www.biostars.org/p/710/#1412
"""

from itertools import groupby

def fasta_iter(fasta_name, noiter=False):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

