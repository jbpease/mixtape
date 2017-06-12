#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA Codon Aligner: Take CDS sequences, translate, protein align, untranslate

MixTAPE: Mix of Tools for Analysis in Phylogenetics and Evolution
http://www.github.org/jbpease/mixtape

@author: James B. Pease

This file is part of MixTAPE.

MixTAPE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MixTAPE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MixTAPE.  If not, see <http://www.gnu.org/licenses/>.


"""
import sys
import argparse
from itertools import groupby
from random import randint
import subprocess


def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


STANDARD_CODON_TABLE = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"
    }


def get_translation(seq):
    """Returns translated amino acids from nucleotides"""
    aa_seq = []
    for i in range(int(len(seq) / 3.0) + 1):
        codon = seq[i * 3:(i + 1)*3].upper().replace('U', 'T')
        if not codon:
            break
        if set(codon) == set(['-']):
            aa_seq.append('-')
        elif (len(set(codon) - set('ATGCU')) > 0) or (len(codon) < 3):
            aa_seq.append('X')
        else:
            aa_seq.append(STANDARD_CODON_TABLE[codon])
    return ''.join(aa_seq)


def untranslate(amino_acids, nucleotides, firststop=False):
    """Converts an aligned amino acid sequence
    back to codons while maintaining gaps"""
    nucleotides = nucleotides.strip().replace('-', '')
    codons = ''
    j = 0
    for i in range(len(amino_acids)):
        if amino_acids[i] == '*' and firststop:
            return codons
        elif amino_acids[i] == '-':
            codons += '---'
        else:
            codons += nucleotides[j:j+3]
            j += 3
    return codons


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile")
    parser.add_argument("outfile")
    args = parser.parse_args()
    tmppath = args.infile + "." + str(randint(0, 1000000)) + ".tmp"
    tmppathout = tmppath + ".out"
    nucseqs = {}
    with open(tmppath, 'w') as tempfile:
        for hdr, seq in fasta_iter(args.infile):
            nucseqs[hdr] = seq[:]
            newseq = get_translation(seq)
            tempfile.write(">{}\n{}\n".format(hdr, newseq))
    mafft_args = ['mafft --auto', tmppath, ">", tmppathout]
    proc = subprocess.Popen(' '.join(mafft_args),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)
    serr, sout = proc.communicate()
    prot_seq = {}
    maxlen = 0
    for hdr, seq in fasta_iter(tmppathout):
        new_seq = untranslate(seq, nucseqs[hdr])
        maxlen = max(maxlen, len(new_seq))
        prot_seq[hdr] = new_seq[:]
    with open(args.outfile, 'w') as outfile:
        for hdr, seq in prot_seq.items():
            outfile.write(">{}\n{}{}\n".format(
                hdr, seq, '-' * (maxlen - len(seq))))
    return ''


if __name__ == "__main__":
    main()
