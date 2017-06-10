#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
#
# @author: James B. Pease
# @version: 1.0
"""
Takes FASTA as input and outputs CSV of pairwise distances.
"""


import sys
import argparse
from itertools import combinations, groupby
from random import sample as rsample

AMBIG = {"A": "A",
         "C": "C",
         "G": "G",
         "T": "T",
         "R": "AG",
         "Y": "CT",
         "W": "AT",
         "S": "CG",
         "M": "AC",
         "K": "GT",
         "B": "CGT",
         "D": "AGT",
         "H": "ACT",
         "V": "ACG",
         "N": "ACGT"}


def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    i = -1
    for header in faiter:
        i += 1
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def calc_pwdist_nuc(seq0, seq1, mode="random2"):
    diff = 0
    total = 0
    for i, base0 in enumerate(seq0):
        base0 = base0.upper()
        base1 = seq1[i].upper()
        if (base0 not in 'RYWSMKBDHVNATGC' or 
                base1 not in 'RYWSMKBDHVNATGC') :
            continue
        if base0 in "BDHVN" or base1 in "BDHVN":
            if mode != "random4":
                continue
        if base0 in "RYWSMK" or base1 in "RYWSMK":
            if mode != "random2":
                continue
        base0 = rsample(AMBIG[base0], 1)[0]
        base1 = rsample(AMBIG[base1], 1)[0]
        total += 1
        if base0 != base1:
            diff +=1
    return float(diff / total)

def calc_pwdist_prot(seq0, seq1):
    diff = 0
    total = 0
    for i, base0 in enumerate(seq0):
        base0 = base0.upper()
        base1 = seq1[i].upper()
        if (base0 not in 'ACDEFGHIKLMNPQRSTVWY' or 
                base1 not in 'ACDEFGHIKLMNPQRSTVWY'):
            continue
        total += 1
        if base0 != base1:
            diff +=1
    return float(pwdist / total)


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile')
    parser.add_argument("-m", "--mode", choices=("strict", "random2",
                                                 "random4"),
                        default="random", 
                        help=("strict=skip ambiguous; "
                             "random2=choose random at biallelic only; "
                             "random4=choose random at any ambiguous"))
    parser.add_argument("-t", "--seqtype", choices=("nuc", "prot"),
                    default="nuc", 
                    help="nucleotide or protein alignment")
    args = parser.parse_args(args=arguments)
    headers = []
    seqs = {}
    distances = {}
    for hdr, seq in fasta_iter(args.infile):
        headers.append(hdr)
        seqs[hdr] = seq
    # Calculate pairwise distances for all pairs, store values under
    # both combinations
    for hdr0, hdr1 in combinations(headers, 2):
        if args.seqtype == "prot":
            pwdist = calc_pwdist_prot(seqs[hdr0], seqs[hdr1])
        else:
            pwdist = calc_pwdist_nuc(seqs[hdr0], seqs[hdr1], 
                                     mode=args.mode)
        distances[(hdr0, hdr1)] = pwdist + 0.0
        distances[(hdr1, hdr0)] = pwdist + 0.0
    # Add zeros on the diagonal
    for hdr in headers:
        distances[(hdr, hdr)] = 0.0
    # Print first header line
    print(",".join([""] + headers))
    # Iterate through the headers (skip the blank)
    for hdr0 in headers:
        print(",".join([hdr0] + [
            str(distances[(hdr1, hdr0)]) for hdr1 in headers]))
    return ''


if __name__ == '__main__':
    main()
