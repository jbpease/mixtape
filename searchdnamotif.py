#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
#
#
# @author: James B. Pease
# @version: 1.0
"""
Search for DNA motif using a simple string with ambiguous nucleotide characters 
"""


import sys
import argparse
import re
from itertools import groupby

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


def process_code(pattern):
    ambig = {"A": "A",
             "C": "C",
             "G": "G",
             "T": "T",
             "R": "[AG]",
             "Y": "[CT]",
             "W": "[AT]",
             "S": "[CG]",
             "M": "[AC]",
             "K": "[GT]",
             "B": "[CGT]",
             "D": "[AGT]",
             "H": "[ACT]",
             "V": "[ACG]",
             "N": "[ACGT]"}
    new_pattern = ['.*']
    pattern = pattern.upper()
    for i in range(len(pattern)):
        new_pattern.append(ambig[pattern[i]])
    new_pattern.append('.*')
    return ''.join(new_pattern)



def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', help="input fasta file", required=True)
    parser.add_argument('-p', '--pattern', required=True,
                        help="pattern with ambiguous nucleotides")
    args = parser.parse_args(args=arguments)
    pattern = re.compile(process_code(args.pattern.upper()))
    for hdr, seq in fasta_iter(args.fasta):
        if re.search(pattern, seq):
            print(">{}\n{}\n".format(hdr,seq))
    return ''


if __name__ == '__main__':
    main()
