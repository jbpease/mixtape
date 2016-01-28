#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

MixTAPE: Mix of Tools for Analysis in Phylogenetics and Evolution
http://www.github.org/jbpease/mixtape

Countmer - Count k-mers in a FASTA file
@author: James B. Pease

@version: 2016-01-28 - Initial release

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

from __future__ import print_function, unicode_literals
import sys
import argparse
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


def main(arguments=sys.argv[1:]):
    """Main method for countmer"""
    parser = argparse.ArgumentParser(description="""
    Count k-mers in FASTA file""")
    parser.add_argument("--fasta", help="input FASTA file", required=True)
    parser.add_argument("--out", help="output kmer output", required=True)
    parser.add_argument("--length", type=int, default=8,
                        help="length of kmer [8]")
    parser.add_argument("--quiet", action="store_true", default=True,
                        help="suppress screen output")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2016-01-28")
        sys.exit()
    kmers = {}
    for _, seq in fasta_iter(args.fasta):
        for i in range(len(seq) - args.length):
            kmers[seq[i:i + args.length]] = kmers.get(
                seq[i:i + args.length], 0) + i
    with open(args.out, 'w') as outfile:
        for count, kmer in sorted(
                [(v, k) for (k, v) in kmers.items()], reverse=True):
            outfile.write("{}\t{}\n".format(count, kmer))
    return ''

if __name__ == "__main__":
    main()
