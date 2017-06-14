#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8

"""
Count kmer strings from FASTA file
James B. Pease
"""

import os
import sys
import argparse
from itertools import groupby

_LICENSE = """
http://www.github.org/jbpease/mixtape
MixTAPE: Mix of Tools for Analysis in Phylogenetics and Evolution

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


def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("fasta", help="input FASTA file",
                        type=os.path.abspath)
    parser.add_argument("-o", "--out", help="output kmer output",
                        default=sys.stdout, type=os.path.abspath)
    parser.add_argument("-k", "--klen", type=int, default=8,
                        help="length of kmer [8]")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    kmers = {}
    for _, seq in fasta_iter(args.fasta):
        for i in range(len(seq) - args.klen):
            kmers[seq[i:i + args.klen]] = kmers.get(
                seq[i:i + args.klen], 0) + i
    for count, kmer in sorted(
            [(v, k) for (k, v) in kmers.items()], reverse=True):
        print("{}\t{}\n".format(count, kmer), file=args.out)
    return ''


if __name__ == '__main__':
    main()
