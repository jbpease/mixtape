#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FASTA Extract: Extract sequences from a list of headers

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
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", help="input fasta file", nargs=1)
    parser.add_argument("-l", "--labels", nargs=1,
                        help="comma-separated list of labels")
    parser.add_argument("-L", "--labelfile", nargs=1,
                        help="file containing labels (one-per-line)")
    parser.add_argument('--version', action='version',
                        version='%(prog)s version 1')
    args = parser.parse_args(args=arguments)
    if args.labels and args.labelfiles:
        raise RuntimeError("Cannot use both -l/--label and -L/--labelfile")
    labels = args.labels[:] if args.labels is not None else []
    for labelfile in args.labelfiles:
        with open(labelfile, 'r') as lblfile:
            for line in lblfile:
                labels.append(line.strip())
        for header, seq in fasta_iter(args.fasta[0]):
            if header in labels:
                print(">{}\n{}".format(header, seq), file=sys.stdout)
    return ''


if __name__ == "__main__":
    main()
