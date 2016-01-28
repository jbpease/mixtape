#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

MixTAPE: Mix of Tools for Analysis in Phylogenetics and Evolution
http://www.github.org/jbpease/mixtape

FASTA Extract: Extract sequences from a list of headers
@author: James B. Pease

@version: 2016-01-05 - Initial release

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
    """Main method for fasta_extract"""
    parser = argparse.ArgumentParser(description="""
    Extract sequences from fasta file based on provided headers""")
    parser.add_argument("--fasta", help="input MVF file", required=True)
    parser.add_argument("--out", help="output FASTA file", required=True)
    parser.add_argument("--labels", nargs='*',
                        help="list of labels")
    parser.add_argument("--labelfiles", nargs='*',
                        help="list of files containing labels")
    parser.add_argument("--quiet", action="store_true", default=True,
                        help="suppress screen output")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2016-01-05")
        sys.exit()
    if args.labels and args.labelfiles:
        raise RuntimeError("Cannot use both --labels and --labelfiles")
    labels = args.labels or []
    for labelfile in args.labelfiles:
        with open(labelfile, 'r') as lblfile:
            for line in lblfile:
                labels.append(line.strip())
    with open(args.out, 'w') as outfile:
        for header, seq in fasta_iter(args.fasta):
            if header in labels:
                outfile.write(">{}\n{}\n".format(
                    header, seq))
    return ''

if __name__ == "__main__":
    main()
