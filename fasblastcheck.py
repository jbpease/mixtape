#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""
Template Script
James B. Pease
"""

import os
import sys
import argparse
from mixcore import fasta_iter, BlastSet

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


def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("--query", type=os.path.abspath,
                        help="input of query sequence FASTA file")
    parser.add_argument("--db", type=os.path.abspath,
                        help="input of db sequences file")
    parser.add_argument("--blast", type=os.path.abspath,
                        help="input BLAST tabular tab-delimted file")
    parser.add_argument('--out', type=os.path.abspath,
                        default=sys.stdout, help="output FASTA file")
    parser.add_argument('--ncol', type=int,
                        default=14, help="number of BLAST columns")
    parser.add_argument('--evalue', type=float,
                        default=1e-10, help="number of BLAST columns")
    parser.add_argument('--pid', type=float,
                        default=0.90, help="number of BLAST columns")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    db_seqs = set([])
    query_seqs = set([])
    for hdr, seq in fasta_iter(args.db):
        hdr = hdr.split()[0]
        if ":" in hdr:
            hdr = hdr[:hdr.find(":")]
        db_seqs.update([hdr])
    for hdr, seq in fasta_iter(args.query):
        hdr = hdr.split()[0]
        if ":" in hdr:
            hdr = hdr[:hdr.find(":")]
        query_seqs.update([hdr])
    blastdata = BlastSet()
    blastdata.parse_blast(args.blast, ncol=args.ncol)
    matches = set([])
    for blasthit in blastdata:
        if blasthit.evalue > args.evalue:
            continue
        matches.update([blasthit.queryseq])
        matches.update([blasthit.subjectseq])
    with open(args.out + ".not_in_db", 'w') as outfile:
        for seqname in db_seqs - matches:
            outfile.write(seqname + "\n")
    with open(args.out + ".not_in_query", 'w') as outfile:
        for seqname in query_seqs - matches:
            outfile.write(seqname + "\n")

    return ''


if __name__ == '__main__':
    main()
