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
from mixcore import fasta_iter, BlastSet, complement

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
    parser.add_argument('--pid', type=float,
                        default=0.90, help="number of BLAST columns")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    seqs = {}
    for hdr, seq in fasta_iter(args.db):
        hdr = hdr.split()[0]
        if ":" in hdr:
            hdr = hdr[:hdr.find(":")]
        seqs[hdr] = seq
    for hdr, seq in fasta_iter(args.query):
        hdr = hdr.split()[0]
        if ":" in hdr:
            hdr = hdr[:hdr.find(":")]
        seqs[hdr] = seq
        query_hdr = hdr[:]
    newseqs = []
    blastdata = BlastSet()
    blastdata.parse_blast(args.blast, ncol=args.ncol)
    seqnames = set([])
    min_coord = 0
    for blasthit in blastdata:
        if blasthit.pident < args.pid:
            continue
        if blasthit.queryseq not in seqnames:
            seqnames.update([blasthit.queryseq])
            newseqs.append((0, blasthit.queryseq,
                            seqs[blasthit.queryseq]))
        if blasthit.subjectseq not in seqnames:
            xseq = seqs[blasthit.subjectseq][:]
            orientation = blasthit.orientation()
            if orientation == -1:
                xseq = complement(xseq)
                start_coord = (blasthit.query_start - len(xseq) +
                               blasthit.subject_start - 1)
                prefix = "-"
            else:
                prefix = "+"
                start_coord = blasthit.query_start - (
                    blasthit.subject_start)
            newseqs.append((start_coord, prefix + blasthit.subjectseq,
                           xseq))
            seqnames.update([blasthit.subjectseq])
            if start_coord < min_coord:
                min_coord = 0 + start_coord
    max_coord = max(x[0] + len(x[2]) for x in newseqs)
    consensus = {}
    with open(args.out, 'w') as outfile, open(
            args.out + ".con.fa", "w") as confile:
        for coord, hdr, seq in sorted(newseqs):
            outfile.write(">{}\n{}{}{}\n".format(
                hdr,
                '-' * (coord - min_coord),
                seq,
                '-' * (max_coord - coord - len(seq))))

            if hdr == query_hdr:
                confile.write(">{}\n{}{}{}\n".format(
                    hdr,
                    '-' * (coord - min_coord),
                    seq,
                    '-' * (max_coord - coord - len(seq))))
                continue
            for i in range(len(seq)):
                xcoord = i + coord - min_coord
                if xcoord not in consensus:
                    consensus[xcoord] = {}
                consensus[xcoord][seq[i]] = (
                    consensus[xcoord].get(seq[i], 0) + 1)
        conseq = []
        print(consensus)
        for i in range(max(consensus) + 1):
            if i not in consensus:
                conseq.append("-")
            else:
                print(list(sorted((v, k) for k, v in
                                  consensus[i].items() if k != '-')))
                conseq.append(
                    list(sorted((v, k) for k, v in
                                consensus[i].items() if k != '-'))[-1][1])
        outfile.write(">CONSENSUS\n{}{}".format(''.join(conseq),
                                                "-"*(max_coord - len(conseq))))
        confile.write(">CONSENSUS\n{}{}".format(''.join(conseq),
                                                "-"*(max_coord - len(conseq))))
    return ''


if __name__ == '__main__':
    main()
