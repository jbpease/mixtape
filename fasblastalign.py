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


class BlastHit(object):
    """Blast hit record"""
    def __init__(self, queryseq=None, subjectseq=None, pident=None,
                 length=None, mismatches=None, gaps=None,
                 query_start=None, query_end=None, subject_start=None,
                 subject_end=None, evalue=None, bitscore=None):
        self.queryseq = queryseq
        self.subjectseq = subjectseq
        self.pident = float(pident) / 100
        self.length = int(length)
        self.mismatches = int(mismatches)
        self.gaps = int(gaps)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.subject_start = int(subject_start)
        self.subject_end = int(subject_end)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)

    def orientation(self):
        if ((self.query_start > self.query_end
             and self.subject_end > self.subject_start) or
            (self.query_end > self.query_start and
             self.subject_start > self.subject_end)):
            return -1
        return 1


class BlastSet(object):
    """Set of blast hits"""
    def __init__(self, hits=None):
        self.hits = [] if hits is None else hits

    def __len__(self):
        return len(self.hits)

    def __getitem__(self, k):
        return self.hits[k]

    def __iter__(self):
        return iter(self.hits)

    def add_tabular(self, arr, ncol=12):
        """Reads a standard 12-column blast input as list"""
        xstart = 0 if ncol == 12 else 2
        self.hits.append(BlastHit(
            queryseq=arr[0],
            subjectseq=arr[1],
            pident=arr[2 + xstart],
            length=arr[3 + xstart],
            mismatches=arr[4 + xstart],
            gaps=arr[5 + xstart],
            query_start=arr[6 + xstart],
            query_end=arr[7 + xstart],
            subject_start=arr[8 + xstart],
            subject_end=arr[9 + xstart],
            evalue=arr[10 + xstart],
            bitscore=arr[11 + xstart]))
        return ''

    def parse_blast(self, blastfile, ncol=12):
        """Parse 12-column tabular BLAST file"""
        with open(blastfile, 'r') as blastfile:
            for line in blastfile:
                if line[0] == '#':
                    continue
                if len(line) < 10:
                    continue
                self.add_tabular(line.rstrip().split(), ncol=ncol)
        return ''

    def remove_low_pident(self, cutoff, verbose=True):
        """Remove hits with low_pident"""
        if cutoff >= 0 and cutoff <= 1:
            remlist = []
            for i in range(len(self.hits)):
                if self.hits[i].pident < cutoff:
                    remlist.append(i)
        remlist.sort(reverse=True)
        for i in remlist:
            if verbose:
                print("{} removed due to low pid").format(
                    self.hits[i].subjectseq)
            self.hits.pop(i)
        return ''

    def remove_superceded(self, verbose=False):
        """Remove hits where the another
           sequence completely contains the region"""
        remlist = []
        for i in range(len(self.hits) - 1):
            if i in remlist:
                continue
            for j in range(i + 1, len(self.hits)):
                if overlap_range([self.hits[i].query_start,
                                  self.hits[i].query_end],
                                 [self.hits[j].query_start,
                                  self.hits[j].query_end]) == 2:
                    remlist.append(i)
                    break
                if overlap_range([self.hits[j].query_start,
                                  self.hits[j].query_end],
                                 [self.hits[i].query_start,
                                  self.hits[i].query_end]) == 2:
                    remlist.append(j)
        remlist.sort(reverse=True)
        for i in remlist:
            if verbose:
                print("{} superceded".format(self.hits[i].subjectseq))
            self.hits.pop(i)
        return ''

    def n_subjectseq(self):
        """Return number of unique database sequences"""
        return len(set([hit.subjectseq for hit in self.hits]))

    def get_subjectseqs(self):
        """Returns a list of database sequence names"""
        retlist = list(set([hit.subjectseq for hit in self.hits]))
        retlist.sort()
        return retlist

    def get_hits(self, subject='', query=''):
        """Search for specific hits by database name or query name"""
        retlist = []
        for hit in self.hits:
            if ((not subject or subject in hit.subjectseq)
                    and (not query or query in hit.queryseq)):
                retlist.append(hit)
        return retlist


def overlap_range(range0, range1, inclusive=True):
    """
    Tests whether two ranges of integers overlap,
    -2 = range0 inside range1
    -1 = the end of range1 overlaps the start of range0
    0 = no overlap
    1 = the end of range0 overlaps the start of range1
    2 = range1 inside range0
    """
    range0 = [min(range0), max(range0)]
    range1 = [min(range1), max(range1)]
    if inclusive:
        range0 = [range0[0] - 1, range0[1] + 1]
        range1 = [range1[0] - 1, range1[1] + 1]
    if min(range0) > min(range1) and min(range0) < max(range1):
        return 2 if max(range0) < max(range1) else -1
    elif max(range0) > min(range1) and max(range0) < max(range1):
        return -2 if min(range0) > min(range1) else 1
    return 0


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


def complement(seq):
    return seq[::-1].lower().replace(
        'a', 'T').replace('t', 'A').replace(
            'c', 'G').replace('g', 'C')


def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument('-f', "--fasta", type=os.path.abspath,
                        help="input FASTA file")
    parser.add_argument('-b', "--blast", type=os.path.abspath,
                        help="input BLAST tabular tab-delimted file")
    parser.add_argument('-o', '--out', type=os.path.abspath,
                        default=sys.stdout, help="output FASTA file")
    parser.add_argument('-n', '--ncol', type=int,
                        default=12, help="number of BLAST columns")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    seqs = {}
    for hdr, seq in fasta_iter(args.fasta):
        if ":" in hdr:
            hdr = hdr[:hdr.find(":")]
        seqs[hdr] = seq
    newseqs = []
    blastdata = BlastSet()
    blastdata.parse_blast(args.blast, ncol=args.ncol)
    seqnames = set([])
    for blasthit in blastdata:
        if blasthit.queryseq not in seqnames:
            seqnames.update([blasthit.queryseq])
            newseqs.append((-20000, blasthit.queryseq,
                            seqs[blasthit.queryseq]))
        if blasthit.subjectseq not in seqnames:
            xseq = seqs[blasthit.subjectseq][:]
            orientation = blasthit.orientation()
            if orientation == 1:
                xseq = ("-" * (blasthit.query_start - 1)) + xseq
                newseqs.append((blasthit.query_start, blasthit.subjectseq,
                                xseq))
                seqnames.update([blasthit.subjectseq])
            elif orientation == -1:
                xseq = ("-" * (blasthit.query_start - 1)) + complement(xseq)
                newseqs.append((blasthit.query_start, blasthit.subjectseq,
                               xseq))
                seqnames.update([blasthit.subjectseq])
    maxlen = max(len(x[2]) for x in newseqs)
    with open(args.out, 'w') as outfile:
        for _, hdr, seq in sorted(newseqs):
            outfile.write(">{}\n{}{}\n".format(hdr, seq, '-'*(maxlen -
                                                              len(seq))))
    return ''


if __name__ == '__main__':
    main()
