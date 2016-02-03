#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

MixTAPE: Mix of Tools for Analysis in Phylogenetics and Evolution
http://www.github.org/jbpease/mixtape

FASTA Stats: Quick FASTA Stats
@author: James B. Pease

@version: 2016-02-03 - Initial release

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
from math import log


def mean(x):
    if len(x) == 0:
        return 0
    else:
        return float(sum(x))/len(x)


def median(x, presorted=False):
    if len(x) == 0:
        return 0
    if not presorted:
        x = sorted(x)
    if len(x) % 2 == 1:
        return x[int((len(x) + 1) / 2)]
    else:
        return float(x[int(len(x) / 2)] + x[int(len(x) / 2) + 1]) / 2.0


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


def mutual_info(elems, subelems):
    """Given an element and subelement of half the size, calculate MI"""
    for x in subelems:
        k = len(x)
        break
    total_elem = float(sum(elems.values()))
    total_subelem = float(sum(subelems.values()))
    mutual_information = 0
    if not total_elem:
        return 0.0
    for elem in elems:
        p_elem = float(elems[elem]) / total_elem
        mutual_information += (p_elem * log(p_elem / (
            (float(subelems[elem[:k]]) / total_subelem) *
            (float(subelems[elem[k:2 * k]]) / total_subelem))))
    return mutual_information


def main(arguments=sys.argv[1:]):
    """Main method for fasta_stats"""
    parser = argparse.ArgumentParser(description="""
    Calculate basic statistics for one or more fasta files""")
    parser.add_argument("fasta", nargs='*', help="input MVF file")
    parser.add_argument("--out", help="output stats file")
    parser.add_argument("--version", help="version information")
    args = parser.parse_args(args=arguments)
    characters = {}
    lengths = []
    dimers = {}
    tetramers = {}
    if args.version:
        print("Version 2016-02-03")
        sys.exit()
    entries = []
    for fastapath in args.fasta:
        entry = {'filename': fastapath}
        for header, seq in fasta_iter(fastapath):
            lengths.append(len(seq))
            for i in range(len(seq)):
                characters[seq[i]] = characters.get(seq[i], 0) + 1
                if len(seq[i:i + 2]) == 2:
                    dimers[seq[i:i + 2]] = dimers.get(seq[i:i + 2], 0) + 1
                if len(seq[i:i + 4]) == 4:
                    tetramers[seq[i:i + 4]] = tetramers.get(
                       seq[i:i + 4], 0) + 1
        entry['mi2'] = mutual_info(tetramers, dimers)
        entry['len.min'] = min(lengths)
        entry['len.max'] = max(lengths)
        entry['len.mean'] = mean(lengths)
        lengths.sort()
        entry['len.median'] = median(lengths, presorted=True)
        lensum = sum(lengths)
        runtotal = 0
        for i in range(len(lengths)):
            runtotal += lengths[i]
            if runtotal >= lensum:
                entry['len.N50'] = lengths[i]
                break
        for char in characters:
                entry['chr{}'.format(char)] = characters[char]
        entries.append(entry)
        print("finished processing {}".format(fastapath))
    headers = ["filename"]
    headers.extend(["chr{}".format(x) for x in sorted(characters)])
    headers.extend(["len.min", "len.mean", "len.max", "len.median", "len.N50"])
    headers.append("mi2")
    print("{}".format("\t".join(headers)))
    for entry in entries:
        print("\t".join([str(entry.get(x, "0")) for x in headers]))

    return ''

if __name__ == "__main__":
    main()
