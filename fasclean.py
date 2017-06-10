#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fasclean.py
Simple Fasta Cleanup for Standardizing Length and Labels,

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


def labelmod(label, args):
    """FASTA Label Modifications"""
    if args.firstlabel:
        label = label.split()[0]
    elif args.underscore:
        label = label.replace(" ", "_").replace("\t", "_")
    if args.remchar:
        charset = ''
        if args.remchar == 'STRICT':
            charset = "!@#$%^&*=~,.:;?'\"/|\\()[]{}<>"
        elif args.remchar == 'SAFE':
            charset = "!@#$%&*;<>\"'"
        else:
            charset = args.remchar
        label = ''.join([x for x in label if x not in charset])
    if args.labellength:
        label = label[:args.labellength]
    return label


def trimseq(seq):
    j = 0
    while j < len(seq):
        if seq[j] not in 'Nn-':
            break
        j += 1
    seq = seq[j:]
    if not seq:
        return seq
    j = len(seq)
    while j > 0:
        if seq[j - 1] not in 'Nn-':
            break
        j -= 1
    seq = seq[:j]
    return seq


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fasta", nargs=1,
                        help="input fasta file name")
    parser.add_argument("-s", "--seqlength", type=int, default=60,
                        help=("add newlines to sequences after this many"
                              " characters, set to 0 for all characters"
                              " on one line"))
    parser.add_argument("-l", "--labellength", type=int,
                        help="truncate label to this number of characters")
    parser.add_argument("-r", "--remchar",
                        help=("remove these special characters"
                              + "enter STRICT or SAFE to remove punctuation"))
    parser.add_argument("-f", "--firstlabel", action="store_true",
                        help="crop label at first whitespace (space/tab)")
    parser.add_argument("-u", "--underscore", action="store_true",
                        help="replace whitespace with underscores")
    parser.add_argument("-p", "--prefix",
                        help="replace labels with prefixes")
    parser.add_argument("-n", "--number", type=int, default=0,
                        help=("use with -p/--prefix for a fixed"
                              " number of digits with leading zeros"))
    parser.add_argument("-t", "--seqtrim", action="store_true",
                        help="trim 'N' and '-' from the ends of sequences")
    parser.add_argument("-L", "--log", default="fasclean.log",
                        help="specify path for log file")
    parser.add_argument('--version', action='version',
                        version='%(prog)s version 1')
    args = parser.parse_args(args=arguments)
    ndex = 0
    if args.prefix and (args.labellength or args.remchar):
        raise SyntaxError("Cannot use --prefix with "
                          + "-l/--labellength or -r/--remchar")
    if args.firstlabel and args.underscore:
        print("-u/--underscore is redundant with -f/--firstlabel",
              file=sys.stderr)
    logfile = open(args.log, 'w')
    with open(args.fasta[0]) as infile:
        line = infile.readline()
        seq = ''
        while line:
            if line[0] == ">":
                if seq:
                    if args.trim:
                        seq = trimseq(seq)
                    j = 0
                    while j < len(seq):
                        print(seq[j:j+args.seqlength], file=sys.stdout)
                        j += args.seqlength
                    seq = ''
                oldlabel = line.lstrip(">").strip()
                if args.prefix:
                    newlabel = ("{!s}_{!s}").format(args.prefix,
                                                    str(ndex).zfill(
                                                        args.number))
                else:
                    newlabel = labelmod(oldlabel, args)
                print(">{}".format(newlabel), file=sys.stdout)
                logfile.write(("{}\t{}\n").format(oldlabel, newlabel))
                ndex += 1
            else:
                seq += line.strip()
            line = infile.readline()
    if seq:
        if args.trim:
            seq = trimseq(seq)
        j = 0
        while j < len(seq):
            print(seq[j:j+args.seqlength] + "\n")
            j += args.seqlength
    seq = ''
    return ''


if __name__ == "__main__":
    main()
