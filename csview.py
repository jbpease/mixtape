#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#
# @author: James B. Pease
# @version: 1.0

import sys
import argparse


def main(arguments=None):
    arguments = arguments or sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Template script for general use')
    parser.add_argument('infile')
    parser.add_argument('-d', '--delim', default=',')
<<<<<<< HEAD
    parser.add_argument('-t', '--tsv',
                        action="store_true",
                        help="for tab-separated value")
    parser.add_argument('-n', '--nlines', type=int)
    parser.add_argument('-l', '--length', default=20,
                        type=int, help="field length")
=======
    parser.add_argument('-n', '--nlines', type=int)
>>>>>>> b0eaf4323e96abb38abce9466633faee9f515086
    parser.add_argument('--no-header', action='store_true')
    args = parser.parse_args(args=arguments)
    nlines = args.nlines + 0 if args.nlines is not None else None
    firstline = True
    if args.delim == 'TAB' or args.tsv is True:
        args.delim = '\t'
    with open(args.infile, 'rt') as infile:
        for line in infile:
            arr = line.rstrip().split(args.delim)
            print("|".join("{}{}".format(
                x, " "*(args.length - len(x)))[:args.length] for x in arr))
            if firstline is True:
                firstline = False
                print("-"*((args.length + 1)*len(arr)))
            if nlines is not None:
                nlines -= 1
                if nlines == 0:
                    break
    return ''


if __name__ == '__main__':
    main()
