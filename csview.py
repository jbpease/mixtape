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
    parser.add_argument('-n', '--nlines', type=int)
    parser.add_argument('--no-header', action='store_true')
    args = parser.parse_args(args=arguments)
    nlines = args.nlines + 0 if args.nlines is not None else None
    firstline = True
    with open(args.infile, 'rt') as infile:
        for line in infile:
            arr = line.rstrip().split(args.delim)
            print("|".join("{}{}".format(x, " "*(20-len(x)))[:20] for x in arr))
            if firstline is True:
                firstline = False
                print("-"*(21*len(arr)))
            if nlines is not None:
                nlines -= 1
                if nlines == 0:
                    break
    return ''


if __name__ == '__main__':
    main()
