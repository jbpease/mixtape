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
    parser.add_argument('infile', type=os.path.abspath,
                        help="input file")
    parser.add_argument('-o', '--out', type=os.path.abspath,
                        default=sys.stdout, help="output file")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    with open(args.infile, 'r') as infile:
        for line in infile:
            print(line.rstrip(), file=args.out)
    return ''


if __name__ == '__main__':
    main()
