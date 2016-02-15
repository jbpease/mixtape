#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
#
# This script is meant as a general template for creating standard scripts
# It also includes features that make it Python 2.6+/3.x compatible.
#
# @author: James B. Pease
# @version: 1.0

from __future__ import print_function, unicode_literals
import sys
import argparse


def main(arguments=None):
    arguments = arguments or sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Template script for general use')
    parser.add_argument('infile')
    parser.add_argument('outfile')
    args = parser.parse_args(args=arguments)
    with open(args.infile, 'rt') as infile, open(
            args.outfile, 'wt') as outfile:
        for line in infile:
            outfile.write(line)

    return ''

if __name__ == '__main__':
    main()
