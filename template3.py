#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
#
# This script is meant as a general template for Py3
#
# @author: James B. Pease
# @version: 1.0

import sys
import argparse


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile')
    parser.add_argument('outfile')
    args = parser.parse_args(args=arguments)
    with open(args.infile, 'r') as infile, open(
            args.outfile, 'w') as outfile:
        for line in infile:
            outfile.write(line)
    return ''


if __name__ == '__main__':
    main()
