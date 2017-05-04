#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This script is meant as a general template for creating standard scripts in Python3.x
#
# @author: James B. Pease
# @version: 1.0

import sys
import argparse
from Bio import Phylo

def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Template script for general use')
    parser.add_argument('infile')
    args = parser.parse_args(args=arguments)
    tree = Phylo.read(args.infile, 'newick')
    Phylo.draw_ascii(tree)
    for clade in tree.get_nonterminals():
        subnodes = ';'.join([x.name for x in clade.get_terminals()])
        print(subnodes, "=", clade.name)

    return ''


if __name__ == '__main__':
    main()
