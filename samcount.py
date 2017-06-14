#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Count SAM file flags
Author: James B. Pease
"""

import sys
import os
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


def explain_sam_flag(flag, code=None):
    """Decode SAM File Flags into plain text"""
    flag = int(flag)
    bin_string = str(bin(flag))[2:]
    bin_string = bin_string[::-1]
    bool_flags = []
    i = 0
    while i < len(bin_string):
        if bin_string[i] == '1':
            bool_flags.append(code[i])
        i += 1
    return bool_flags


def generate_argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("samfile",
                        help="SAM file to analyze")
    parser.add_argument("-o", "--out", type=os.path.abspath,
                        default=sys.stdout,
                        help="write to output file")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    sam_flag_code = {0:  "read_paired",
                     1:  "read_mapped_in_proper_pair",
                     2:  "read_unmapped",
                     3:  "mate_unmapped",
                     4:  "read_reverse_strand",
                     5:  "mate_reverse_strand",
                     6:  "first_in_pair",
                     7:  "second_in_pair",
                     8:  "not_primary_alignment",
                     9:  "read_fails_qc",
                     10: "read_is_PCR_or_optical_duplicate",
                     11: "supplementary_alignment"}
    flag_count = {}
    flag_bits = {}
    with open(args.samfile, 'r') as samfile:
        for line in samfile:
            if line[0] == '@':
                continue
            row = line.split()
            flag_count[row[1]] = flag_count.get(row[1], 0) + 1
    flag_codes = sorted([(v, k) for (k, v) in flag_count.items()],
                        reverse=True)
    total = sum(flag_count.values())
    if total < 1:
        raise SyntaxError("Nothing found in SAM file!")
    print("count\tflag\tpct\tinfo", file=args.outfile)
    for (count, flag) in flag_codes:
        print(("{}\t{}\t{}%\t{}"
               ).format(count, flag, round(float(count)*100/total, 1),
                        ';'.join(explain_sam_flag(flag, code=sam_flag_code))),
              file=args.outfile)
        bincode = bin(int(flag))[2:].zfill(11)[::-1]
        for j in range(len(bincode)):
            if bincode[j] == '1':
                flag_bits[j] = flag_bits.get(j, 0) + count
    print("", file=args.outfile)
    print("\t".join(["flagbit", "desc", "count", "pct"]),
          file=args.outfile)
    for samcode in sorted(flag_bits.keys()):
        print(("{}\t{}\t{}\t{}%"
               ).format(samcode, sam_flag_code[samcode], flag_bits[samcode],
                        round(float(flag_bits[samcode])*100/total, 1)),
              file=args.outfile)
    return ''


if __name__ == "__main__":
    main()
