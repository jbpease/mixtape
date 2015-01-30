# -*- coding: utf-8 -*-
"""

mixTAPE: mix of Tools of Analysis of Phylogenetic and Evolution
count_sam.py
Basic analysis of SAM file to determine mapping efficiency.

Version 2014-07-25: Initial Release

"""

from __future__ import print_function
import sys, os, argparse

def progress_meter(i, target):
    """Tracks Progress on a Numeric Scale"""
    tlen = len(str(target))
    pct = str(((i) * 100)/target)
    msg = (str(i).zfill(tlen) + "/" + str(target) + "=" + pct.zfill(3)  + '%')
    if i > 0 and i < target:
        sys.stderr.write('\b' * len(msg) + msg)
    elif i == 0:
        sys.stderr.write(msg)
    elif i == target - 1 or i == target:
        sys.stderr.write("\b" * len(msg) + msg + " complete.\n")
    else:
        sys.stderr.write("ERROR")
    return ''

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

def main(arguments=sys.argv[1:]):
    """Main count_sam function"""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("samfile",
                        help="SAM file to analyze")
    parser.add_argument("--out",
                        help="write to output file (default=print to stdout)")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
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
    if not args.quiet:
        file_size = os.stat(args.samfile).st_size
    with open(args.samfile, 'r') as samfile:
        if not args.quiet:
            j = 1000000
            progress_meter(0, file_size)
        line = samfile.readline()
        while line:
            if not args.quiet:
                j -= 1
                if not j:
                    j = 1000000
                    progress_meter(samfile.tell(), file_size)
            if line[0] == '@':
                line = samfile.readline()
                continue
            arr = line.split()
            flag_count[arr[1]] = flag_count.get(arr[1], 0) + 1
            line = samfile.readline()
    if not args.quiet:
        progress_meter(file_size, file_size)
    flag_codes = sorted([(v, k) for (k, v) in flag_count.items()],
                        reverse=True)
    total = sum(flag_count.values())
    if not total:
        raise SyntaxError("Nothing found in SAM file!")
    if not args.out:
        outfile = sys.stdout
    else:
        outfile = open(args.out, 'w')
    print("Count\tFlag\tPct\tInfo", file=outfile)
    for (count, flag) in flag_codes:
        print(("{!s}\t{!s}\t{!s}%\t{!s}"
              ).format(count, flag, round(float(count)*100/total, 1),
                       ';'.join(explain_sam_flag(flag, code=sam_flag_code))),
              file=outfile)
        bincode = bin(int(flag))[2:].zfill(11)[::-1]
        for j in xrange(len(bincode)):
            if bincode[j] == '1':
                flag_bits[j] = flag_bits.get(j, 0) + count
    print("", file=outfile)
    print("\t".join(["FlagBit", "Explanation", "Count", "Pct"]),
          file=outfile)
    for samcode in sorted(flag_bits.keys()):
        print(("{!s}\t{!s}\t{!s}\t{!s}%"
              ).format(samcode, sam_flag_code[samcode], flag_bits[samcode],
                       round(float(flag_bits[samcode])*100/total, 1)),
              file=outfile)
    outfile.close()
    return ''


if __name__ == "__main__":
    main()
