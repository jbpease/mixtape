# -*- coding: utf-8 -*-
"""

mixTAPE: mix of Tools of Analysis of Phylogenetic and Evolution
sam2seq.py
Converts SAM to Fasta or Fastq

Version 2014-09-09: Initial Release

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

def main(arguments=sys.argv[1:]):
    """Main count_sam function"""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("samfile",
                        help="SAM file to analyze")
    parser.add_argument("--out",
                        help="write to output file (default=print to stdout)")
    parser.add_argument("--mode", choices=('fq', 'fqpair', 'fa'), default='fq',
                        help="FastQ, FastQ paired (separate files), or FastA")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
    parser.add_argument("--notflag", type=int, nargs='*',
                        help="output sequences CANNOT have this flag")
    parser.add_argument("--hasflag", type=int, nargs='*',
                        help="output sequences must have this flag")
    args = parser.parse_args(args=arguments)

    if not args.quiet:
        file_size = os.stat(args.samfile).st_size
    filtering = False
    if args.hasflag or args.notflag:
        filtering = True
    if args.mode == 'fqpair':
        outfile = open(args.out + ".p1", 'w')
        outfile2 = open(args.out + ".p2", 'w')
        outfile_un = open(args.out + ".s", 'w')
    else:
        outfile = open(args.out, 'w')
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
            output_read = True
            if filtering:
                flag = int(arr[1])
                if any(x & flag == x for x in args.notflag):
                    output_read = False
                if output_read:
                    if any(x & flag != x for x in args.hasflag):
                        output_read = False
            if output_read:
                if args.mode == 'fq':
                    outfile.write("@{!s}\n{!s}\n+\n{!s}\n".format(arr[0],
                                                                  arr[9],
                                                                  arr[10]))
                elif args.mode == 'fqpair':
                    flag = int(arr[1])
                    if 64 & flag == 64:
                        outfile.write("@{!s}/1\n{!s}\n+\n{!s}\n".format(arr[0],
                                                                        arr[9],
                                                                        arr[10]))
                    elif 128 & flag == 128:
                        outfile2.write("@{!s}/2\n{!s}\n+\n{!s}\n".format(arr[0],
                                                                        arr[9],
                                                                        arr[10]))
                    else:
                        outfile_un.write("@{!s}\n{!s}\n+\n{!s}\n".format(arr[0],
                                                                         arr[9],
                                                                         arr[10]))
                elif args.mode == 'fa':
                    outfile.write(">{!s}\n{!s}\n".format(arr[0], arr[9]))
            line = samfile.readline()
    outfile.close()
    if args.mode == 'fqpaired':
        outfile2.close()
    if not args.quiet:
        progress_meter(file_size, file_size)
    return ''


if __name__ == "__main__":
    main()
