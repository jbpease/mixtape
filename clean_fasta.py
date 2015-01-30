# -*- coding: utf-8 -*-
"""

mixTAPE: mix of Tools of Analysis of Phylogenetic and Evolution
James B. Pease
http://www.bitbucket.org/jbpease/mixtape

clean_fasta.py: Simple tool to cleanup a FASTA file:
- Remove problematic puctuation symbols from the headers
- Remove additional words from the header lines
- Standardize length of sequence lines

Version 2014-07-25: Initial Release

"""

import sys, argparse

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

def main(arguments=sys.argv[1:]):
    """Main clean_fasta process"""
    parser = argparse.ArgumentParser(description="""
    Simple Fasta Cleanup for Standardizing Length and Labels,
    use --remchar to remove symbols, --labellength to crop the label to
    a maximum lengthg, --firstlabe""")
    parser.add_argument("infasta",
                        help="input fasta file name")
    parser.add_argument("outfasta",
                        help="output fasta name")
    parser.add_argument("-s", "--seqlength", type=int, default=60)
    parser.add_argument("-l", "--labellength", type=int)
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
                        help="use with prefix for a fixed number of digits")
    parser.add_argument("-t", "--trim", action="store_true",
                        help="trim 'N' and '-' from the ends of sequences")
    parser.add_argument("--log",
                        help="specify path for log file (defaut=outfasta.log)")

    args = parser.parse_args(args=arguments)
    ndex = 0
    outfile = open(args.outfasta, 'w')
    if args.prefix and (args.labellength or args.remchar):
        raise SyntaxError("Cannot use --prefix with "
                          + "--labellength or --remchar")
    if args.firstlabel and args.underscore:
        raise SyntaxError("--underscore is redundant with --firstlabel")
    if not args.log:
        args.log = args.outfasta + ".log"
    logfile = open(args.log, 'w')
    with open(args.infasta) as infile:
        line = infile.readline()
        seq = ''
        while line:
            if line[0] == ">":
                if seq:
                    if args.trim:
                        seq = trimseq(seq)
                    j = 0
                    while j < len(seq):
                        outfile.write(seq[j:j+args.seqlength] + "\n")
                        j += args.seqlength
                    seq = ''
                oldlabel = line.lstrip(">").strip()
                if args.prefix:
                    newlabel = ("{!s}_{!s}").format(args.prefix,
                                                    str(ndex).zfill(
                                                        args.number))
                else:
                    newlabel = labelmod(oldlabel, args)
                outfile.write(">{!s}\n".format(newlabel))
                logfile.write(("{!s}\t{!s}\n").format(oldlabel, newlabel))
                ndex += 1
            else:
                seq += line.strip()
            line = infile.readline()
    if seq:
        if args.trim:
            seq = trimseq(seq)
        j = 0
        while j < len(seq):
            outfile.write(seq[j:j+args.seqlength] + "\n")
            j += args.seqlength
    seq = ''
    outfile.close()
    return ''

if __name__ == "__main__":
    main()
