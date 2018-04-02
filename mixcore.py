# -*- coding: utf-8 -*-

import os
from itertools import groupby


class CSVfile(object):
    "field-delimited data file"

    def __init__(self, fpath, delimiter=",",
                 colheaders=True):
        self.path = os.path.abspath(fpath)
        self.delimiter = delimiter
        self.colheaders = colheaders
        if self.colheaders is True:
            with open(self.path) as csvfile:
                self.colheaders = csvfile.readline().split(
                    self.delimiter)

    def __iter__(self):
        firstline = True
        with open(self.path) as csvfile:
            for line in csvfile:
                if self.colheaders is not False and firstline is True:
                    firstline = False
                    continue
                yield line.rstrip().split(self.delimiter)


STANDARD_CODON_TABLE = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"
    }


def complement(seq):
    return seq[::-1].lower().replace(
        'a', 'T').replace('t', 'A').replace(
            'c', 'G').replace('g', 'C')


def translate(seq):
    """Returns translated amino acids from nucleotides"""
    aa_seq = []
    for i in range(int(len(seq) / 3.0) + 1):
        codon = seq[i * 3:(i + 1)*3].upper().replace('U', 'T')
        if not codon:
            break
        if all(x == "-" for x in codon):
            aa_seq.append('-')
        elif (len(set(codon) - set('ATGCU')) > 0) or (len(codon) < 3):
            aa_seq.append('X')
        else:
            aa_seq.append(STANDARD_CODON_TABLE[codon])
    return ''.join(aa_seq)


def untranslate(amino_acids, nucleotides, firststop=False):
    """Converts an aligned amino acid sequence
    back to codons while maintaining gaps"""
    nucleotides = nucleotides.strip().replace('-', '')
    codons = ''
    j = 0
    for i in range(len(amino_acids)):
        if amino_acids[i] == '*' and firststop:
            return codons
        elif amino_acids[i] == '-':
            codons += '---'
        else:
            codons += nucleotides[j:j+3]
            j += 3
    return codons


def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

def fasta_text_iter(filehandler):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

def overlap_range(range0, range1, inclusive=True):
    """
    Tests whether two ranges of integers overlap,
    -2 = range0 inside range1
    -1 = the end of range1 overlaps the start of range0
    0 = no overlap
    1 = the end of range0 overlaps the start of range1
    2 = range1 inside range0
    """
    range0 = [min(range0), max(range0)]
    range1 = [min(range1), max(range1)]
    if inclusive:
        range0 = [range0[0] - 1, range0[1] + 1]
        range1 = [range1[0] - 1, range1[1] + 1]
    if min(range0) > min(range1) and min(range0) < max(range1):
        return 2 if max(range0) < max(range1) else -1
    elif max(range0) > min(range1) and max(range0) < max(range1):
        return -2 if min(range0) > min(range1) else 1
    return 0


class BlastHit(object):
    """Blast hit record"""
    def __init__(self, queryseq=None, subjectseq=None, pident=None,
                 length=None, mismatches=None, gaps=None,
                 query_start=None, query_end=None, subject_start=None,
                 subject_end=None, evalue=None, bitscore=None):
        self.queryseq = queryseq
        self.subjectseq = subjectseq
        self.pident = float(pident) / 100
        self.length = int(length)
        self.mismatches = int(mismatches)
        self.gaps = int(gaps)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.subject_start = int(subject_start)
        self.subject_end = int(subject_end)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)

    def orientation(self):
        if ((self.query_start > self.query_end
             and self.subject_end > self.subject_start) or
            (self.query_end > self.query_start and
             self.subject_start > self.subject_end)):
            return -1
        return 1


class BlastSet(object):
    """Set of blast hits"""
    def __init__(self, hits=None):
        self.hits = [] if hits is None else hits

    def __len__(self):
        return len(self.hits)

    def __getitem__(self, k):
        return self.hits[k]

    def __iter__(self):
        return iter(self.hits)

    def add_tabular(self, arr, ncol=12):
        """Reads a standard 12-column blast input as list"""
        print(arr)
        try:
            xstart = 0 if ncol == 12 else 2
            self.hits.append(BlastHit(
                queryseq=arr[0],
                subjectseq=arr[1],
                pident=arr[2 + xstart],
                length=arr[3 + xstart],
                mismatches=arr[4 + xstart],
                gaps=arr[5 + xstart],
                query_start=arr[6 + xstart],
                query_end=arr[7 + xstart],
                subject_start=arr[8 + xstart],
                subject_end=arr[9 + xstart],
                evalue=arr[10 + xstart],
                bitscore=arr[11 + xstart]))
        except:
            pass
        return ''

    def parse_blast(self, blastfile, ncol=12):
        """Parse 12-column tabular BLAST file"""
        with open(blastfile, 'r') as blastfile:
            for line in blastfile:
                if line[0] == '#':
                    continue
                if len(line) < 10:
                    continue
                self.add_tabular(line.rstrip().split(), ncol=ncol)
        return ''

    def remove_low_pident(self, cutoff, verbose=True):
        """Remove hits with low_pident"""
        if cutoff >= 0 and cutoff <= 1:
            remlist = []
            for i in range(len(self.hits)):
                if self.hits[i].pident < cutoff:
                    remlist.append(i)
        remlist.sort(reverse=True)
        for i in remlist:
            if verbose:
                print("{} removed due to low pid").format(
                    self.hits[i].subjectseq)
            self.hits.pop(i)
        return ''

    def remove_superceded(self, verbose=False):
        """Remove hits where the another
           sequence completely contains the region"""
        remlist = []
        for i in range(len(self.hits) - 1):
            if i in remlist:
                continue
            for j in range(i + 1, len(self.hits)):
                if overlap_range([self.hits[i].query_start,
                                  self.hits[i].query_end],
                                 [self.hits[j].query_start,
                                  self.hits[j].query_end]) == 2:
                    remlist.append(i)
                    break
                if overlap_range([self.hits[j].query_start,
                                  self.hits[j].query_end],
                                 [self.hits[i].query_start,
                                  self.hits[i].query_end]) == 2:
                    remlist.append(j)
        remlist.sort(reverse=True)
        for i in remlist:
            if verbose:
                print("{} superceded".format(self.hits[i].subjectseq))
            self.hits.pop(i)
        return ''

    def n_subjectseq(self):
        """Return number of unique database sequences"""
        return len(set([hit.subjectseq for hit in self.hits]))

    def get_subjectseqs(self):
        """Returns a list of database sequence names"""
        retlist = list(set([hit.subjectseq for hit in self.hits]))
        retlist.sort()
        return retlist

    def get_hits(self, subject='', query=''):
        """Search for specific hits by database name or query name"""
        retlist = []
        for hit in self.hits:
            if ((not subject or subject in hit.subjectseq)
                    and (not query or query in hit.queryseq)):
                retlist.append(hit)
        return retlist
