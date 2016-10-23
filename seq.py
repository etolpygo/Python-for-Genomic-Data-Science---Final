#!/usr/local/bin/python3

import sys
import getopt
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 0.01

seq = 'TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG'


result_handle = NCBIWWW.qblast("blastn", "nt", seq)
# print(result_handle)
blast_record = NCBIXML.read(result_handle)
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print('**alignment**')
            print('sequence: ', alignment.title)
            print('length: ', alignment.length)
            print('e value: ', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)