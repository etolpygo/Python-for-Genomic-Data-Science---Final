#!/usr/local/bin/python3
def usage():
    print("""
builds a dictionary with all sequences from a FASTA file, whose name is passed in as an argument.

analyze.py [-h] [-l <length>] <filename>

-h              print this message
-l <length>     filter out sequences with length smaller than <length>
<filename>      the file in FASTA format
""")

import sys
import getopt
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

o, a = getopt.getopt(sys.argv[1:], 'l:h')
opts = {}
seqlen = 0
E_VALUE_THRESH = 0.01

# for k,v in o:
#     opts[k] = v
# if 'h' in opts.keys():
#     usage(); exit()
# if len(a) < 1:
#     usage(); exit("input file is missing.")
# if 'l' in opts.keys():
#     if opts['l'] < 0:
#         exit("length of sequence should be positive.");
#     seqlen = opts['l']

fname = a[0]

start_codon = 'ATG'
stop_codons = ['TAA', 'TAG', 'TGA']

try:
    f = open(fname)
except IOError:
    print("File %s does not exist." %fname)
    exit()

seqs = {}
for line in f:
    line = line.rstrip()
    if line.startswith('>'):
        words = line.split()
        name = words[0][1:]
        seqs[name] = ''
    else:
        seqs[name] = seqs[name] + line

f.close()

print("there are ", len(seqs), " sequences.")

longest = ''
for name,seq in seqs.items():
    if len(seq) > len(longest):
        longest = seq
print("length of the longest one is ", len(longest))        


shortest = longest
for name,seq in seqs.items():
    if len(seq) < len(shortest):
        shortest = seq
print("length of the shortest one is ", len(shortest))        


ORFs = {}
frame = 2
for name,seq in seqs.items():
    started = False
    started_pos = None
    stopped_pos = None
    for i in range(frame, len(seq), 3):
        codon = seq[i:i+3].upper()
        if started == False and codon == start_codon:
            started = True
            started_pos = i
        if started == True and codon in stop_codons:
            stopped_pos = i+3
            ORF = seq[started_pos:stopped_pos]
            ORFs[(name,started_pos)] = ORF
            started = False
            started_pos = None
            stopped_pos = None

longest_ORF = ''
longest_ORF_pos = 0
for (name,i),ORF in ORFs.items():
    if len(ORF) > len(longest_ORF) and name == 'gi|142022655|gb|EQ086233.1|16':
        longest_ORF = ORF
        longest_ORF_pos = i
        
print("longest ORF in frame ", frame+1, " is: ", longest_ORF)
print("its length is: ", len(longest_ORF))
print("its starting position is: ", longest_ORF_pos)

rep_length = 12
repeats = {}
for name,seq in seqs.items():
    for i in range(0, len(seq)):
        string = seq[i:i+rep_length].upper()
        if string in repeats:
            repeats[string] = repeats[string] + 1
        else:
            repeats[string] = 1

most_repeats = 0
for (repeat,times) in repeats.items():
    if times > most_repeats:
        most_repeats = times
        
print("most common repeat of ", rep_length, " length occurs ", most_repeats, "times")    

longest_repeats = []  
for (repeat,times) in repeats.items():
    if times == most_repeats:
        longest_repeats.append(repeat)
        
print(longest_repeats)   


lst = ['CGCGCCG', 'TGCGCGC', 'CATCGCC', 'GCGCGCA']
counts = {}
for name,seq in seqs.items():
    for sq in lst:
        count_seq = seq.count(sq)    
        if sq in counts:
            counts[sq] = counts[sq]+count_seq
        else:
            counts[sq] = count_seq
            
            
print(counts)


# for name,seq in seqs.items():
#     # print(name,seq)
#     result_handle = NCBIWWW.qblast("blastn", "nt", seq)
#     # print(result_handle)
#     blast_record = NCBIXML.read(result_handle)
#     for alignment in blast_record.alignments:
#         for hsp in alignment.hsps:
#             if hsp.expect < E_VALUE_THRESH:
#                 print('**alignment**')
#                 print('sequence: ', alignment.title)
#                 print('length: ', alignment.length)
#                 print('e value: ', hsp.expect)
#                 print(hsp.query)
#                 print(hsp.match)
#                 print(hsp.sbjct)