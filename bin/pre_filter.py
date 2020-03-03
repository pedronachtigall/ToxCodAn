#!/usr/bin/env python3
'''
parse blast result in the pre-filter step of venomancer
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''
import sys
from Bio import SeqIO

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

#tblastn -query ../blastDB/FinalToxinDB.fa -out pre_filter_blast.out -db blastDB/ChorrB_contigs -num_threads 6 -qcov_hsp_perc 70 -max_target_seqs 10 -outfmt '6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore'
def _ParseBlast_(blast):
    final = set([])
    a = open(blast,"r")
    for line in a:
        line1 = line.rstrip().split("\t")
        final.add(line1[2])
    a.close()
    return list(final)

def _GetOutput_(blast, fasta, output):
    BLAST = _ParseBlast_(blast)
    FASTA = _ParseFasta_(fasta)
    OUT = open(output,"w")
    for i in BLAST:
        SEQ = "\n".join([FASTA[i][n:n+100] for n in range(0, len(FASTA[i]), 100)])
        OUT.write(">"+i+"\n"+SEQ+"\n")
    OUT.close()

def _main_():
    if len(sys.argv) != 4:
        print("Basic usage: pre_filter.py blast_out transcripts.fa out.fa")
        print("\t> blast_out: blast output in tab-delimited format")
        print("\t> transcript.fa: input transcripts in fasta format")
        print("\t> out.fa: output filtered sequences")
        quit()

    blast = sys.argv[1]
    fasta = sys.argv[2]
    output = sys.argv[3]
    _GetOutput_(blast, fasta, output)

_main_()

#END
