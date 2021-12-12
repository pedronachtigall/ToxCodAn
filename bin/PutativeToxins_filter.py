#!/usr/bin/env python3
'''
parse blast result in the PutativeToxin filter step of venomancer
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''
import sys

import os
from Bio import SeqIO
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = Nonefrom collections import defaultdict

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _ParseBlast_(blast):
    final = {}
    rejected = set([])
    a = open(blast,"r")
    for line in a:
        line1 = line.rstrip().split("\t")
        ID = line1[0]
        sname = line1[2]
        slen = int(line1[3])
        align = int(line1[5])
        sst = int(line1[9])
        sen = int(line1[10])
        evalue = float(line1[11])
        score = float(line1[12])
        final[ID] = [ID,sname,evalue,score]
    return final

def _GetOutput_(blast, fasta, folder, prefix):
    BLAST = _ParseBlast_(blast)
    FASTA = _ParseFasta_(fasta)

    ADDED = []

    PTcds = open(folder+prefix+"_PutativeToxins_cds_SPfiltered.fasta","w")
    PTcontigs = open(folder+prefix+"_PutativeToxins_contigs_SPfiltered.fasta","w")
    PTpep = open(folder+prefix+"_PutativeToxins_pep_SPfiltered.fasta","w")
    for k in BLAST.keys():
        for j in FASTA.keys():
            if k in j and k not in ADDED:
                j1 = j.split("||_")
                j2 = j1[-1].split("_")
                st = int(j2[0])
                end = int(j2[1])
                SEQ = FASTA[j]
                hit = BLAST[k][1]
                if end+3 <= len(SEQ):
                    end = end+3
                #write contigs
                seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                PTcontigs.write(">"+k+"||"+hit+"\n"+seq+"\n")
                #write cds
                CDS = SEQ[st-1:end]
                #print(j, CDS)
                cds = "\n".join([CDS[n:n+100] for n in range(0, len(CDS), 100)])
                PTcds.write(">"+k+"||"+hit+"\n"+cds+"\n")
                #translate CDSs
                if generic_dna:
                    CodingDna = Seq(SEQ, generic_dna)
                else:
                    CodingDna = Seq(SEQ)
                pep = str(CodingDna.translate(table=1))
                PTpep.write(">"+k+"||"+hit+"\n"+pep+"\n")
                ADDED.append(k)
    PTcds.close()
    PTcontigs.close()
    PTpep.close()

def _main_():
    if len (sys.argv) != 5:
        print("Basic usage: toxin_filter.py blast.out transcripts.fa output_folder_prefix")
        print("\t> blast.out: blast output file in tab-delimited format")
        print("\t> transcript.fa: PutativeToxin transcripts")
        print("\t> output_folder: folder to output results")
        print("\t> prefix: output prefix to be used in files")
        quit()

    blast = sys.argv[1]
    fasta = sys.argv[2]
    folder = sys.argv[3]
    prefix = sys.argv[4]
    _GetOutput_(blast, fasta, folder, prefix)

_main_()

#END
