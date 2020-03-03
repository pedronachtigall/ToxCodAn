#!/usr/bin/env python3
'''
parse signalP result to filter transcripts with no signalP predicted
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''
import sys
import os
from Bio import SeqIO

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _ParseSignalP_(signalP):
    final = []
    s = open(signalP,"r")
    for line in s:
        if not line.startswith("#"):
            line1 = line.rstrip().split("\t")
            id = line1[0]
            final.append(id)
    s.close()
    return final

def _SignalPfilter_(signalP, prefix):
    CDS = _ParseFasta_(prefix+"_Toxins_cds.fasta")
    CONTIGS = _ParseFasta_(prefix+"_Toxins_contigs.fasta")
    PEP = _ParseFasta_(prefix+"_Toxins_pep.fasta")
    SIGNAL = _ParseSignalP_(signalP)
    REDUNDANCY = _ParseFasta_(prefix+"_Toxins_cds_RedundancyFiltered.fasta")

    OUTcds = open(prefix+"_Toxins_cds_SPfiltered.fasta","w")
    OUTcontigs = open(prefix+"_Toxins_contigs_SPfiltered.fasta","w")
    OUTpep = open(prefix+"_Toxins_pep_SPfiltered.fasta","w")
    OUTredundancy = open(prefix+"_Toxins_cds_SPfiltered_RedundancyFiltered.fasta","w")
    for i in SIGNAL:
        SEQ = CDS[i]
        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        OUTcds.write(">"+i+"\n"+seq+"\n")
        i1 = i.split("||")
        id = i1[0]
        SEQ = CONTIGS[id]
        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        OUTcontigs.write(">"+i+"\n"+seq+"\n")
        SEQ = PEP[i]
        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        OUTpep.write(">"+i+"\n"+seq+"\n")
        if i in REDUNDANCY.keys():
            SEQ = REDUNDANCY[i]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            OUTredundancy.write(">"+i+"\n"+seq+"\n")
    OUTcds.close()
    OUTcontigs.close()
    OUTpep.close()
    OUTredundancy.close()

def _main_():
    if len (sys.argv) != 3:
        print("Basic usage: signalp_filter.py signalp.gff folder+prefix")
        print("\t> signalp: signalp result in gff format")
        print("\t> folder+prefix: input and output prefix to be used in files")
        quit()

    signalp = sys.argv[1]
    prefix = sys.argv[2]
    _SignalPfilter_(signalp, prefix)

_main_()

#END
