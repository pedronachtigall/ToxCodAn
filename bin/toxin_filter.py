#!/usr/bin/env python3
'''
parse blast result in the toxin filter step of venomancer
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''
import sys

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _ParseGTF_(gtf, gene):
    final = []
    a = open(gtf, "r")
    for line in a:
        line1 = line.split("\t")
        if line1[0] == gene:
            final.append(line1)
    a.close()
    return final

def _ParseBlast_(blast, model):
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
        if align > slen*0.8 and sen-sst > slen*0.8 and (evalue <= 0.0 or "e-" in str(evalue)):
            final.setdefault(ID,[])
            if len(final[ID]) > 0:
                SCORE = final[ID][3]
                if score > SCORE:
                    final[ID] = [ID,sname,evalue,score,model]
            else:
                final[ID] = [ID,sname,evalue,score,model]
        else:
            rejected.add(ID)
    a.close()
    return final

def _GetOutput_(folder, fasta, output):
    M = ["TOXIN","TOXIN_ACHE","TOXIN_HYAL","TOXIN_MYO","TOXIN_PDE","TOXIN_VEGF",
        "TOXIN_3FTx","TOXIN_CTL","TOXIN_KUN","TOXIN_NUC","TOXIN_PLA2"]
    dd = defaultdict(list)
    for m in M:
        blast = folder+"BLAST_"+m+".out"
        a = _ParseBlast_(blast, m)
        for d in (dd, a):
            for key, value in d.items():
                dd[key].append(value)

    final = {}
    for k in sorted(dd.keys()):
        for i in dd[k]:
            if i[0] == k:
                final.setdefault(k,[])
                if len(final[k]) == 0:
                    final[k] = i
                else:
                    alfa = final[k]
                    if i[3] > alfa[3]:
                        final[k] = i

    #generating output files
    #toxins
    cds = open(folder+output+"_Toxins_cds.fasta","w")
    gtf = open(folder+output+"_Toxins_annotation.gtf","w")
    ToxinsPEP = open(folder+output+"_Toxins_pep.fasta","w")
    for k in sorted(final.keys()):
        alfa = final[k]
        hit = alfa[1]
        m = alfa[-1]
        FASTA = _ParseFasta_(folder+m+"_out/ORF_sequences.fasta")
        SEQ = FASTA[k]
        cdsSeq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        cds.write(">"+k+"||"+hit+"\n"+cdsSeq+"\n")
        #translate toxins CDSs
        CodingDna = Seq(SEQ, generic_dna)
        pep = str(CodingDna.translate(table=1))
        ToxinsPEP.write(">"+k+"||"+hit+"\n"+pep+"\n")
        #GTF
        GTF = _ParseGTF_(folder+m+"_out/annotation.gtf", k)
        for i in GTF:
            gtf.write("\t".join(i))
    cds.close()
    gtf.close()
    ToxinsPEP.close()

    #putative toxins
    cds = open(folder+output+"_PutativeToxins_cds.fasta","w")
    for m in M:
        FASTA = _ParseFasta_(folder+m+"_out/ORF_sequences.fasta")
        for k in sorted(FASTA.keys()):
            if k not in final.keys():
                SEQ = FASTA[k]
                cdsSeq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                cds.write(">"+k+"||"+m+"\n"+cdsSeq+"\n")
    cds.close()

def _main_():
    if len (sys.argv) != 4:
        print("Basic usage: toxin_filter.py folder transcripts.fa output_prefix")
        print("\t> folder: folder with all blast output files in tab-delimited format")
        print("\t> transcript.fa: pre-filtered toxin transcripts")
        print("\t> output_prefix: output prefix to be used in files")
        quit()

    folder = sys.argv[1]
    fasta = sys.argv[2]
    output = sys.argv[3]
    _GetOutput_(folder, fasta, output)

_main_()

#END
