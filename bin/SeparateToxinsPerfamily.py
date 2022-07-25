#!/usr/bin/env python3
'''
SeparateToxinsPerfamily - Script to generate specific toxin files from ToxCodAn's output
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import sys
import os
from Bio import SeqIO

def _GetToxins_(fasta, folder):
    final = {}
    F = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        if "_PARTIAL" in str(record.id):
            toxin = ID.split("_")[-2].replace(":","_")
        if "_PARTIAL" not in str(record.id):
            toxin = ID.split("_")[-1].replace(":","_")
        if "SVMP" in toxin:
            toxin = "SVMP"
        if "SNACLEC" in toxin:
            toxin = "CTL"
        if "LAAO" in toxin:
            toxin = "LAO"
        if "HOM_Toxin" in toxin:
            toxin = "HOM"
        final.setdefault(toxin, [])
        final[toxin].append(str(record.id))
        F[ID] = str(record.seq)
    for toxin in final.keys():
        fileo = fasta.split("/")[-1]
        OUT = open(folder+toxin+"_"+fileo, "w")
        for id in final[toxin]:
            OUT.write(">"+id+"\n"+F[id]+"\n")
        OUT.close()

def __main__():

    if len (sys.argv) != 3:
        print("Script designed to generate specific toxin files from ToxCodAn's output")
        print("Basic usage: SeparateToxinsPerfamily.py Toxins.fa output_folder")
        print("\t> Toxins.fa: input file in fasta format [Toxins and/or PutativeToxins output by ToxCodAn]")
        print("\t> output_folder: folder to output the separated files")
        quit()

    fasta = sys.argv[1]
    folder = sys.argv[2]

    if not folder.endswith("/"):
        folder += "/"
    if os.path.isdir(folder) == False:
        os.mkdir(folder)

    _GetToxins_(fasta, folder)

if __name__ == '__main__':
    __main__()

#END
