#!/usr/bin/env python3
'''
ToxcodanCleaner - Script to clean Toxin Families output by ToxCodAn and eliminate putative false-positives
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

##modules
import os
import sys
import subprocess
import datetime as dt
from optparse import OptionParser
from collections import Counter
try:
    from Bio import SeqIO
except:
    print('''ToxcodanCleaner was not able to run due to the ERROR below:
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()


##functions
def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq).upper()
        final[ID] = SEQ
    return final

def _CalculateIdentity_(seqA, seqB):
    sa, sb, sl = seqA, seqB, len(seqA)
    matches = [sa[i] == sb[i] for i in range(sl)]
    seqID = (100 * sum(matches)) / sl
    #gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
    #gapID = (100 * sum(matches)) / gapless_sl
    return seqID

def _AlignSequences_(seqA, seqB):
    from Bio import pairwise2
    alns = pairwise2.align.globalxx(seqA, seqB)
    best_aln = alns[0]
    alignedA, alignedB, score, begin, end = best_aln
    seqID = _CalculateIdentity_(alignedA, alignedB)
    return seqID

def _CLUSTER_(F, cluster):
    if cluster <= 1.0:
        cluster = cluster*100

    S = {}
    for k in F.keys():
        SEQ = F[k]
        S.setdefault(len(SEQ),[])
        S[len(SEQ)].append(k)

    WRITE = {}
    REPORT = {}
    for size in S.keys():
        if len(S[size]) != 1:
            CLUSTERED = {}
            for n in range(len(S[size])):
                CLUSTERED.setdefault(S[size][n], [])
                for m in range(n+1, len(S[size])):
                    seqREF = F[S[size][n]]
                    seqNEX = F[S[size][m]]
                    seqID = _AlignSequences_(seqREF, seqNEX)
                    if seqID > cluster:
                        CLUSTERED[S[size][n]].append([S[size][m], seqID])
            TOREMOVE = []
            for j in CLUSTERED.keys():
                for c in CLUSTERED[j]:
                    TOREMOVE.append(c[0])
            for j in TOREMOVE:
                CLUSTERED.pop(j, None)
            for j in CLUSTERED.keys():
                REPORT[j] = CLUSTERED[j]
                WRITE[j] = F[j]
        if len(S[size]) == 1:
            WRITE[S[size][0]] = F[S[size][0]]
            REPORT[S[size][0]] = []

    return WRITE, REPORT

def _ToxcodanCleaner_(cds, output, minS, maxS, partial, cluster, cpu):
    REPORT = open(output.split(".")[0]+"_REPORT.txt", "w")
    F = _ParseFasta_(cds)
    #filter step (size and partial)
    OUTfilt = open(cds.split(".")[0]+"_FILTERED.fasta","w")
    for k in F.keys():
        if partial in ["True", "TRUE"]:
            if len(F[k]) % 3 == 0 and F[k][-3:] in ["TAA", "TAG", "TGA"]:
                if minS <= len(F[k]) <= maxS:
                    OUTfilt.write(">"+k+"\n"+F[k]+"\n")
                    continue
                else:
                    REPORT.write(k+" -> filtered [size filter]\n")
            if len(F[k]) % 3 == 0 and F[k][-3:] not in ["TAA", "TAG", "TGA"]:
                REPORT.write(k+" -> filtered [partial filter]\n")
            if len(F[k]) % 3 != 0:
                REPORT.write(k+" -> filtered [partial filter]\n")
        if partial in ["False", "FALSE"]:
            if minS <= len(F[k]) <= maxS:
                OUTfilt.write(">"+k+"\n"+F[k]+"\n")
                continue
            else:
                REPORT.write(k+" -> filtered [size filter]\n")
    OUTfilt.close()
    #alignment step
    subprocess.call("mafft --thread -"+cpu+" --auto "+cds.split(".")[0]+"_FILTERED.fasta > "+cds.split(".")[0]+"_ALIGNED.fasta",shell=True)
    #clean data based on msa
    alignment = cds.split(".")[0]+"_ALIGNED.fasta"
    A = _ParseFasta_(alignment)
    OUTclean = open(output,"w")
    correctATG = []
    for k in A.keys():
        correctATG.append(A[k].find("ATG"))
    mainATG = Counter(correctATG).most_common(1)[0][0]
    RemoveRedundancy = []
    for k in A.keys():
        if A[k].find("ATG") > mainATG:
            REPORT.write(k+" -> removed; ATG in the middle of a correct CDS\n")
        if A[k].find("ATG") < mainATG:
            REPORT.write(k+" -> trimmed; upstream ATG used and it is not correct")
            if A[k][mainATG:].replace("-","") in RemoveRedundancy:
                REPORT.write("; but it turned into 100% identical to another CDS, removed to reduce redundancy\n")
            if A[k][mainATG:].replace("-","") not in RemoveRedundancy:
                OUTclean.write(">"+k+"\n"+A[k][mainATG:].replace("-","")+"\n")
                RemoveRedundancy.append(A[k][mainATG:].replace("-",""))
                REPORT.write("\n")
        if A[k].find("ATG") == mainATG:
            if A[k].replace("-","") in RemoveRedundancy:
                REPORT.write(k+" -> kept; but it is 100% identical to another CDS, removed to reduce redundancy\n")
            if A[k].replace("-","") not in RemoveRedundancy:
                OUTclean.write(">"+k+"\n"+A[k].replace("-","")+"\n")
                RemoveRedundancy.append(A[k].replace("-",""))
                REPORT.write(k+" -> kept\n")

    OUTclean.close()
    REPORT.close()

    #clustering step
    F = _ParseFasta_(output)
    if cluster not in ["False", "FALSE"]:
        write, report = _CLUSTER_(F, cluster)
        OUTclst = open(output.split(".")[0]+"_CLST.fasta","w")
        for k in write.keys():
            OUTclst.write(">"+k+"\n"+write[k]+"\n")
        OUTclst.close()
        REPORT = open(output.split(".")[0]+"_CLST_REPORT.txt","w")
        n = 1
        for k in report.keys():
            REPORT.write(">Cluster_"+str(n)+"\n")
            REPORT.write("\t"+k+" -> *\n")
            n += 1
            if report[k] != []:
                for i in report[k]:
                    REPORT.write("\t"+i[0]+" -> "+str(i[1])+"\n")
        REPORT.close()


##options
def __main__():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="cds", help="Mandatory - CDSs of a toxin family output by ToxCodAn in FASTA format, /path/to/TOXIN_cds.fasta", metavar="fasta", default=None)
    parser.add_option("-m", "--minsize", dest="minS", help="Optional - threshold value used as the minimum size of the CDSs of the toxin family to be filter [default=200]", metavar="int", default="200")
    parser.add_option("-M", "--maxsize", dest="maxS", help="Optional - threshold value used as the maximum size of the CDSs of the toxin family to be filter [default=4000]", metavar="int", default="4000")
    parser.add_option("-p", "--partial", dest="partial", help="Optional - turn on/off the partial filtering step, use True to turn on or False to turn off [default=False]", metavar="boolean value", default="False")
    parser.add_option("-C", "--cluster", dest="cluster", help="Optional - turn on/off the cluster step, set any threshold [e.g., 0.99 to cluster 99% similar sequences with similar size] to turn on [default=False]", metavar="float", default="False")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used [default=1]", metavar="int", default="1")

    (options, args) = parser.parse_args()

    if options.cds == None:
        print(
        """

>>>> ToxcodanCleaner v1.0 May 2022 <<<<
      ****Use -h for help!****

USAGE:
ToxcodanCleaner.py -i TOXIN_cds.fasta
        """)
        quit()

    output = str(options.cds).split(".")[0]+"_CLEAN.fasta"

    if options.cds != None:
        #print("""""")
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting ToxcodanCleaner (v1.0 May 2022)...")
        CWD = os.getcwd()
        print("\ttoxin CDS file ->", options.cds)
        print("\tminimum size filter ->", options.minS)
        print("\tmaximum size filter -> ", options.maxS)
        print("\tpartial filter -> ", options.partial)
        print("\tcluster threshold -> ", options.cluster)
        print("\tNumber of threads ->", options.cpu)

        _ToxcodanCleaner_(options.cds,
                    output,
                    int(options.minS),
                    int(options.maxS),
                    str(options.partial),
                    float(options.cluster),
                    str(options.cpu))

        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> ToxcodanCleaner finished!")

        print('''
        Description of the output files:
            _CLEAN.fasta -> final file with the false-positive sequences removed and trimmed
            _CLEAN_REPORT.txt -> report about the processing of each inputted sequence
            _FILTERED.fasta -> file with the filtered CDSs [partial and size filter steps] used to align sequences to be cleaned [trim sequences and remove spurious sequences]
            _ALIGNED.fasta -> file with the alignment used to clean the toxin CDSs
            _CLEAN_CLST.fasta -> file with clustered CDSs (if cluster step was turned on)
            _CLEAN_CLST_REPORT.txt -> report about clustering step output (if cluster was turned on)
        ''')

if __name__ == '__main__':
    __main__()

#END
