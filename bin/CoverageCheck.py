#!/usr/bin/env python

# Additional software necessary to run this:
# (1) bwa
# (2) samtools
# (3) bedtools
# (4) biopython
# (5) pandas
# (6) dfply

import argparse
import sys, os, shutil
import subprocess as sp
import datetime as dt
import csv
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from dfply import *

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='Check read coverage within transcripts to remove transcripts with low coverage across a certain percentage of the contig/transcript and transcripts containing only multi-mapped reads. This may be indicitive of a misassembly or false-positive annotation. Default is to check for transcripts with <5x coverage for >10% of the total contig/transcript length and remove sequences with 0.0 TPM/FPKM expression level.')
parser.add_argument("-i","--input",
					type=str,
					help="Fasta file of contigs/transcripts to check. Use only CODING SEQUENCES.")
parser.add_argument("-r","--reads",
					type=str,
					help="Fastq file of UNPAIRED reads.")
parser.add_argument("-c","--cov",
					type=int,
					default=5,
					help="Coverage threshold. (default: %(default)s)")
parser.add_argument("-p","--prop",
					type=int,
					default=10,
					help="Proportion of the contig/transcript. (default: %(default)s)")
parser.add_argument("-e","--exp",
					type=float,
					default=0.0,
					help="Expression value (in TPM) to be used as threshold to filter low expressed contigs/transcripts, which may be a result of only multi-mapped reads in the contig/transcript. (default: %(default)s)")
parser.add_argument("-m","--mr",
					type=float,
					default=0.1,
					help="The maximum mismatch rate allowed to be set in RSEM parameter \"--bowtie2-mismatch-rate\" (to be used in Bowtie2 mapping), which by default is 0.1 (as specified by RSEM manual). (default: %(default)s)")
parser.add_argument("-t","--threads",
					type=int,
					default=1,
					help="Number of processing threads. (default: %(default)s)")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

input = os.path.abspath(args.input)
input_name = os.path.basename(input)
input_name2 = input_name.split(".fasta")[0]
reads = os.path.abspath(args.reads)
cov = args.cov
prop = args.prop
exp = args.exp
mr = args.mr
threads = args.threads

print("""

CoverageCheck - filtering spurious sequences based on read coverage

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
print("\tInput -> "+ input_name)
print("\tReads -> "+ os.path.basename(reads))
print("\tRemoving transcripts with...")
print("\t\t< " + str(cov) + "x coverage")
print("\t\tfor > " + str(prop) + "% of the total transcript length")
print("\tExpression level threshold (RSEM) --> "+ str(exp))
print("\tMismatch rate (RSEM) -> "+ str(mr))
print("\tThreads -> " + str(threads))

########################################
################# CODE #################
########################################

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Generating bwa alignment :::")
sp.call("bwa index " + input, shell=True)
sp.call("bwa mem -t " + str(threads) + " " + input + " " + reads + " | samtools sort -@ " + str(threads) + " > " + input_name2 + ".bam", shell=True)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Calculating coverage :::")
sp.call("bedtools genomecov -d -ibam " + input_name2 + ".bam > coverage.txt", shell=True)
sp.call("rm " + input + ".* " + input_name2 + ".bam", shell=True)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Assessing coverage for Presence/Absence :::")
coverage = pd.read_csv("coverage.txt",sep='\t',names=['transcript','pos','cov'])
seqs = list(SeqIO.parse(input,'fasta'))
print("\t>>> Input transcripts = " + str(len(seqs)))
print("\t>>> Coverage data for = " + str(len(set(coverage['transcript']))))

tmp1 = coverage >>  group_by(X.transcript) >> summarize(length=n(X.pos))
tmp2 = coverage[coverage['cov'] < cov]
if(len(tmp2)!=0):
	tmp2 = tmp2 >> group_by(X.transcript) >> summarize(lowcov_count=n(X.pos))
	tmp3 = tmp2 >> right_join(tmp1, by='transcript')
	tmp4 = tmp3 >> mutate(lowcov_prop = (X.lowcov_count/X.length)*100)
	tmp4.to_csv("coverage_low.txt", sep="\t", index=False)
	tmp5 = tmp4[tmp4['lowcov_prop'] > prop]
	T = tmp5['transcript'].tolist()
else:
	T = list()

print(str(len(T)) + " transcripts with less than " + str(cov) + "x coverage for greater than " + str(prop) + "% of total transcript length")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Removing \"Absent\" sequences :::")

#check for any sites with no coverage in the CDS (putative chimeras)
zeroBad = []
for seq in seqs:
	zeros = list(coverage[coverage['transcript'] == seq.name]['cov']).count(0)
	if zeros :
		zeroBad.append(seq.name)

#print(len(zeroBad),"-> transcripts with zero coverage in one position")

good_seqs = []
for seq in seqs:
	if seq.name not in T and seq.name not in zeroBad:
		good_seqs.append(seq)

print("\t>>>",len(good_seqs),"transcripts with good coverage")

output_handle=open(input_name2 + "_GOODCOVERAGE.fasta", "w")
writer = FastaIO.FastaWriter(output_handle, wrap=None)
writer.write_file(good_seqs)
output_handle.close()

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Removing sequences with TPM <= "+str(exp)+" :::")

#run RSEM
sp.call("rsem-prepare-reference --bowtie2 " + input_name2 + "_GOODCOVERAGE.fasta rsem_reference", shell=True)
sp.call("rsem-calculate-expression -p " + str(threads) + " --bowtie2 --bowtie2-mismatch-rate " + str(mr) + " " + reads + " rsem_reference rsem", shell=True)

# Remove the contigs with low TPM sites from the temp list
dataFrameRSEM = pd.read_csv("rsem.isoforms.results", sep="\t")
zeroTPM = set(list(dataFrameRSEM[dataFrameRSEM["TPM"] <= float(exp)]["transcript_id"]))
goodrsem = []
for seq in good_seqs:
	if seq.name not in zeroTPM:
		goodrsem.append(seq)

print("\t>>>",len(goodrsem),"transcripts with TPM value higher than "+str(exp))

output_handle=open(input_name2 + "_FILTERED.fasta", "w")
writer = FastaIO.FastaWriter(output_handle, wrap=None)
writer.write_file(goodrsem)
output_handle.close()

#generating report about unique and multi-mapped reads
print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Generating report about mapped reads :::")

sp.call("samtools view -@ " + str(threads) + " -h -o rsem.sam rsem.transcript.bam",shell=True)

def _CheckReads_(sam, out):
	countUNMAPPED = 0
	countUR = 0
	countMR = 0
	ReadMR = {}
	RefM = {}
	primary = {}
	secondary = {}
	F = []
	a = open(sam,"r")
	for line in a:
		if line.startswith("@SQ"):
			header = line.strip().split()[1].replace("SN:","")
			F.append(header)
		if not line.startswith("@"):
			line1 = line.strip().split("\t")
			read = line1[0]
			flag = line1[1]
			ref = line1[2]
			if flag != "4": #mapped reads - 4 indicates an unmapped read
				ReadMR.setdefault(read, set([]))
				ReadMR[read].add(ref)
				if flag == "0" or flag == "16": #primary alignments (independently of multimapping status)
					primary.setdefault(read, set([]))
					primary[read].add(ref)
				if flag == "256" or flag == "272": #secondary alignments of multi-mapped reads
					secondary.setdefault(read, set([]))
					secondary[read].add(ref)
			if flag == "4":
				countUNMAPPED += 1
	a.close()
	for k in ReadMR.keys():
		if len(ReadMR[k]) > 1:
			countMR += 1
			for ref in ReadMR[k]: #multi-mapped
				RefM.setdefault(ref, [0, 0, 0, 0])
				RefM[ref][1] += 1
			if k in primary.keys():
				for ref in primary[k]:
					RefM.setdefault(ref, [0, 0, 0, 0])
					RefM[ref][2] += 1
			if k in secondary.keys():
				for ref in secondary[k]:
					RefM.setdefault(ref, [0, 0, 0, 0])
					RefM[ref][3] += 1
		if len(ReadMR[k]) == 1: #unique mapped
			countUR += 1
			ref = list(ReadMR[k])[0]
			RefM.setdefault(ref, [0, 0, 0, 0])
			RefM[ref][0] += 1

	OUT = open(out,"w")
	OUT.write("id\tUniqueMapped\tMultiMapped\tprimary\tsecondary\n")
	for k in F:
		if k in RefM.keys():
			OUT.write(k+"\t"+"\t".join([str(x) for x in RefM[k]])+"\n")
		if k not in RefM.keys():
			OUT.write(k+"\t0\t0\n")
	OUT.close()
	print("\tReport about reads:")
	print("\tTotal Reads ->",countUNMAPPED+countUR+countMR)
	print("\t\tUnMapped ->",countUNMAPPED)
	print("\t\tMapped ->",countUR+countMR)
	print("\t\t\t- UniqueMapped ->",countUR)
	print("\t\t\t- MultiMapped -> ",countMR)

_CheckReads_("rsem.sam", "ReadsReport.txt")

sp.call("rm rsem.sam",shell=True)
sp.call("rm rsem_reference.*",shell=True)
sp.call("rm rsem.transcript.b*",shell=True)
sp.call("rm coverage.txt",shell=True)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Finished! :::")

#END
