#!/usr/bin/env python
"""
parse blast result of the PutativeToxins generated by venomancer
Author: Schyler Ellsworth & Pedro G. Nachtigall
"""
import sys
import re
import os
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from random import *
from Bio.Blast import NCBIXML
from Bio import SeqIO

#Reverse transcribes the sequence
def reverse_sequence(s):
	"""Return a reverse sequence of `s`"""
	s = s.replace("G","B")
	s = s.replace("C","G")
	s = s.replace("B","C")
	s = s.replace("T","U")
	s = s.replace("A","T")
	s = s.replace("U","A")
	return(s)

#Translates Sequence
def translate(seq):

	ambiguities=["Y","R","W","S","K","M","N","D","V","H","B"]
	table = {
		"ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
		"ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
		"AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
		"AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
		"CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
		"CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
		"CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
		"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
		"GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
		"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
		"GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
		"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
		"TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
		"TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
		"TAC":"Y", "TAT":"Y", "TAA":"*", "TAG":"*",
		"TGC":"C", "TGT":"C", "TGA":"*", "TGG":"W",
	}
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3):
			codon = str(seq[i:i + 3]).upper()
			if any(base in codon for base in ambiguities):
				protein+="X"
			else:
				protein+= table[codon]
	return protein

def signalp(p, folder):
	f = open(folder+"signalp.fasta","w+")
	f.write(">P"+"\n"+p[:-1])
	f.close()
	sp = "signalp -u 0.34 -U 0.34 "+folder+"signalp.fasta  > "+folder+"signalp.txt"
	os.system(sp)
	cut = 0
	file = open(folder+"signalp.txt", "rU")
	for line in file:
		if line[0] == "#":
			line = line
		else:
			line = line.strip("\n")
			signalre = r"([NY])"
			signal = re.search(signalre,line)
			sig = signal.group(1)
			if sig == "Y":
				line = line.replace(" ","l")
				cut = line[34:36]
	return(sig,cut)

def annotate(a, folder):
	bestre = r"(M\w+)"
	best = re.search(bestre,a)
	#print(a, best)
	if best != None:
		best = best.group(1)
		sig = signalp(best, folder)
	else:
		best = "T"
		sig = signalp(best, folder)
	return(best,sig)

def startpos(P, n):
    x = P.find("M")
    while x >= 0 and n > 1:
        x = P.find("M", x+1)
        n -= 1
    return x

if len(sys.argv) != 4:
	print("""
Usage: BlastNamer.py file.fasta file.xml folder > folder/report.log
	file.fasta - fasta file to have the CDS annotated
	file.xml - xml file with the BLAST search results
	folder - output folder
	""")
else:
	FASTA = sys.argv[1] # reads txt
	xml = sys.argv[2] #reads fasta
	folder = sys.argv[3] #folder to output annotations
	blast = open(xml,"r")
	blast_records = NCBIXML.parse(blast)
	List=[]
	for rec in blast_records:
		for seq_record in SeqIO.parse(FASTA, "fasta"):
			if rec.query.split()[0] == seq_record.id:
				aligncount = 1
				for alignment in rec.alignments:
					if aligncount == 1:
						aligncount +=1
						hspcount = 1
						for hsp in alignment.hsps:
							#print(hsp)
							if hspcount == 1:
								hspcount +=1
								count=0
								if hsp.frame[0] < 0 :
									start = hsp.query_end
									end = hsp.query_start
									Tstart = start
									Tend= end
								else :
									start = hsp.query_start
									end = hsp.query_end
									Tstart = start
									Tend= end
								if start > end:
									start = (len(str(seq_record.seq))-start+1)
									start = (start - (int(start/3)*3))
									if start == 0:
										start = 3
									end = (len(str(seq_record.seq))-end+4)
									seq = str(reverse_sequence(str(seq_record.seq[::-1])))
								else:
									start = (start - (int(start/3)*3))
									seq = seq_record.seq
									if start == 0:
										start = 3
								if int(Tend) < int(Tstart):
									Tstart = (len(str(seq_record.seq))-Tstart+1)
									Tend = (len(str(seq_record.seq))-Tend+4)
								if len(seq[start-1:]) % 3 == 0:
									CDS = seq[start-1:]
								elif len(seq[start-1:-1]) % 3 == 0:
									CDS = seq[start-1:-1]
								elif len(seq[start-1:-2]) % 3 == 0:
									CDS = seq[start-1:-2]
								prot = translate(str(CDS))
								prot2 = prot
								while prot2.count("*") >= 1:
									num = prot2.find("*")
									prot3 = prot2[:num+1]
									if prot3.count("M") == 0:
										prot2 = prot2[num+1:]
									else:
										while prot3.count("M") >=1:
											mcount= 0
											while count == 0:
												count = 0
												#print(prot3)
												if prot3.count("M") == 0:
													break
												if prot3 == "M":
													break
												if prot3.count("M") ==1 and prot3[-2] =="M":
													break
												pep = annotate(prot3, folder)
												prot3 = pep[0]
												sig = pep[1]
												mcount+=1
												mstart = startpos(prot2,mcount)
												new=prot2[mstart:]
												b = (len(prot[:])-len(new))*3
												if sig[0] =="Y" and Tstart <= ((len(prot3)*3)+start+b-1) and Tstart >= (start+b) and Tend <= ((len(prot3)*3)+start+b-1):
													n= open(folder+"PutativeToxins_annotation.fasta","a+")
													n.write(">"+rec.query.split()[0]+"||_"+str(start+b)+"_"+str((len(prot3)*3)+start+b-1)+"_sp"+str(int(sig[1])-1)+"\n"+str(seq)+"\n")
													n.close()
													m = open(folder+"PutativeToxins_annotation_prot.fasta","a+")
													m.write(">"+rec.query.split()[0]+"\n"+prot3+"\n")
													m.close()
													#print(">"+rec.query.split()[0]+"||_"+str(start+b)+"_"+str((len(prot3)*3)+start+b-1)+"_sp"+str(int(sig[1])-1))
													count = 1
												else:
													prot3 = prot3[1:]
											if prot3 == "M":
												break
											if prot3.count("M") ==1 and prot3[-2] =="M":
												break
											if count == 1:
												break
										prot2 = prot2[num+1:]
	os.system("rm "+folder+"signalp.*")
#END
