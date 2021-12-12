#!/usr/bin/env python3
#Script designed to parse BLAST results and annotate sequences
#Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except:
    print('''ERROR!!!
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None
import os
import datetime as dt
from optparse import OptionParser

#>>>>Translate CDSs
def _translateCDS_(nonannotated):
    translated = {}
    for k in nonannotated.keys():
        if generic_dna:
            cds = Seq(nonannotated[k], generic_dna)
        else:
            cds = Seq(nonannotated[k])
        protein = str(cds.translate())
        translated[k] = protein
    return translated

#>>>>Parse Fasta files
def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

#>>>>Parse Blast/Diamond results
def _ParseBlast_(blast, report):
    final = {}
    count = 0
    a = open(blast,"r")
    for line in a:
        line1 = line.strip().split("\t")
        qID = line1[0]
        qlen = int(line1[1])
        sID = line1[2]
        slen = int(line1[3])
        ident = float(line1[4])
        alen = int(line1[5])
        description = line1[-1]
        if alen >= slen * 0.90 and ident >= 60.0:
            count += 1
            final.setdefault(qID, [])
            final[qID].append((sID, description))
    a.close()
    b = open(report,"w")
    for k in final.keys():
        b.write(k+" -> ")
        for i in final[k]:
            b.write(i[0]+"||"+i[1]+" <> ")
        b.write("\n")
    b.close()
    return final

##finish this
#>>>>Perform BLAST/DIAMOND search
def _BLASTsearch_(transcripts, outF, blastdb, coverage, evalue, cpu, search):
    annotated = {}
    nonannotated = {}
    #if search engine is blast
    if search == "blast":
        if "," in blastdb: #if it has more than one blastDB
            CDSs = _ParseFasta_(transcripts)
            temp = open(outF+"temp.fa","w")
            for k in CDSs.keys():
                temp.write(">"+k+"\n"+CDSs[k]+"\n")
            temp.close()
            dblist = blastdb.split(",")
            for n in range(0, len(dblist)):
                print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Performing BLAST search using DB"+str(n+1)+"!!!")
                temp = outF+"temp.fa"
                OUT = outF+"blast_"+str(n+1)+".out"
                os.system("blastx -qcov_hsp_perc "+coverage+" -num_threads "+cpu+" -evalue "+evalue+" -max_target_seqs 5 -db "+dblist[n]+" -query "+temp+" -out "+OUT+" -outfmt \'6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore stitle\'")
                report = outF+"blast_"+str(n+1)+"_report.txt"
                parsed = _ParseBlast_(OUT, report)

                temp = open(outF+"temp.fa","w")
                for k in CDSs.keys():
                    if k in parsed.keys():
                        annotated[k] = [CDSs[k]]
                        for i in parsed[k]:
                            annotated[k].append(i[0]+" "+i[1])
                    else:
                        temp.write(">"+k+"\n"+CDSs[k]+"\n")
                temp.close()

            for k in CDSs.keys():
                if k not in annotated.keys():
                    nonannotated[k] = CDSs[k]
            os.remove(outF+"temp.fa")

        else: # if has only one blastDB
            print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Performing BLAST search!!!")
            OUT = outF+"blast.out"
            os.system("blastx -qcov_hsp_perc "+coverage+" -num_threads "+cpu+" -evalue "+evalue+" -max_target_seqs 5 -db "+blastdb+" -query "+transcripts+" -out "+OUT+" -outfmt \'6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore stitle\'")
            report = outF+"blast_report.txt"
            parsed = _ParseBlast_(OUT, report)

            CDSs = _ParseFasta_(transcripts)
            for k in CDSs.keys():
                if k in parsed.keys():
                    annotated[k] = [CDSs[k]]
                    for i in parsed[k]:
                        annotated[k].append(i[0]+" "+i[1])
                else:
                    nonannotated[k] = CDSs[k]

    #if search engine is diamond
    if search == "diamond":
        if "," in blastdb: #if it has more than one DB
            CDSs = _ParseFasta_(transcripts)
            temp = open(outF+"temp.fa","w")
            for k in CDSs.keys():
                temp.write(">"+k+"\n"+CDSs[k]+"\n")
            temp.close()
            dblist = blastdb.split(",")
            for n in range(0, len(dblist)):
                print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Performing DIAMOND search using DB"+str(n+1)+"!!!")
                temp = outF+"temp.fa"
                OUT = outF+"blast_"+str(n+1)+".out"
                os.system("diamond blastx --query-cover "+coverage+" --threads "+cpu+" --evalue "+evalue+" --max-target-seqs 5 --strand plus --ultra-sensitive -d "+dblist[n]+" -q "+transcripts+" -o "+OUT+" --outfmt 6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore stitle")
                report = outF+"blast_"+str(n+1)+"_report.txt"
                parsed = _ParseBlast_(OUT, report)

                temp = open(outF+"temp.fa","w")
                for k in CDSs.keys():
                    if k in parsed.keys():
                        annotated[k] = [CDSs[k]]
                        for i in parsed[k]:
                            annotated[k].append(i[0]+" "+i[1])
                    else:
                        temp.write(">"+k+"\n"+CDSs[k]+"\n")
                temp.close()

            for k in CDSs.keys():
                if k not in annotated.keys():
                    nonannotated[k] = CDSs[k]
            os.remove(outF+"temp.fa")

        else: # if has only one DB
            print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Performing DIAMOND search!!!")
            OUT = outF+"blast.out"
            #diamond blastx --query-cover 90 --threads 2 --max-target-seqs 5 --strand plus --ultra-sensitive -d diamondDB/pepDB -q Mfulv_NT_codan_RR.fa -o diamond_out.txt --outfmt 6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore stitle
            os.system("diamond blastx --query-cover "+coverage+" --threads "+cpu+" --evalue "+evalue+" --max-target-seqs 5 --strand plus --ultra-sensitive -d "+blastdb+" -q "+transcripts+" -o "+OUT+" --outfmt 6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore stitle")
            report = outF+"blast_report.txt"
            parsed = _ParseBlast_(OUT, report)

            CDSs = _ParseFasta_(transcripts)
            for k in CDSs.keys():
                if k in parsed.keys():
                    annotated[k] = [CDSs[k]]
                    for i in parsed[k]:
                        annotated[k].append(i[0]+" "+i[1])
                else:
                    nonannotated[k] = CDSs[k]

    return annotated, nonannotated

#>>>>Perform hmm search - busco
def _BUSCOsearch_(nonannotated, outF, busco, cpu):
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Performing HMM search with BUSCO hmm!!!")
    final = {}
    translated = _translateCDS_(nonannotated)
    hmmout = outF+"hmm/"
    if os.path.isdir(hmmout) == False:
        os.mkdir(hmmout)
    temp = open(outF+"HMMsearch_pep.fa", "w")
    for k in translated.keys():
        temp.write(">"+k+"\n"+translated[k]+"\n")
    temp.close()
    peptide = outF+"HMMsearch_pep.fa"
    if busco.endswith("/"):
        lineagepath = busco+"hmms/"
    else:
        lineagepath = busco+"/hmms/"
    for hmm in os.listdir(lineagepath):
        os.system("hmmsearch --cpu "+cpu+" --domtblout "+hmmout+hmm.replace(".hmm",".blout")+" -o "+hmmout+hmm.replace(".hmm", ".out")+" "+lineagepath+hmm+" "+peptide)
        if os.path.isfile(hmmout+hmm.replace(".hmm",".blout")) == True:
            a = open(hmmout+hmm.replace(".hmm",".blout"),"r")
            for line in a:
                if not line.startswith("#"):
                    line1 = line.strip().split()
                    id = line1[0]
                    qlen = int(line1[2])
                    slen = int(line1[5])
                    evalue = float(line1[6])
                    if qlen > 0.5*slen and evalue < 0.001:
                        final[id] = line1[3]
            a.close()
    b = open(outF+"hmm_report.txt","w")
    for k in final.keys():
        b.write(k+" -> "+final[k]+"\n")
    b.close()
    return final

#Perform hmm scan - pfam
def _PFAMsearch_(nonannotated, outF, pfam, cpu):
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Performing HMM scan with PFAM hmm!!!")
    final = {}
    translated = _translateCDS_(nonannotated)
    temp = open(outF+"HMMsearch_pep.fa", "w")
    for k in translated.keys():
        temp.write(">"+k+"\n"+translated[k]+"\n")
    temp.close()
    pfamout = outF+"pfam.blout"
    peptide = outF+"HMMsearch_pep.fa"
    os.system("hmmscan --cpu "+cpu+" -E 1e-23 --domE 0.2 --domtblout "+pfamout+" "+pfam+" "+peptide+" > "+outF+"pfam.log")
    a = open(pfamout, "r")
    for line in a:
        if not line.startswith("#"):
            line1 = line.strip().split()
            tID = line1[0]
            qID = line1[3]
            acc = float(line1[21])
            tdes = line1[22]
            if acc > 0.90:
                final.setdefault(qID,[])
                final[qID].append(tID+" "+tdes)
    a.close()
    b = open(outF+"pfam_report.txt","w")
    for k in final.keys():
        b.write(k+" -> "+" <> ".join(final[k])+"\n")
    b.close()
    return final

# open file with 'a' to append data at the end (I will do it to the final annotated fasta file!!!)

#>>>>Perform annotation
def _Annotation_(transcripts, outF, blastdb, coverage, evalue, busco, pfam, cpu, search):
    annotated, nonannotated = _BLASTsearch_(transcripts, outF, blastdb, coverage, evalue, cpu, search)
    if busco != None:
        busco_annotated = _BUSCOsearch_(nonannotated, outF, busco, cpu)
    else:
        busco_annotated = {}
    if pfam != None:
        pfam_annotated = _PFAMsearch_(nonannotated, outF, pfam, cpu)
    else:
        pfam_annotated = {}

    #writing output
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Generating annotated fasta file!!!")
    OUT = open(outF+"annotated.fa","w")
    OUTNA = open(outF+"uncharacterized.fa","w")
    for k in annotated.keys():
        SEQ = annotated[k][0]
        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        OUT.write(">"+k+"||"+annotated[k][1]+"\n"+seq+"\n")
    for k in nonannotated.keys():
        if k in busco_annotated.keys() and k in pfam_annotated.keys():
            SEQ = nonannotated[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            OUT.write(">"+k+"||"+busco_annotated[k]+"__"+"_".join(pfam_annotated[k])+"\n"+seq+"\n")
        if k in busco_annotated.keys() and k not in pfam_annotated.keys():
            SEQ = nonannotated[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            OUT.write(">"+k+"||"+busco_annotated[k]+"\n"+seq+"\n")
        if k not in busco_annotated.keys() and k in pfam_annotated.keys():
            SEQ = nonannotated[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            OUT.write(">"+k+"||"+"_".join(pfam_annotated[k])+"\n"+seq+"\n")
        if k not in busco_annotated.keys() and k not in pfam_annotated.keys():
            SEQ = nonannotated[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            OUTNA.write(">"+k+"\n"+seq+"\n")
    OUT.close()
    OUTNA.close()

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-t", "--transcripts", dest="transcripts", help="Mandatory - transcripts in FASTA format, /path/to/transcripts.fasta", metavar="fasta", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - output folder, /path/to/output_folder; if not defined, the output folder will be set in the current directory [Annotation_output]", metavar="folder", default=None)
    parser.add_option("-d", "--db", dest="blastdb", help="Mandatory - path to database generated with protein sequences, /path/to/proteinDB; the user can set more than one proteinDB (i.e., from one to \"n\" databases) separeted by a comma (e.g., path/to/proteinDB1,path/to/proteinDB2, ... ,path/to/proteinDBn)", metavar="path", default=None)
    parser.add_option("-s", "--searchengine", dest="search", help="Optional - set the search engine used: blast or diamond [default=blast]", metavar="str", default="blast")
    parser.add_option("-f", "--coveragefilter", dest="coverage", help="Optional - threshold value used as the minimum coverage in the BLAST/DIAMOND search [default=90]", metavar="int", default="90")
    parser.add_option("-e", "--evalue", dest="evalue", help="Optional - e-value threshold used in the BLAST/DIAMOND search [default=0.001]", metavar="str", default="0.001")
    parser.add_option("-b", "--busco", dest="busco", help="Optional - path to BUSCO lineage to perform hmmsearch and annotated BUSCO genes not detected by the BLAST search", metavar="folder", default=None)
    parser.add_option("-p", "--pfam", dest="pfam", help="Optional - path to PFAM hmm models to perform hmmscan and annotated genes not detected by the BLAST search", metavar="folder", default=None)
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")

    (options, args) = parser.parse_args()

    if options.transcripts == None or options.blastdb == None:
        print('''Script designed to perform automated annotation of coding sequences (CDSs) predicted by CodAn software! (Author: Pedro G. Nachtigall [pedronachtigall@gmail.com])

    BASIC USAGE: NonToxins_annotator.py -t transcripts.fa -d blastDB

    Optional options:
    \'-b\', which will perform search in the CDSS using the BUSCO hmm models.
    \'-p\', which will perform hmm scan using the pfam hmm models.
    \'-f\', set the query coverage filter threshold considered in the BLAST/DIAMOND search [default=90].
    \'-e\', set the evalue threshold considered in the BLAST/DIAMOND search [default=0.001].
    \'-s\', set the search tool used in the annotation pipeline: blast or diamond [default=blast].
    \'-c\', set the number of CPUs to be used in the analysis.

    ***Please notice that the fasta file should have the CDS sequence of transcripts.

    ***Please notice that the \'-d\' option expects a protein DB.
        - Depending on the search tool you set to be used, the DB must be designed specifically to each tool.
        - Check the BLAST instructions at https://www.ncbi.nlm.nih.gov/books/NBK279671/ and the DIAMOND instructions at https://github.com/bbuchfink/diamond/wiki

    ***The user can specify more than one blast db separated by a comma (e.g. path/to/db1,path/to/db2).

    ***Make sure that BLAST (v2.9 or higher) is installed and properly working.

    ***Make sure that DIAMOND (v2.0.6 or higher) is installed and properly working.

    ***If the BUSCO option is specified, make sure that hmmsearch is installed and properly working.

    ***If the PFAM option is specified, make sure that hmmscan is installed and properly working.''')
        quit()

    if options.output != None:
        if not options.output.endswith("/"):
            options.output += "/"
    if options.output == None:
        CWD = os.getcwd()
        options.output = CWD+"/Annotation_output/"

    if os.path.isdir(options.output) == False:
        os.mkdir(options.output)

    if options.transcripts != None and options.blastdb != None:

        print('''
    Script designed to perform automated annotation of coding sequences (CDSs) predicted by CodAn software! (Author: Pedro G. Nachtigall [pedronachtigall@gmail.com])''' )
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting Annotation...")
        if options.search == "blast":
            if "," in options.blastdb:
                dbs = options.blastdb.count(",")+1
                print("\t>>> "+str(dbs)+" databases where set to be used:")
                dblist = options.blastdb.split(",")
                for db in dblist:
                    print("\t\t- "+db)
                    if os.path.isfile(db+".phr") == False:
                        print('''ERROR!!!

        The blastDB indicated is not a valid protein DB.
        Please, indicate a valid protein DB to the \"-d\" option.
        Check the instructions at https://www.ncbi.nlm.nih.gov/books/NBK279671/ to generate your own protein DB, or donwload a valid protein DB from https://ftp.ncbi.nlm.nih.gov/blast/db/''')

            else:
                if os.path.isfile(options.blastdb+".phr") == False:
                    print('''ERROR!!!

        The blastDB indicated is not a valid protein DB.
        Please, indicate a valid protein DB to the \"-d\" option.
        Check the instructions at https://www.ncbi.nlm.nih.gov/books/NBK279671/ to generate your own protein DB, or download a valid protein DB from https://ftp.ncbi.nlm.nih.gov/blast/db/''')
                    quit()

        if options.search == "diamond":
            if "," in options.blastdb:
                dbs = options.blastdb.count(",")+1
                print("\t>>> "+str(dbs)+" databases where set to be used:")
                dblist = options.blastdb.split(",")
                for db in dblist:
                    print("\t\t- "+db)
                    if os.path.isfile(db+".dmnd") == False:
                        print('''ERROR!!!

        The diamondDB indicated is not a valid protein DB.
        Please, indicate a valid protein DB to the \"-d\" option.
        Check the instructions at https://github.com/bbuchfink/diamond/wiki to generate your own protein DB to be used with diamond''')

            else:
                if os.path.isfile(options.blastdb+".dmnd") == False:
                    print('''ERROR!!!

        The diamondDB indicated is not a valid protein DB.
        Please, indicate a valid protein DB to the \"-d\" option.
        Please, indicate a valid protein DB to the \"-d\" option.
        Check the instructions at https://github.com/bbuchfink/diamond/wiki to generate your own protein DB to be used with diamond''')
                    quit()

        _Annotation_(options.transcripts,
                    options.output,
                    options.blastdb,
                    str(options.coverage),
                    str(options.evalue),
                    options.busco,
                    options.pfam,
                    str(options.cpu),
                    options.search
        )

        print('''
    <><><><><><><><><><><><><><><><><><><><><><><><><><><>
    <> The final annotated fasta file is \'annotated.fa\' <>
    <><><><><><><><><><><><><><><><><><><><><><><><><><><>''')

        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Finished!!!")

if __name__ == '__main__':
    __main__()

#END
