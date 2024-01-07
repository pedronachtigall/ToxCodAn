#!/usr/bin/env python3
'''
ToxCodAn - Toxin Annotator of Transcriptome Data
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

##modules
import os
import sys
import datetime as dt
from optparse import OptionParser
try:
    from Bio import SeqIO
except:
    print('''ToxCodAn was not able to run due to the ERROR below:
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()

def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        if "||" in ID:
            id1 = ID.split("||")
            ID = id1[0]
            final[ID] = [SEQ, id1[1]]
        else:
            final[ID] = SEQ
    return final

def _ParseFastaRR_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _ParseFastaInv_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        if SEQ in final.keys():
            final[SEQ].append(ID)
        if SEQ not in final.keys():
            final[SEQ] = [ID]
    return final

def _ParseFastaPARTIAL_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

def _RemoveRedundancy_(fasta, outF):
    FASTAinv = _ParseFastaInv_(fasta)
    FASTA = _ParseFastaRR_(fasta)
    report = open(outF+"RemoveRedundancy.log","w")
    output = fasta.replace("_Toxins_cds.fasta","_Toxins_cds_RedundancyFiltered.fasta")
    OUT = open(output,"w")
    countSEQ = 0
    for k in FASTAinv.keys():
        alfa = FASTAinv[k]
        ID = alfa[0]
        SEQ = FASTA[ID]
        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
        OUT.write(">"+ID+"\n"+seq+"\n")
        report.write("//Cluster_"+str(countSEQ)+"\n")
        countID = 0
        for i in alfa:
            report.write(str(countID)+"\t"+i+"\n")
            countID += 1
        countSEQ += 1
    OUT.close()
    report.close()

def _PartialFilter_(sample, output, partial):
    CDSs = [sample+"_Toxins_cds.fasta", sample+"_PutativeToxins_cds.fasta" , sample+"_PutativeToxins_cds_SPfiltered.fasta"]
    if partial == "False":
        PARTIAL = []
        for i in CDSs:
            F = _ParseFastaPARTIAL_(output+i)
            OUT = open(output+i,"w")
            for k in F.keys():
                SEQ = F[k]
                if SEQ[-3:] in ["TAA", "TAG", "TGA"]:
                    seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                    OUT.write(">"+k+"\n"+seq+"\n")
                if SEQ[-3:] not in ["TAA", "TAG", "TGA"]:
                    seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                    OUT.write(">"+k+"_PARTIAL\n"+seq+"\n")
                    PARTIAL.append(k)
            OUT.close()
        for i in os.listdir(output):
            if "_pep" in i:
                F = _ParseFastaPARTIAL_(output+i)
                OUT = open(output+i,"w")
                for k in F.keys():
                    SEQ = F[k]
                    if k in PARTIAL:
                        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                        OUT.write(">"+k+"_PARTIAL\n"+seq+"\n")
                    if k not in PARTIAL:
                        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                        OUT.write(">"+k+"\n"+seq+"\n")
                OUT.close()
    if partial != "False":
        removed = []
        for i in CDSs:
            F = _ParseFastaPARTIAL_(output+i)
            OUT = open(output+i,"w")
            for k in F.keys():
                SEQ = F[k]
                if SEQ[-3:] in ["TAA", "TAG", "TGA"]:
                    seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                    OUT.write(">"+k+"\n"+seq+"\n")
                if SEQ[-3:] not in ["TAA", "TAG", "TGA"]:
                    removed.append(k)
            OUT.close()
        for i in os.listdir(output):
            if ("_contigs" in i or "_pep" in i) and "_NonToxins_" not in i:
                F = _ParseFastaPARTIAL_(output+i)
                OUT = open(output+i,"w")
                for k in F.keys():
                    if k not in removed:
                        SEQ = F[k]
                        seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
                        OUT.write(">"+k+"\n"+seq+"\n")
                OUT.close()

def _ToxCodAn_(sample, transcripts, output, model, signalp, N, covpre, covtoxin, partial):
    blastDB = model+"blastDB/"
    ##pre-filter step
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> performing pre-filter step...")
    #generating blastDB for pre-filter step
    os.system("mkdir "+output+"blastDB/")
    os.system("makeblastdb -dbtype nucl -in "+transcripts+" -out "+output+"blastDB/WC")
    #blast CDSs of toxins against the transcripts
    os.system("tblastn -query "+model+"seqs/FinalToxinDB.fa -out "+output+"PreFilterBlast.out -db "+output+"blastDB/WC -num_threads "+N+" -qcov_hsp_perc "+covpre+" -max_target_seqs 12 -outfmt \'6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore\'")
    #pre-filter
    os.system("pre_filter.py "+output+"PreFilterBlast.out "+transcripts+" "+output+"PreFilter_WC.fa")

    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> performing toxin prediction...")

    ##prediction step
    #run prediction using several models
    M = ["TOXIN","TOXIN_ACHE","TOXIN_HYAL","TOXIN_MYO","TOXIN_PDE","TOXIN_VEGF",
        "TOXIN_3FTx","TOXIN_CTL","TOXIN_KUN","TOXIN_NUC","TOXIN_PLA2"]
    for tm in M:
        os.system("codan.py -t "+output+"PreFilter_WC.fa -o "+output+tm+"_out -m "+model+tm+" -c "+N)

    ##toxin filtering step
    #blast against specific toxinDB to filter toxin preditions
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> filtering toxin predictions...")

    #add blast search parameters
    for tm in M:
        os.system("blastx -query "+output+tm+"_out/ORF_sequences.fasta -out "+output+"BLAST_"+tm+".out -db "+blastDB+"TOXINS -num_threads "+N+" -qcov_hsp_perc "+covtoxin+" -max_target_seqs 10 -strand plus -outfmt \'6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore\'")

    #parse blast and filter the toxins
    os.system("toxin_filter.py "+output+" "+transcripts+" "+sample)

    ##generate final results
    WholeContigs = _ParseFasta_(transcripts)
    PutativeToxins = _ParseFasta_(output+"PreFilter_WC.fa")
    CDStoxins = _ParseFasta_(output+sample+"_Toxins_cds.fasta")

    Toxins = open(output+sample+"_Toxins_contigs.fasta","w")
    UnknownToxins = open(output+sample+"_PutativeToxins_contigs.fasta","w")
    NonToxins = open(output+sample+"_NonToxins_contigs.fasta","w")
    for k in sorted(WholeContigs.keys()):
        if k in CDStoxins.keys():
            SEQ = WholeContigs[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            hit = CDStoxins[k][1]
            Toxins.write(">"+k+" "+hit+"\n"+seq+"\n")
        if k in PutativeToxins.keys() and k not in CDStoxins.keys():
            SEQ = WholeContigs[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            UnknownToxins.write(">"+k+"\n"+seq+"\n")
        if k not in PutativeToxins.keys() and k not in CDStoxins.keys():
            SEQ = WholeContigs[k]
            seq = "\n".join([SEQ[n:n+100] for n in range(0, len(SEQ), 100)])
            NonToxins.write(">"+k+"\n"+seq+"\n")
    Toxins.close()
    UnknownToxins.close()
    NonToxins.close()

    ##PutativeToxins filtering step
    #blast search to detective CDSs in the PutativeToxins
    os.system("blastx -query "+output+sample+"_PutativeToxins_contigs.fasta -out "+output+"PutativeToxinsBlast.xml -db "+blastDB+"TOXINS -num_threads "+N+" -max_target_seqs 10 -outfmt 5 -evalue 0.0001")
    os.system("BlastNamer.py "+output+sample+"_PutativeToxins_contigs.fasta "+output+"PutativeToxinsBlast.xml "+output)
    #blast search on the PutativeToxins
    os.system("blastp -query "+output+"PutativeToxins_annotation_prot.fasta -db "+blastDB+"TOXINS -qcov_hsp_perc 80 -out "+output+"PutativeToxins_blast.out -max_target_seqs 1 -outfmt \'6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore\'")
    #cleaning the PutativeToxins
    os.system("PutativeToxins_filter.py "+output+"PutativeToxins_blast.out "+output+"PutativeToxins_annotation.fasta "+output+" "+sample)

    ##Partial rename and filtering
    _PartialFilter_(sample, output, partial)

    ##Redundancy filtering
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> removing redundancy from the toxin predictions with 100% identity...")
    _RemoveRedundancy_(output+sample+"_Toxins_cds.fasta", output)

    ##signalP filtering
    if signalp == "True":
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> filtering toxin predictions with no signalp...")
        os.system("signalp -u 0.34 -U 0.34 -n "+output+"signalp_annotation.gff "+output+sample+"_Toxins_pep.fasta > "+output+"signalp_summary.txt")
        os.system("signalp_filter.py "+output+"signalp_annotation.gff "+output+sample)
        os.remove(output+"signalp_summary.txt")

    ##Clean the output folder
    try:
        os.system("rm "+output+"PreFilterBlast.out")
        os.system("rm "+output+"*.out")
        os.remove(output+"PreFilter_WC.fa")
        os.remove(output+"PutativeToxinsBlast.xml")
        os.remove(output+"PutativeToxins_annotation.fasta")
        os.remove(output+"PutativeToxins_annotation_prot.fasta")
    except:
        do=0

def _NonToxinsAnnotation_(sample, output, model, N):
    nontoxins = output+sample+"_NonToxins_contigs.fasta"
    #run CodAn to predict the CDS
    os.system("codan.py -t "+output+sample+"_NonToxins_contigs.fasta -o "+output+"NonToxins_CodAn -m "+model+"VERT_full/ -c "+N)
    #run annotation on predicted CDSs using the pepDB
    os.system("NonToxinAnnotation.py -t "+output+"NonToxins_CodAn/ORF_sequences.fasta -d "+model+"pepDB/pepDB"+" -c "+N+" -o "+output+"NonToxins_Annotation")

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-s", "--sample", dest="sample", help="Optional - sample ID to be used in the output files [default=toxcodan]", metavar="string", default="toxcodan")
    parser.add_option("-t", "--transcripts", dest="transcripts", help="Mandatory - transcripts in FASTA format, /path/to/transcripts.fasta", metavar="fasta", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - output folder, /path/to/output_folder; if not defined, the output folder will be set in the current directory [ToxCodAn_output]", metavar="folder", default=None)
    parser.add_option("-m", "--model", dest="model", help="Mandatory - path to model folder, /path/to/models", metavar="path", default=None)
    parser.add_option("-p", "--signalp", dest="signalp", help="Optional - turn on/off the signalP filtering step, use True to turn on or False to turn off [default=True]", metavar="boolean value", default="True")
    parser.add_option("-P", "--partial", dest="partial", help="Optional - turn on/off the partial filtering step, use True to turn on or False to turn off [default=False]", metavar="boolean value", default="False")
    parser.add_option("-n", "--nontoxinannotation", dest="nontoxins", help="Optional - path to folder containing the protein DB and CodAn model to be used in the NonToxin Annotation pipeline [default=None]", metavar="path", default="False")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")
    parser.add_option("-f", "--covprefilter", dest="covprefilter", help="Optional - threshold value used as the minimum coverage in the pre-filter step [default=90]", metavar="int", default="90")
    parser.add_option("-F", "--covtoxinfilter", dest="covtoxinfilter", help="Optional - threshold value used as the minimum coverage in the toxin filter step [default=80]", metavar="int", default="80")

    (options, args) = parser.parse_args()

    if options.transcripts == None or options.model == None:
        print(
        """

▄▄▄█████▓ ▒█████  ▒██   ██▒ ▄████▄   ▒█████  ▓█████▄  ▄▄▄       ███▄    █
▓  ██▒ ▓▒▒██▒  ██▒▒▒ █ █ ▒░▒██▀ ▀█  ▒██▒  ██▒▒██▀ ██▌▒████▄     ██ ▀█   █
▒ ▓██░ ▒░▒██░  ██▒░░  █   ░▒▓█    ▄ ▒██░  ██▒░██   █▌▒██  ▀█▄  ▓██  ▀█ ██▒
░ ▓██▓ ░ ▒██   ██░ ░ █ █ ▒ ▒▓▓▄ ▄██▒▒██   ██░░▓█▄   ▌░██▄▄▄▄██ ▓██▒  ▐▌██▒
  ▒██▒ ░ ░ ████▓▒░▒██▒ ▒██▒▒ ▓███▀ ░░ ████▓▒░░▒████▓  ▓█   ▓██▒▒██░   ▓██░
  ▒ ░░   ░ ▒░▒░▒░ ▒▒ ░ ░▓ ░░ ░▒ ▒  ░░ ▒░▒░▒░  ▒▒▓  ▒  ▒▒   ▓▒█░░ ▒░   ▒ ▒
    ░      ░ ▒ ▒░ ░░   ░▒ ░  ░  ▒     ░ ▒ ▒░  ░ ▒  ▒   ▒   ▒▒ ░░ ░░   ░ ▒░
  ░      ░ ░ ░ ▒   ░    ░  ░        ░ ░ ░ ▒   ░ ░  ░   ░   ▒      ░   ░ ░
             ░ ░   ░    ░  ░ ░          ░ ░     ░          ░  ░         ░
                           ░                  ░



>>>> ToxCodAn v1.0 January 2021 <<<<
      ****Use -h for help!****

USAGE:
toxcodan.py -t transcripts.fa -m path/to/models
        """)
        quit()

    if options.output != None:
        if not options.output.endswith("/"):
            options.output += "/"
        #if options.output.endswith("/"):
            #options.output += "ToxCodAn_output/"
    if options.output == None:
        CWD = os.getcwd()
        options.output = CWD+"/ToxCodAn_output/"


    if os.path.isdir(options.output) == False:
        os.mkdir(options.output)

    if options.model != None and not options.model.endswith("/"):
        options.model += "/"

    if options.transcripts != None and options.model != None:
        print("""

▄▄▄█████▓ ▒█████  ▒██   ██▒ ▄████▄   ▒█████  ▓█████▄  ▄▄▄       ███▄    █
▓  ██▒ ▓▒▒██▒  ██▒▒▒ █ █ ▒░▒██▀ ▀█  ▒██▒  ██▒▒██▀ ██▌▒████▄     ██ ▀█   █
▒ ▓██░ ▒░▒██░  ██▒░░  █   ░▒▓█    ▄ ▒██░  ██▒░██   █▌▒██  ▀█▄  ▓██  ▀█ ██▒
░ ▓██▓ ░ ▒██   ██░ ░ █ █ ▒ ▒▓▓▄ ▄██▒▒██   ██░░▓█▄   ▌░██▄▄▄▄██ ▓██▒  ▐▌██▒
  ▒██▒ ░ ░ ████▓▒░▒██▒ ▒██▒▒ ▓███▀ ░░ ████▓▒░░▒████▓  ▓█   ▓██▒▒██░   ▓██░
  ▒ ░░   ░ ▒░▒░▒░ ▒▒ ░ ░▓ ░░ ░▒ ▒  ░░ ▒░▒░▒░  ▒▒▓  ▒  ▒▒   ▓▒█░░ ▒░   ▒ ▒
    ░      ░ ▒ ▒░ ░░   ░▒ ░  ░  ▒     ░ ▒ ▒░  ░ ▒  ▒   ▒   ▒▒ ░░ ░░   ░ ▒░
  ░      ░ ░ ░ ▒   ░    ░  ░        ░ ░ ░ ▒   ░ ░  ░   ░   ▒      ░   ░ ░
             ░ ░   ░    ░  ░ ░          ░ ░     ░          ░  ░         ░
                           ░                  ░


        """)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting ToxCodAn (v1.0 January 2021)...")
        CWD = os.getcwd()
        print("\ttranscripts file ->", options.transcripts)
        print("\tOutput files and folders ->", options.output)
        print("\tModels folder ->", options.model)
        print("\tsignalP filter -> ", options.signalp)
        print("\tpartial filter -> ", options.partial)
        print("\tNumber of threads ->", options.cpu)

        #check if fasta is lowercase
        case = []
        a = open(options.transcripts, "r")
        for line in a:
            if not line.startswith(">"):
                if line.strip().isupper():
                    case.append("UPPER")
                elif line.strip().islower():
                    case.append("LOWER")
                else:
                    case.append("NONE")
        a.close()
        if "LOWER" in case or "NONE" in case:
            a = open(options.transcripts, "r")
            upper = open(options.output+"transcripts.fasta", "w")
            for line in a:
                if line.startswith(">"):
                    upper.write(line)
                if not line.startswith(">"):
                    upper.write(line.upper())
            a.close()
            upper.close()
            options.transcripts = options.output+"transcripts.fasta"

        _ToxCodAn_(options.sample,
                    options.transcripts,
                    options.output,
                    options.model,
                    str(options.signalp),
                    str(options.cpu),
                    str(options.covprefilter),
                    str(options.covtoxinfilter),
                    str(options.partial))

        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> ToxCodAn finished!")
        print("\tCheck the final results at", options.output)
        print("\t\t"+options.sample+"_Toxins_cds.fasta")
        print("\t\t"+options.sample+"_Toxins_pep.fasta")
        print("\t\t"+options.sample+"_Toxins_annotation.gtf")
        print("\t\t"+options.sample+"_Toxins_contigs.fasta")
        print("\t\t"+options.sample+"_PutativeToxins_cds.fasta")
        print("\t\t"+options.sample+"_PutativeToxins_contigs.fasta")
        print("\t\t"+options.sample+"_NonToxins_contigs.fasta")

        if str(options.signalp) == "True":
            print("\t\t"+options.sample+"_Toxins_cds_SPfiltered.fasta")
            print("\t\t"+options.sample+"_Toxins_pep_SPfiltered.fasta")
            print("\t\t"+options.sample+"_Toxins_contigs_SPfiltered.fasta")
            print("\t\t"+options.sample+"_Toxins_cds_SPfiltered_RedundancyFiltered.fasta")

        print('''
        Description of the output files:
            cds -> coding sequence of the predicted toxins
            pep -> protein sequence of the predicted toxins
            contigs -> whole contigs containing the predicted CDSs
            Toxins -> sequences with very high probability of being toxins
            PutativeToxins -> sequences with medium/high probability of being toxins
            NonToxins -> sequences with very low probability of being toxins
            RedundancyFiltered -> CDSs with 100% identity filtered
            SPfiltered -> signalP filtered sequences (optional step)
        ''')

        if options.nontoxins != "False":
            print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Starting NonToxins Annotation!")

            _NonToxinsAnnotation_(options.sample,
                    options.output,
                    options.nontoxins,
                    str(options.cpu))

            print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> NonToxins Annotation finished!")

if __name__ == '__main__':
	__main__()

#END
