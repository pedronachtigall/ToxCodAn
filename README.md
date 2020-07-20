![venomancer_logo](/venomancer_logo.png)

Venomancer
==========
<!---[![Latest GitHub release](https://img.shields.io/github/release/pedronachtigall/Venomancer.svg)](https://github.com/pedronachtigall/Venomancer/releases/latest) -->
<!---[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403273.svg)](https://doi.org/10.5281/zenodo.3403273) -->
<!---[![Published in Genome Biology](https://img.shields.io/badge/published%20in-Genome%20Biology-blue.svg)](https://doi.org/10.1101/gr.214270.116) -->

**Venomancer** is a computational tool designed to detect and annotate toxin genes in transcriptome assembly.

Getting Started
=================

# Installation

Download the master folder and follow the steps below:

```
unzip Venomancer-master.zip
export PATH=$PATH:path/to/Venomancer-master/bin/
```
OR git clone the Venomancer respository and add the bin folder into your PATH:
```
git clone https://github.com/pedronachtigall/Venomancer.git
export PATH=$PATH:path/to/Venomancer/bin/
```

# Requirements

- [Python3](https://www.python.org/) and [Biopython](https://biopython.org/wiki/Download)
  - ```apt-get install python3-biopython```
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE) (libmce-perl)
  - ```apt-get install bioperl libmce-perl```
- [CodAn](https://github.com/pedronachtigall/CodAn/)
- [NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (v2.9.0 or above)
- [SignalP-4.1](http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+4.1)

Ensure that all requirements are working properly.

:warning: If the user wants to install Venomancer and all dependencies using [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html), follow the steps below:
- Create the environment:
    - ```conda create -n venomancer_env python=3.7 biopython perl perl-bioperl perl-mce blast```
- Git clone the Venomancer repository and add to your PATH:
    - ```git clone https://github.com/pedronachtigall/Venomancer.git```
    - ```export PATH=$PATH:path/to/Venomancer/bin/```
- Git clone the CodAn repository and add to your PATH:
    - ```git clone https://github.com/pedronachtigall/CodAn.git```
    - ```export PATH=$PATH:path/to/CodAn/bin/```
- Download the SignalP-4.1, decompress and add it to your PATH:
    - ```tar -xzf signalp-4.1g.Linux.tar.gz```
    - ```export PATH=$PATH:path/to/signalp-4.1/```
- Then, run Venomancer as described in the "Usage" section.
- To activate the environment to run Venomancer just use the command: ```conda activate venomancer_env```
- To deactivate the environment just use the command: ```conda deactivate```

# Models

The model folder contains specific gHMM models and the toxinDB used in the Venomancer pipeline.

Download the [models.zip](https://github.com/pedronachtigall/Venomancer/blob/master/models.zip) file, uncompress (```unzip models.zip```) and specify it to the ```-m``` option of Venomancer command line (```-m path/to/models/```).

# Usage

```
Usage: venomancer.py [options]

Options:
  -h, --help            show this help message and exit
  -s string, --sample=string
                        Optional - sample ID to be used in the output files
                        [default=venomancer]
  -t fasta, --transcripts=fasta
                        Mandatory - transcripts in FASTA format,
                        /path/to/transcripts.fasta
  -o folder, --output=folder
                        Optional - output folder, /path/to/output_folder; if
                        not defined, the output folder will be set in the
                        current directory [Venomancer_output]
  -m path, --model=path
                        Mandatory - path to model folder, /path/to/models
  -p boolean value, --signalp=boolean value
                        Optional - turn on/off the signalP filtering step, use
                        True to turn on or False to turn off [default=True]
  -n path, --nontoxinannotation=path
                        Optional - path to folder containing the protein DB
                        and CodAn model to be used in the NonToxin Annotation
                        pipeline [default=None]
  -c int, --cpu=int     Optional - number of threads to be used in each step
                        [default=1]
  -f int, --covprefilter=int
                        Optional - threshold value used as the minimum
                        coverage in the pre-filter step [default=90]
  -F int, --covtoxinfilter=int
                        Optional - threshold value used as the minimum
                        coverage in the toxin filter step [default=80]
```

Basic usage:
```
venomancer.py -t transcripts.fa -m path/to/models
```

Check our [tutorial](https://github.com/pedronachtigall/Venomancer/tree/master/tutorial) to learn how to use Venomancer.

# Inputs
Venomancer has the following inputs as mandatory:
 - Transcripts in [fasta format](https://en.wikipedia.org/wiki/FASTA_format) through the ```-t``` option.
 - The uncompressed [models](https://github.com/pedronachtigall/Venomancer/blob/master/models.zip) folder through the ```-m``` option

# Outputs

Venomancer outputs the following files:
```
SampleID_Toxins_cds.fasta
SampleID_Toxins_pep.fasta
SampleID_Toxins_annotation.gtf
SampleID_Toxins_contigs.fasta
SampleID_PutativeToxins_cds.fasta
SampleID_PutativeToxins_contigs.fasta
SampleID_NonToxins_contigs.fasta

SampleID_Toxins_cds_SPfiltered.fasta (optional step)
SampleID_Toxins_pep_SPfiltered.fasta (optional step)
SampleID_Toxins_contigs_SPfiltered.fasta (optional step)
SampleID_Toxins_cds_SPfiltered_RedundancyFiltered.fasta (optional step)

signalp_annotation.gff (optional step)
RemoveRedundancy.log
```

Description of the output files:
```
cds -> coding sequence of the predicted toxins
pep -> protein sequence of the predicted toxins
contigs -> whole contigs containing the predicted CDSs
Toxins -> sequences with very high probability of being toxins
PutativeToxins -> sequences with medium/high probability of being toxins
NonToxins -> sequences with very low probability of being toxins
RedundancyFiltered -> CDSs with 100% identity filtered
SPfiltered -> signalP filtered sequences (optional step)
```

Annotation of Non Toxin transcripts
===================================

The user can take advantages of a simple script designed to annotate Non Toxin transcripts named **NonToxinAnnotation.py**. Follow the steps below:
- First, perform the CDS prediction with the "VERT_full" model using [CodAn](https://github.com/pedronachtigall/CodAn) designed by [Nachtigall et al. (2020)](https://doi.org/10.1093/bib/bbaa045)
    - ```codan.py -t path/to/NonToxins_contigs.fasta -m path/to/VERT_full/ -o path/to/output/NonToxins_codan/ -c N```
    - We have a copy of the "VERT_full" in the "non_toxin_models" folder: ```cd path/to/non_toxin_models/``` and ```gzip -d VERT_full```
- Then, use the ```NonToxinAnnotation.py``` on the predicted CDSs.
- This script performs `blast` search (mandatory) and hmm search using `BUSCO` and `Pfam` models (optional).
- The use of a protein DB pre-compiled or designed with `makeblastdb` can be set with the `-b` option.
    - The user can use a DB such as Swissprot and/or the designed protein DB available at the "non_toxin_models" folder (just uncompress the DB ```tar xjf pepDB.tar.bz2```).
    - The user can set one or more DBs by using a comma "," among DBs, which can be any number (from 1 to n).
- Optionally, the user can set any of the [BUSCO models](https://busco.ezlab.org/busco_v4_data.html) to perform hmm search by using the option `-b`.
- Optionally, the user can set the [Pfam models](https://pfam.xfam.org/) to perform hmm search by using the option `-p`. (link for download the pfam.hmm: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
- This script takes advantage of MultiThreading by using the option `-c`.
- Usage: ```NonToxins_Annotator.py -t predicted_CDS.fasta -b path/to/db1,...,path/to/dbn -b path/to/busco/odb -p path/to/pfam.hmm -c N```

:warning: Alternatively, if the user wants to directly perform the NonToxins annotation within the Venomancer pipeline just follow the steps below:
- Enter in the ["non_toxin_models"](https://github.com/pedronachtigall/Venomancer/tree/master/non_toxin_models)
    - ```cd path/to/venomancer/non_toxin_models/```
- Uncompress the proteinDB (```tar xjf pepDB.tar.bz2```) and the CodAn model for Vertebrates (```gzip -d VERT_full.zip```)
- Then, use the option ```-n``` in the venomancer command line to automatically perform the NonToxin annotation:
    - ```venomancer.py -s sampleID -t assembly.fasta -o out_venomancer -m /path/to/models -c 4 -n path/to/non_toxin_models/```

Reference
=========

If you use or discuss **Venomancer**, please cite:

Nachtigall et al., under review

License
=======

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)

Contact
=======
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

Frequently Asked Questions (FAQ)
================================

**[Q1]** What Operation System (OS) do I need to use Venomancer?
  - We tested Venomancer in Linux Ubuntu 16 and 18, and macOS Mojave and Catalina. However, we believe that Venomancer should work on any UNIX OS able to have all dependencies of Venomancer.
 
**[Q2]** How long will take to Venomancer finish the analysis?
  - We tested Venomancer using a personal computer (6-Core i7 with 16Gb memory) and 6 threads (```-c 6```), it took only 55 minutes to finish the analysis by using a de novo dataset with 146,077 sequences. If the user has more threads available for use, the running time will decrease.

**[Q3]** Is Venomancer only available for snake species? :snake:
  - Unfortunately, we only acquired sufficient trainning data for snake toxins. But we are working to get more training data to other venomous taxa and make them available soon. Stay tune!

**[Q4]** When was the Databases in the Venomancer last updated?
  - Our models and databases used in the annoations were last updated in September 2019.
