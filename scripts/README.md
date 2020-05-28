Scripts
=======

Some scripts to help analyze the toxin genes identified by Venomancer.

- **BreakLines.py**: script designed to break huge sequences in a single line to 100 nucleotides length per line.
    - Usage: ```BreakLines.py input.fasta output_breaklines.fasta```

- **RemoveRedundancy.py**: script designed to remove redundancy from datasets (cluster sequences 100% similar).
    - Usage: ```RemoveRedundancy.py input.fa output_RedundancyRemoved.fa report.txt```

<!---
- **rps2gff.py**: script designed to convert rps-blast result in tabular format into gff annotation.
    - First, download the pfam domain database (the uncompressed file has more than 3.0Gb): ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz
    - Second, run RPS-Blast as follow: ```rpstblastn -query IN.fa -db path/to/Pfam_LE/Pfam -out blast_out.out -evalue 0.001 -max_target_seqs 12 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send slen evalue bitscore salltitles'```
    - Then, run "rps2gff" script: ```python rps2gff.py < blast_out.txt > domain_annotation.gff```
--->

- **NonToxins_Annotator.py**: script designed to annotate the NonToxins detected by the Venomancer pipeline.
    - This script performs `blast` search (mandatory) and hmm search using `BUSCO` and `Pfam` models (optional).
    - The use of a protein DB pre-compiled or designed with `makeblastdb` can be set with the `-b` option.
        - The user can set one or more DBs by using a comma "," among DBs, which can be any number (from 1 to n).
    - Optionally, the user can set any of the [BUSCO models](https://busco.ezlab.org/busco_v4_data.html) to perform hmm search by using the option `-b`.
    - Optionally, the user can set the [Pfam models](https://pfam.xfam.org/) to perform hmm search by using the option `-p`. (link for download the pfam.hmm: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
    - This script takes advantage of MultiThreading by using the option `-c`.
    - Usage: ```NonToxins_Annotator.py -t predicted_CDS.fasta -b path/to/db1,...,path/to/dbn -b path/to/busco/odb -p path/to/pfam.hmm -c N```
