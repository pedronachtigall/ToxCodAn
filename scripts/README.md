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

- **ChimeraKiller**: script designed to detect chimeras and remove them from the CDSs identified in the Transcriptome Assembly.
    - Usage: ```ChimeraKiller.py```
