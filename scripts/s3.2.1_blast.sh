#!/usr/bin/sh
cat ./raw_data/Araport11.pep.fasta ./raw_data/GLYMA2.pep.fasta > ./processed_data/ATHGMA.pep.fasta


### Build BLAST database  ===============================================================
cd processed_data
makeblastdb -in ATHGMA.pep.fasta \
            -out ATHGMA.pep.blastdb \
            -dbtype prot \
            -logfile makeblastdb.log

### Perform BLAST analysis  ===============================================================
blastp -evalue 0.00001 \
  -outfmt 6 -db ATHGMAX.pep.blastdb \
  -query ATHGMA.pep.fasta > ATHGMA.pep.blastout 




