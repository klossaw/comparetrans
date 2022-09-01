#!/usr/bin/bash
### Identify k-best-hits genes between two species using BLAST results

cd processed_data

###  k = 5, blast results with Arabidopsis query - soybean DB
python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 5 ARATH2GLYMA.BLAST_ARATH2GLYMA.subset.txt

###  k = 5, blast results with soybean query - Arabidopsis DB
python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 5 ARATH2GLYMA.BLAST_GLYMA2ARATH.subset.txt


