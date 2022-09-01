#!/usr/bin/bash
WORKDIR=$(pwd)
OUTDIR=$2

mkdir -p $WORKDIR/raw_data/fastq/$OUTDIR
mkdir -p $WORKDIR/logs/fastq

# export PATH=$WORKDIR/software/sratoolkit.2.8.2-1-centos_linux64/bin:$PATH

echo "---------[Download SRR.fastq]---------"
for SRR in `cat $1`
do
	echo -e "`date +"%b%d|%T"`\t\c"; echo "$SRR"
	fastq-dump --gzip --split-3 $SRR --outdir $WORKDIR/raw_data/fastq/$OUTDIR >$WORKDIR/logs/fastq/$SRR.log 2>$WORKDIR/logs/fastq/$SRR.err
done
