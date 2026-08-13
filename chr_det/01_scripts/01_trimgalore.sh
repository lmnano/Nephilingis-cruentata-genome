#!/bin/bash

# Create variables

INPUT= ls ../00_data/sex_determination_wgrs/01.RawData
OUTPUT=../00_data/01_trimgalore

rmdir $OUTPUT
mkdir $OUTPUT

SAMPLES=$(ls $INPUT | cat)
CURRENTDIR=$(pwd)

for s in $SAMPLES;
do
#	cd $INPUT/"$s"
	INPUTFWD=$(ls $INPUT/"$s"/*1.fq.gz)
	INPUTREV=$(ls $INPUT/"$s"/*2.fq.gz)

	OUTDIR=$OUTPUT/"$s"
	mkdir $OUTDIR

	trim_galore \
		-o $OUTDIR \
		-j 4 \
		--phred33 \
		-q 5 \
		--stringency 1 \
		-e 0.1 \
		--length 75 \
		--trim-n \
		--paired \
		--retain_unpaired \
		--gzip \
		$INPUTFWD \
		$INPUTREV

done;


