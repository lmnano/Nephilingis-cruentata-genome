#!/bin/bash

GENOMES=../00_data/01_repeatDB/genomeList.txt
OUT=../00_data/01_repeatDB/


REPEATDB=$(ls ../00_data/01_repeatDB/*/*-families.fa)

for r in $REPEATDB;
do

	GENNAME=$(basename "$r" -families.fa)
	GENOME=$(cat $GENOMES | grep $GENNAME)

        singularity exec -e -B $(realpath ..) ../dfam-tetools-latest.sif RepeatMasker \
		-lib "$r" \
		-xsmall \
		-nolow \
		-pa 16 \
		$GENOME


done



