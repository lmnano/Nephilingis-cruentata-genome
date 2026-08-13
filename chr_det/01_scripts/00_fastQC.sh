#!/bin/bash

#input sapmles
DATA=../00_data/sex_determination_wgrs/01.RawData
#output folder
OUTPUT=../00_data/00_fastQC

#remove old output and create new one
rm -r $OUTPUT
mkdir $OUTPUT

#get .fq.gz files
SAMPLES=$(ls $DATA/*/*fq.gz | cat)

for s in $SAMPLES;
do
	fastqc -o $OUTPUT -t 4 "$s"

done


