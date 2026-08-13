#!/bin/bash

INDEX=../00_data/genome/F11-0418.Chr.fasta
INPUT=../00_data/01_trimgalore

OUTPUT=../00_data/03_bwa_mapping

rm -r $OUTPUT
mkdir $OUTPUT

SAMPLES=$(ls $INPUT)

for s in $SAMPLES;
do
        INPUTFWD=$(ls $INPUT/"$s"/*val_1.fq.gz)
        INPUTREV=$(ls $INPUT/"$s"/*val_2.fq.gz)

	bwa mem -t 16 $INDEX $INPUTFWD $INPUTREV | samtools sort -@ 16 -o $OUTPUT/"$s".sorted.bam
	samtools index $OUTPUT/"$s".sorted.bam
done
