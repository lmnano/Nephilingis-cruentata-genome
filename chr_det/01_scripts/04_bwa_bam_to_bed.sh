#!bin/bash

DATA=../00_data/03_bwa_mapping
OUT=../00_data/04_bwa_bam_to_bed

mkdir $OUT

BAM=$(ls $DATA/*.bam | cat)


for b in $BAM;
do
	NAME=$(basename "$b" .bam)
	bedtools bamtobed -i "$b" > $OUT/$NAME.bed

done
