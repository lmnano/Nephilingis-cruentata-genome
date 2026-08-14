#!/bin/bash

GENOMES=../00_data/01_repeatDB/genomeList.txt
OUT=../00_data/01_repeatDB/

FILES=$(cat $GENOMES)

for f in $FILES;
do
	name=$(basename "$f" .fa)
	mkdir $OUT/$name
	singularity exec -e -B $(realpath ..) ../dfam-tetools-latest.sif \
		 BuildDatabase -name $name "$f"

	mv $name* $OUT/$name


done
