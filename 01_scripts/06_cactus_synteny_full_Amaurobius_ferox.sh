#!/bin/bash

ALIGNMENT=../00_data/05.1_cactusFull/05.1_cactusFull.hal
TARGET=Amaurobius_ferox

QUERYLIST=06_synteny_list_full_$TARGET.txt
QUERIES=$(cat $QUERYLIST)

OUT=../00_data/06.1_cactus_synteny_full

mkdir $OUT
mkdir $OUT/$TARGET

for q in $QUERIES;
do
	OUTFILE=06.1_synteny_$TARGET\_"$q".psl

	singularity exec -e -B $(realpath ..) ../cactus_v3.1.3.sif halSynteny \
		--queryGenome "$q" \
		--targetGenome $TARGET \
		$ALIGNMENT \
		$OUT/$TARGET/$OUTFILE

done
