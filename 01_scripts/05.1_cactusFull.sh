#!/bin/bash

SEQFILE=05_cactus_full_input.txt
WORKDIR=../00_data/05.1_cactusWorkDir
JOBSTORE=../00_data/05.1_cactusFullJobStore
OUT=../00_data/05.1_cactusFull

mkdir $WORKDIR

mkdir $OUT

singularity exec -e -B $(realpath ..) ../cactus_v3.1.3.sif \
	cactus \
	--workDir $WORKDIR \
	--maxCores 32 \
	--logInfo \
	$JOBSTORE \
	$SEQFILE \
	$OUT/05.1_cactusFull.hal
