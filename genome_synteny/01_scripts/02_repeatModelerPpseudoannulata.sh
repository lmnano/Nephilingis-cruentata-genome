#!/bin/bash

GENOMES=../00_data/01_repeatDB/Pardosa_pseudoannulata_GCA_032207245.1_NJAU_Ppse_V1_genomic/Pardosa_pseudoannulata_GCA_032207245.1_NJAU_Ppse_V1_genomic
OUT=../00_data/01_repeatDB/


singularity exec -e -B $(realpath ../..) tetools_latest.sif \
	RepeatModeler -database $GENOMES -threads 32


