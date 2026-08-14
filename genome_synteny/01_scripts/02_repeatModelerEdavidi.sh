#!/bin/bash

GENOMES=../00_data/01_repeatDB/Ectatosticta_davidi.genome/Ectatosticta_davidi.genome
OUT=../00_data/01_repeatDB/


singularity exec -e -B $(realpath ../..) tetools_latest.sif \
	RepeatModeler -database $GENOMES -threads 32


