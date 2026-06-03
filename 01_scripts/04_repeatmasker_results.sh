#!/bin/bash

RESULTS=../00_data/genomes
OUT=../00_data/04_repeatmodeler_result_summary

mkdir $OUT
echo -e "Genome\tGenome_bp\tMasked_%\tInterspersed_%\tRetro_%\tDNA_%\tRC_%\tUnclassified_%" > $OUT/repeatmasker_summary.tsv

TABLES=$(ls $RESULTS/*.tbl)

for tbl in $TABLES;
do
    genome=$(basename "$tbl" .tbl)

    genome_bp=$(grep "total length:" "$tbl" | awk '{print $3}')
    masked_pct=$(grep "bases masked:" "$tbl" | awk '{print $6}')
    interspersed_pct=$(grep "Total interspersed repeats:" "$tbl" | awk '{print $5}')
    retro_pct=$(grep "^Retroelements" "$tbl" | awk '{print $5}')
    dna_pct=$(grep "^DNA transposons" "$tbl" | awk '{print $5}')
    rc_pct=$(grep "^Rolling-circles" "$tbl" | awk '{print $5}')
    unclass_pct=$(grep "^Unclassified:" "$tbl" | awk '{print $5}')

    echo -e "${genome}\t${genome_bp}\t${masked_pct}\t${interspersed_pct}\t${retro_pct}\t${dna_pct}\t${rc_pct}\t${unclass_pct}" \
        >> $OUT/repeatmasker_summary.tsv
done
