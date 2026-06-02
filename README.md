# Readme

## 01 Repeat DB creation

First databases were created using BuildDatabase command from dfam-tetools-latest.sif container. These were used for input for RepeatModeler in the next step. Input is a list of all the genomes with paths to their location. Databases were built for each genome separately, making one database with all the genomes used didn't work. This was done with 01_repeatDB.sh script.

## 02 RepeatModeler

Next step was running RepeatModeler, from tetools_latest.sif container. RepeatModeler was used to create a repeat database for each genome, -LTRStruct option was used. In order to make it run faster this step was done in 4 batches with 4 separate scripts, each covering a few genomes. Inputs are separate lists of genomes.
The scripts used were 02_repeatModeler01.sh and the rest with numered up to 04.
Two specise did not run with -LTRStruct, due to the software getting stuck, those were run separately using 02_repeatModelerEdavidi.sh and 02_repeatModelerPpseudoannulata.sh.

## 03 RepeatMasker

RepeatMasker was used next, again from dfam-tetools-latest.sif container. Inputs are list of genomes and a list of RepeatModeler outputs. Options used are -xsmall and -nolow, Cactus should be able to figure out low complexity repeats. 
Script used was 03_repeatmasker.sh.

## 04 RepeatMaskaer Result Summary

This was in intermediate step to check the results of RepeatMasker run using 04_repeatmodeler_results.sh. It uses the RepeatMasker .tbl results to make the summary and outputs a .tsv table.

## 05 Cactus Whole Genome Alignment

Next whole genome alignments were made with Progressive Cactus. This was done twice, once with all the species (05.1_cactusFull.sh), full analysis and once with a selection of species that had their sex chromosomes identified not by syteny (05.2_cactusMin.sh). Inputs are the Cactus seqfile containing the phylogenetic tree for the species in question, and directories for job store, work dir and final output. It was run with the --logInfo option. The tree for the limited number of species analysis was made with 05.1_tree_for_cactus.R R script, using the ape package.

## 06 Synteny analysis

Synteny analysis was done in the same way for both analyses from the prevoius step. Scripts used were 06_cactus_synteny_full_*.sh for the full analysis and 06_cactus_synteny_*sh for the limited one. Species in the name represents the target species.
Software used was halSynteny that comes packaged in the Progressive Cactus container. The main input is the .hal file, the result from the previous step. Since halSynteny does only pairwise comparisons for synteny the script takes one species from the .hal file as the target species and loops through all the other species in a for loop. As a result of this it requires a list of all other species to loop through.

These lists are provided as the 06_synteny_list_full_*.txt and 06_synteny_list_*.txt files where the species in the name represents the target species. They were made using the 06.1_synteny_sp_list_for_full.R script, which takes the all_sp.xlsx table as the input.

The final result of this step result is a .psl table.

## 07 Synteny results analysis

Synteny results were analysed in the 07_synteny_results_analysis.R R script. The input the script needs is the all_sp.xlsx table with input chromosome information, and the output tableFile. Outputs are tables with best matching chromosome matches for a pairwise comparison between two species.
