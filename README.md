# Readme - work in progress

## Chromosome determination

Finding and determening sex chromosomes is in the folder chr_det. (NEED TO UPLOAD)

### 00 QC

First QC was performed on the raw data using FastQC with the 00_fastQC.sh script. Input data is the folder with raw reads (not available on this repository).

### 01 Read trimming

Read trimming was done with TrimGalore!, using the 01_trimgalore.sh script. The inputs are raw reads.

### 02 Indexing genome

In this step index of the genome was created with bwa index, for later alignment with bwa mem. The script used was 02_bwa_index.sh, input is the genome (not available in this repository).

### 03 Read mapping

Here trimmed reads were mapped to the genome using bwa mem with the script 03_bwa_mapping.sh. Inputs for the script are the genome index and the trimmed reads. Script pipes bwa mem output to samtools sort to create sorted bam files and later uses samtools index for indexing the sorted bam files.

### 04 Converting to bed

This step converts result from the previous step to bed format using bedtools bamtobed using the 04_bwa_bam_to_bed.sh script. Inputs are bed files from the previous step.

## Synteny analysis

Synteny analsis is in the folder genome_synteny (NEED TO MOVE)

### 01 Repeat DB creation

First databases were created using BuildDatabase command from dfam-tetools-latest.sif container. These were used for input for RepeatModeler in the next step. Input is a list of all the genomes with paths to their location. The genomes are not available in this repository and should be downloaded separately from databases and put in an appropriate folder. Databases were built for each genome separately, making one database with all the genomes used didn't work. This was done with 01_repeatDB.sh script.

### 02 RepeatModeler

Next step was running RepeatModeler, from tetools_latest.sif container. RepeatModeler was used to create a repeat database for each genome, -LTRStruct option was used. Input is a list of genomes.
Two species did not run with -LTRStruct, due to the software getting stuck, those were run separately using 02_repeatModelerEdavidi.sh and 02_repeatModelerPpseudoannulata.sh.

### 03 RepeatMasker

RepeatMasker was used next, again from dfam-tetools-latest.sif container. Inputs are list of genomes and a list of RepeatModeler outputs. Options used are -xsmall and -nolow, Cactus should be able to figure out low complexity repeats. 
Script used was 03_repeatmasker.sh.

### 04 RepeatMasker Result Summary

This was in intermediate step to check the results of RepeatMasker run using 04_repeatmodeler_results.sh. It uses the RepeatMasker .tbl results to make the summary and outputs a .tsv table.

### 05 Cactus Whole Genome Alignment

Next whole genome alignments were made with Progressive Cactus. This was done with 05.1_cactusFull.sh script. Inputs are the Cactus seqfile containing the phylogenetic tree for the species in question with the locations for the masked genomes, and directories for job store, work dir and final output. It was run with the --logInfo option.

### 06 Synteny analysis

Synteny analysis was done in the same way for both analyses from the prevoius step. Scripts used were 06_cactus_synteny_full_*.sh. Species in the name represents the target species.
Software used was halSynteny that comes packaged in the Progressive Cactus container. The main input is the .hal file, the result from the previous step. Since halSynteny does only pairwise comparisons for synteny the script takes one species from the .hal file as the target species and loops through all the other species in a for loop. As a result of this, it requires a list of all other species to loop through.

These lists are provided as the 06_synteny_list_full_*.txt files where the species in the name represents the target species. They were made using the 06.1_synteny_sp_list_for_full.R script, which takes the all_sp.xlsx table as the input.

The final result of this step result is a .psl table.

### 07 Synteny results analysis

Synteny results were analysed in the 07_synteny_results_analysis.R R script. The input the script needs is the all_sp.xlsx table with input chromosome information, and the output tableFile. Outputs are tables with best matching chromosome matches for a pairwise comparison between two species.
