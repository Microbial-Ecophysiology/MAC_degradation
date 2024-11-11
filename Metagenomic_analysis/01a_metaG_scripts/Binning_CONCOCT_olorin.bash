#!/bin/bash

assembly=$1
bins=$2
readmap=$3
SID=$4
THRESH=$5
cutup=$6
overlap=$7
NCPU=$8

#CONCOCT Binning SPAdes

mkdir -p $bins/concoct_spades/
cut_up_fasta.py $assembly/SPAdes/contigs_fixed.fa -c $cutup -o $overlap --merge_last -b $bins/concoct_spades/contigs_10K.bed > $bins/concoct_spades/contigs_10K.fa

concoct_coverage_table.py $bins/concoct_spades/contigs_10K.bed $readmap/bowtie_spades/*.bam > $bins/concoct_spades/coverage_table.tsv

concoct --composition_file $bins/concoct_spades/contigs_10K.fa --coverage_file $bins/concoct_spades/coverage_table.tsv -b $bins/concoct_spades/concoct_output/ -l $THRESH -t $NCPU

merge_cutup_clustering.py $bins/concoct_spades/concoct_output/clustering_gt2000.csv > $bins/concoct_spades/concoct_output/clustering_merged.csv

mkdir $bins/concoct_spades/concoct_output/fasta_bins

extract_fasta_bins.py $assembly/SPAdes/contigs_fixed.fa $bins/concoct_spades/concoct_output/clustering_merged.csv --output_path $bins/concoct_spades/concoct_output/fasta_bins


#CONCOCT Binning Megahit

mkdir -p $bins/concoct_megahit/
cut_up_fasta.py $assembly/Megahit/contigs_fixed.fa -c $cutup -o $overlap --merge_last -b $bins/concoct_megahit/contigs_10K.bed > $bins/concoct_megahit/contigs_10K.fa

concoct_coverage_table.py $bins/concoct_megahit/contigs_10K.bed $readmap/bowtie_megahit/*.bam > $bins/concoct_megahit/coverage_table.tsv

concoct --composition_file $bins/concoct_megahit/contigs_10K.fa --coverage_file $bins/concoct_megahit/coverage_table.tsv -b $bins/concoct_megahit/concoct_output/ -l $THRESH -t $NCPU

merge_cutup_clustering.py $bins/concoct_megahit/concoct_output/clustering_gt2000.csv > $bins/concoct_megahit/concoct_output/clustering_merged.csv

mkdir $bins/concoct_megahit/concoct_output/fasta_bins

extract_fasta_bins.py $assembly/Megahit/contigs_fixed.fa $bins/concoct_megahit/concoct_output/clustering_merged.csv --output_path $bins/concoct_megahit/concoct_output/fasta_bins
