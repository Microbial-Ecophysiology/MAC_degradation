#!/bin/bash

assembly=$1
bins=$2
readmap=$3
SID=$4
THRESH=$5
NCPU=$6

#Metabat2 Binning SPAdes

mkdir -p $bins/metabat_spades/
jgi_summarize_bam_contig_depths --outputDepth $readmap/bowtie_spades/depth.txt $readmap/bowtie_spades/*.bam
metabat2 -i $assembly/SPAdes/contigs_fixed.fa -a $readmap/bowtie_spades/depth.txt -o $bins/metabat_spades/${SID} -m $THRESH -v -t $NCPU 


#Metabat2 Binning Megahit

mkdir -p $bins/metabat_megahit/
jgi_summarize_bam_contig_depths --outputDepth $readmap/bowtie_megahit/depth.txt $readmap/bowtie_megahit/*.bam
metabat2 -i $assembly/Megahit/contigs_fixed.fa -a $readmap/bowtie_megahit/depth.txt -o $bins/metabat_megahit/${SID} -m $THRESH -v -t $NCPU 

