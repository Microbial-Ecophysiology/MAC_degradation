#!/bin/bash

IID=$1
SID=$2
NCPU=$3


bowtie2 --threads $NCPU -x ${IID}/intermediate_results/Mapping/bowtie_megahit/contigs -1 ${SID}/intermediate_results/clean_reads/final_pure_reads_1.fq -2 ${SID}/intermediate_results/clean_reads/final_pure_reads_2.fq -S ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".sam" 

samtools sort -o ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".bam" ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".sam" 

samtools index ${IID}/intermediate_results/Mapping/bowtie_megahit/${SID}".bam" 



bowtie2 --threads $NCPU -x ${IID}/intermediate_results/Mapping/bowtie_spades/contigs -1 ${SID}/intermediate_results/clean_reads/final_pure_reads_1.fq -2 ${SID}/intermediate_results/clean_reads/final_pure_reads_2.fq -S ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".sam" 

samtools sort -o ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".bam" ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".sam" 

samtools index ${IID}/intermediate_results/Mapping/bowtie_spades/${SID}".bam" 
