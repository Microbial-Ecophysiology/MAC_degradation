# Script to combine annotation results 

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Manuscript_MAC/Metagenomics/annotation/")

# load annotation ####
#load("annotation_parsed.RData")

## kegg ----
kegg <- read.csv("full/concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg <- subset(kegg, kegg_sseqid != "sseqid")

## NR ----
NR <- read.csv("full/concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR <- subset(NR, NR_sseqid != "sseqid")


# Parse annotation ----
# Combine annotation data of kegg, NR, dbCAN and merops
annotation <- full_join(kegg, NR, by = c("qseqid")) %>%
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")


write.table(annotation, "annotation_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)
