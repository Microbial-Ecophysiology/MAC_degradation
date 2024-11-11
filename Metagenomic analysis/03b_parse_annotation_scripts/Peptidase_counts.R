#Script to combine annotation for peptidase and extracellular peptidase heatmaps

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Manuscript_MAC/Metagenomics/annotation/")

# Peptidase heatmap ----
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

## merops ----
merops <- read.csv("Peptidases/concat_merops_annotation_best.txt", 
                   sep="\t", 
                   header=TRUE, 
                   fill=TRUE, 
                   quote="",
                   col.names=c("qseqid", "merops_sseqid","merops_pident","merops_length","merops_mismatch","merops_gapopen","merops_qstart",
                               "merops_qend","merops_sstart","merops_send","merops_evalue","merops_bitscore","merops_bsr","merops_name","merops_id","merops_subfamily","merops_unit","merops_source"))

merops <- subset(merops, merops_sseqid != "sseqid")


## signalp ----
#first line of output and # in second line need to be removed before loading into R
signalp <- read.csv("Peptidases/concat_signalp_annotation.txt", 
                    sep="\t", 
                    header=TRUE, 
                    fill=TRUE, 
                    quote="")


# combine table with Peptidases in first few columns
Peptidases_table <- full_join(merops, kegg, by = c("qseqid")) %>%
  full_join(., NR, by = c("qseqid")) %>%
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")



# remove all rows from dataframe which don't have any annotation through NR or kegg
Peptidases_short <- Peptidases_table[complete.cases(Peptidases_table[, 3:19]), ]
# 1776 Peptidases remain

Peptidase_Clostridia <- subset(Peptidases_short, 
                               Peptidases_short$genome %in% c("E4_d458_spades_bin_53_orig-contigs",
                                             "G3_d458_megahit_bin_18_orig-contigs",
                                             "G3_d458_spades_bin_18_indstr_refined-contigs",
                                             "G3_d458_spades_bin_33_strict_refined-contigs",
                                             "L3_d458_spades_bin_6_orig_refined-contigs",
                                             "E4_d458_spades_bin_64_permissive_refined-contigs",
                                             "L3_d458_spades_bin_64_permissive_refined-contigs",
                                             "G3_d458_megahit_bin_21_strict_refined-contigs",
                                             "L3_d458_spades_bin_34_indper_refined-contigs",
                                             "Van_BES_megahit_bin_42_indstr_refined-contigs",
                                             "G3_d458_spades_bin_19_indper_refined-contigs",
                                             "Co_bin_102_permissive_refined-contigs",
                                             "E4_d458_spades_bin_8_permissive_refined-contigs",
                                             "G3_d458_megahit_bin_13_indstr_refined-contigs",
                                             "L3_d458_megahit_bin_13_permissive_refined-contigs",
                                             "L3_d458_spades_bin_21_orig_refined-contigs",
                                             "L3_d458_spades_bin_8_indstr_refined-contigs",
                                             "M4_d458_megahit_bin_49_permissive_refined-contigs",
                                             "M4_d458_spades_bin_15_permissive_refined-contigs",
                                             "Co_bin_88_orig_refined-contigs",
                                             "L3_d458_megahit_bin_58_indper-contigs",
                                             "E4_d458_megahit_bin_7_strict_refined-contigs",
                                             "M4_d458_megahit_bin_14_orig_refined-contigs",
                                             "M4_d458_spades_bin_28_indstr_refined-contigs",
                                             "Co_bin_49_orig_refined-contigs",
                                             "L3_d458_spades_bin_31_orig_refined-contigs"))


# Create a matrix of Merops Peptidase families against genomes
matrix_Pep <- table(Peptidase_Clostridia$genome, Peptidase_Clostridia$merops_subfamily)
# melt matrix into tabular view
matrix_Pep_melt <- melt(matrix_Pep)
colnames(matrix_Pep_melt) <- c("genome","merops_subfamily","value")

matrix_Pep_melt$genome <- factor(matrix_Pep_melt$genome, levels = c("E4_d458_spades_bin_53_orig-contigs",
                                 "G3_d458_megahit_bin_18_orig-contigs",
                                 "G3_d458_spades_bin_18_indstr_refined-contigs",
                                 "G3_d458_spades_bin_33_strict_refined-contigs",
                                 "L3_d458_spades_bin_6_orig_refined-contigs",
                                 "E4_d458_spades_bin_64_permissive_refined-contigs",
                                 "G3_d458_megahit_bin_21_strict_refined-contigs",
                                 "L3_d458_spades_bin_34_indper_refined-contigs",
                                 "Van_BES_megahit_bin_42_indstr_refined-contigs",
                                 "G3_d458_spades_bin_19_indper_refined-contigs",
                                 "Co_bin_102_permissive_refined-contigs",
                                 "E4_d458_spades_bin_8_permissive_refined-contigs",
                                 "L3_d458_megahit_bin_13_permissive_refined-contigs",
                                 "G3_d458_megahit_bin_13_indstr_refined-contigs",
                                 "L3_d458_spades_bin_8_indstr_refined-contigs",
                                 "M4_d458_spades_bin_15_permissive_refined-contigs",
                                 "L3_d458_spades_bin_21_orig_refined-contigs",
                                 "M4_d458_megahit_bin_49_permissive_refined-contigs",
                                 "Co_bin_88_orig_refined-contigs",
                                 "L3_d458_megahit_bin_58_indper-contigs",
                                 "E4_d458_megahit_bin_7_strict_refined-contigs",
                                 "M4_d458_megahit_bin_14_orig_refined-contigs",
                                 "M4_d458_spades_bin_28_indstr_refined-contigs",
                                 "Co_bin_49_orig_refined-contigs",
                                 "L3_d458_spades_bin_31_orig_refined-contigs",
                                 "L3_d458_spades_bin_64_permissive_refined-contigs"),
                                 labels = c("Aceti_1","Aceti_2", "Aceti_3","Aceti_4", "Aceti_5",
                                            "Ch29_1",
                                            "Eubac_1","Eubac_2","Eubac_3","Eubac_4",
                                            "Oscillo_1","Oscillo_2","Oscillo_3","Oscillo_4", "Oscillo_5","Oscillo_6","Oscillo_7","Oscillo_8",
                                            "Pepto_1","Pepto_2","Pepto_3","Pepto_4","Pepto_5",
                                            "Tissi_1","Tissi_2",
                                            "Clostridia_1"),
                                 ordered =T)

# plot heatmap
# plot heatmap
ggplot(matrix_Pep_melt, aes(x = genome, y = merops_subfamily, fill = value)) +
  geom_tile()+
  geom_text(size = 2, aes(label = ifelse(value !=0, value, "")))+
  theme_bw() +
  xlab("") + ylab("merops subfamily")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=7),
        axis.text.y =element_text(color='black',size=6),
        legend.position = 'right',
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(limits=rev)+
  scale_fill_gradientn("counts",colours = c("white","#e0f3db","#7fcdbb","#1d91c0","#253494","#081d58"),
                       na.value = "white",
                       breaks = c(0, 2, 4, 6, 8, 10, 12))

ggsave("Peptidases/Peptidase_families_overview.png",  width = 18, height = 21, units = "cm", device = "png" )
ggsave("Peptidases/Peptidase_families_overview.pdf",  width = 18, height = 21, units = "cm", device = "pdf" )
ggsave("Peptidases/Peptidase_families_overview.tiff",  width = 18, height = 21, units = "cm", device = "tiff" )

# save for SI files
Pep_wide <- matrix_Pep_melt %>%
  pivot_wider(names_from = genome, values_from = value, values_fill = list(value = 0))

write.table(Pep_wide, "Peptidases/Peptidases_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)


# Signal P Heatmap ----
# Search different signal peptides 
# (signal peptides are short amino acid sequences, which control protein secretion and translocation)

#rename ID in signalp table
signalp$ID <- gsub(" #.*","", signalp$ID ) 

#combine signalp and merops tables
SPs <- full_join(signalp, merops, by = c("ID"="qseqid" ))

# filter only those rows in which there is data for signalp and merops predicition
SPs_prediction <- SPs[complete.cases(SPs$merops_id), ]
unique(SPs_prediction$Prediction)
# "OTHER" "SP" "LIPO"
SPs_families <- SPs_prediction[SPs_prediction$Prediction =="SP",]
SPs_families <- separate(SPs_families, ID, into = c("genome","cluster"), sep = "___")

SPs_Clostridia <- subset(SPs_families, 
                         SPs_families$genome %in% c("E4_d458_spades_bin_53_orig-contigs",
                                                              "G3_d458_megahit_bin_18_orig-contigs",
                                                              "G3_d458_spades_bin_18_indstr_refined-contigs",
                                                              "G3_d458_spades_bin_33_strict_refined-contigs",
                                                              "L3_d458_spades_bin_6_orig_refined-contigs",
                                                              "E4_d458_spades_bin_64_permissive_refined-contigs",
                                                              "L3_d458_spades_bin_64_permissive_refined-contigs",
                                                              "G3_d458_megahit_bin_21_strict_refined-contigs",
                                                              "L3_d458_spades_bin_34_indper_refined-contigs",
                                                              "Van_BES_megahit_bin_42_indstr_refined-contigs",
                                                              "G3_d458_spades_bin_19_indper_refined-contigs",
                                                              "Co_bin_102_permissive_refined-contigs",
                                                              "E4_d458_spades_bin_8_permissive_refined-contigs",
                                                              "G3_d458_megahit_bin_13_indstr_refined-contigs",
                                                              "L3_d458_megahit_bin_13_permissive_refined-contigs",
                                                              "L3_d458_spades_bin_21_orig_refined-contigs",
                                                              "L3_d458_spades_bin_8_indstr_refined-contigs",
                                                              "M4_d458_megahit_bin_49_permissive_refined-contigs",
                                                              "M4_d458_spades_bin_15_permissive_refined-contigs",
                                                              "Co_bin_88_orig_refined-contigs",
                                                              "L3_d458_megahit_bin_58_indper-contigs",
                                                              "E4_d458_megahit_bin_7_strict_refined-contigs",
                                                              "M4_d458_megahit_bin_14_orig_refined-contigs",
                                                              "M4_d458_spades_bin_28_indstr_refined-contigs",
                                                              "Co_bin_49_orig_refined-contigs",
                                                              "L3_d458_spades_bin_31_orig_refined-contigs"))

# create a matrix of all signalPs.
matrix_SP <- table(SPs_Clostridia$genome, SPs_Clostridia$merops_subfamily)
matrix_SP_melt <- melt(matrix_SP)
colnames(matrix_SP_melt) <- c("genome","merops_subfamily","value")

matrix_SP_melt$genome <- factor(matrix_SP_melt$genome, levels = c("E4_d458_spades_bin_53_orig-contigs",
                                                                  "G3_d458_megahit_bin_18_orig-contigs",
                                                                  "G3_d458_spades_bin_18_indstr_refined-contigs",
                                                                  "G3_d458_spades_bin_33_strict_refined-contigs",
                                                                  "L3_d458_spades_bin_6_orig_refined-contigs",
                                                                  "E4_d458_spades_bin_64_permissive_refined-contigs",
                                                                  "G3_d458_megahit_bin_21_strict_refined-contigs",
                                                                  "L3_d458_spades_bin_34_indper_refined-contigs",
                                                                  "Van_BES_megahit_bin_42_indstr_refined-contigs",
                                                                  "G3_d458_spades_bin_19_indper_refined-contigs",
                                                                  "Co_bin_102_permissive_refined-contigs",
                                                                  "E4_d458_spades_bin_8_permissive_refined-contigs",
                                                                  "L3_d458_megahit_bin_13_permissive_refined-contigs",
                                                                  "G3_d458_megahit_bin_13_indstr_refined-contigs",
                                                                  "L3_d458_spades_bin_8_indstr_refined-contigs",
                                                                  "M4_d458_spades_bin_15_permissive_refined-contigs",
                                                                  "L3_d458_spades_bin_21_orig_refined-contigs",
                                                                  "M4_d458_megahit_bin_49_permissive_refined-contigs",
                                                                  "Co_bin_88_orig_refined-contigs",
                                                                  "L3_d458_megahit_bin_58_indper-contigs",
                                                                  "E4_d458_megahit_bin_7_strict_refined-contigs",
                                                                  "M4_d458_megahit_bin_14_orig_refined-contigs",
                                                                  "M4_d458_spades_bin_28_indstr_refined-contigs",
                                                                  "Co_bin_49_orig_refined-contigs",
                                                                  "L3_d458_spades_bin_31_orig_refined-contigs",
                                                                  "L3_d458_spades_bin_64_permissive_refined-contigs"),
                                labels = c("Aceti_1","Aceti_2", "Aceti_3","Aceti_4", "Aceti_5",
                                           "Ch29_1",
                                           "Eubac_1","Eubac_2","Eubac_3","Eubac_4",
                                           "Oscillo_1","Oscillo_2","Oscillo_3","Oscillo_4", "Oscillo_5","Oscillo_6","Oscillo_7","Oscillo_8",
                                           "Pepto_1","Pepto_2","Pepto_3","Pepto_4","Pepto_5",
                                           "Tissi_1","Tissi_2",
                                           "Clostridia_1"),
                                ordered =T)

# plot heatmap for extracellular peptidases
ggplot(matrix_SP_melt, aes(x = genome, y = merops_subfamily, fill = value)) +
  geom_tile()+
  geom_text(size = 2, aes(label = ifelse(value !=0, value, "")))+
  theme_bw() +
  xlab("") + ylab("merops subfamily")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=7),
        axis.text.y =element_text(color='black',size=6),
        legend.position = 'right',
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(limits=rev)+
  scale_fill_gradientn("counts",colours = c("white","#e0f3db","#7fcdbb","#1d91c0"),
                       na.value = "white",
                       breaks = c(0, 2, 4, 6, 8, 10, 12))

ggsave("Peptidases/Peptidase_SPs_families_overview.png",  width = 18, height = 7, units = "cm", device = "png" )
ggsave("Peptidases/Peptidase_SPs_families_overview.pdf",  width = 18, height = 7, units = "cm", device = "pdf" )
ggsave("Peptidases/Peptidase_SPs_families_overview.tiff",  width = 18, height = 7, units = "cm", device = "tiff" )


# save for SI files
SP_wide <- matrix_SP_melt %>%
  pivot_wider(names_from = genome, values_from = value, values_fill = list(value = 0))

write.table(SP_wide, "Peptidases/Peptidases_SP_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

