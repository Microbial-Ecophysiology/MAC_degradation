# Script to check Transporter in MAGs

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Manuscript_MAC/Metagenomics/annotation/")

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

## TCDB
TCDB <- read.csv("Transporter/concat_TCDB_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "TCDB_sseqid","TCDB_pident","TCDB_length","TCDB_mismatch","TCDB_gapopen","TCDB_qstart",
                             "TCDB_qend","TCDB_sstart","TCDB_send","TCDB_evalue","TCDB_bitscore","TCDB_bsr","TCDB_tc_id","TCDB_tc_fam","TCDB_family","TCDB_go","TCDB_pfam","TCDB_chebi"))
TCDB <- subset(TCDB, TCDB_sseqid != "sseqid")


# heatmap of Transporters
# combine table with Transporters in first few columns
Transporters_table <- full_join(TCDB, kegg, by = c("qseqid")) %>%
  full_join(., NR, by = c("qseqid")) %>%
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")

# remove all rows from dataframe which don't have any annotation through NR (or kegg)as NR has annotation for all
Transporters_short <- Transporters_table[complete.cases(Transporters_table[, 3:17]) & complete.cases(Transporters_table[, 37]), ]
# 5402 Transporters remain

Transporters_Clostridia <- subset(Transporters_short, 
                                  Transporters_short$genome %in% c("E4_d458_spades_bin_53_orig-contigs",
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

# 3520 Transporter


## Create a matrix of TCDB families against genomes ####
matrix_TCDB_short <- table(Transporters_Clostridia$genome, Transporters_Clostridia$TCDB_tc_fam, Transporters_Clostridia$TCDB_family)
# melt matrix into tabular view
matrix_TCDB_melt_short <- melt(matrix_TCDB_short)
colnames(matrix_TCDB_melt_short) <- c("genome","tc_fam","TCDB_family","value")


matrix_TCDB_melt_short$genome <- factor(matrix_TCDB_melt_short$genome, levels = c("E4_d458_spades_bin_53_orig-contigs",
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

## plot heatmap for all TCDB families ####
ggplot(matrix_TCDB_melt_short, aes(x = genome, y = TCDB_family, fill = value)) +
  geom_tile()+
  geom_text(size = 2, aes(label = ifelse(value !=0, value, "")))+
  theme_bw() +
  xlab("") + ylab("TCDB family")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=7),
        axis.text.y =element_text(color='black',size=4),
        legend.position = 'right',
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(limits=rev)+
  scale_fill_gradientn("counts",colours = c("white","#e0f3db","#7fcdbb","#1d91c0","#253494","#081d58"),
                       na.value = "white",
                       breaks = c(0,20,40,60,80))

ggsave("Transporter/TCDB_families_overview.png", plot = last_plot(), width = 25, height = 30, units = "cm", device = "png" )


# save for SI files
TCDB_wide <- matrix_TCDB_melt_short %>%
  pivot_wider(names_from = genome, values_from = value, values_fill = list(value = 0))
df_filtered <- TCDB_wide[rowSums(TCDB_wide[, 3:28]) != 0, ]

write.table(df_filtered, "Transporter/Transporters_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)
