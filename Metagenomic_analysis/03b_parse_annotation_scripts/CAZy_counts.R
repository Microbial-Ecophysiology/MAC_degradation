# Script to investigate CAZY counts

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

## dbCAN ----
dbCAN <- read.csv("CAZymes/concat_dbCAN_annotation.txt", 
                  sep="\t", 
                  header=TRUE, 
                  fill=TRUE, 
                  quote="",
                  col.names=c("Gene.ID","dbCAN_EC","dbCAN_HMMER","dbCAN_eCAMI","dbCAN_DIAMOND","dbCAN_no.Tools"))
dbCAN <- subset(dbCAN, dbCAN_EC != "EC#")


# combine table with Cazymes in first few columns
CAZy <- full_join(dbCAN, kegg, by = c("Gene.ID"="qseqid")) %>%
  full_join(., NR, by = c("Gene.ID"="qseqid")) %>%
  separate(col = Gene.ID, into = c("genome","cluster"), sep = "___")


CAZy_Clostridia <- subset(CAZy, 
                          CAZy$genome %in% c("E4_d458_spades_bin_53_orig-contigs",
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



# remove all rows from dataframe which don't have any annotation through NR and were only annotated by one of the tools
CAZy_short <- CAZy_Clostridia[complete.cases(CAZy_Clostridia[, 3:7]) & complete.cases(CAZy_Clostridia[, 24:41]) & CAZy_Clostridia$dbCAN_no.Tools > 1 , ]
# 855 CAZymes remain

# Remove everything in parentheses in the dbCAN_HMMER column
CAZy_short$dbCAN_HMMER <- gsub("\\(.*?\\)", "", CAZy_short$dbCAN_HMMER)

# Function to find the common annotation across relevant columns
find_common_annotation <- function(row) {
  # Replace "-" with NA to handle missing values
  row <- gsub("-", NA, row)
  
  # Remove NA values from the row (treated as missing)
  row <- row[!is.na(row)]
  
  # Split annotations by "+" in each non-NA column
  all_annotations <- unlist(lapply(row, function(x) strsplit(x, "\\+")))
  
  # Use table to count the occurrences of each annotation
  annotation_counts <- table(all_annotations)
  
  # Keep annotations that appear at least twice
  common_annotations <- names(annotation_counts[annotation_counts >= 2])
  
  # If common annotations exist, return them
  if(length(common_annotations) > 0) {
    return(paste(common_annotations, collapse = "+"))
  } else {
    return(NA)  # No common annotations found, meaning different annotators give different annotations
  }
}

# Apply the function row-wise to the specific columns (dbCAN_HMMER, annotator2, annotator3)
CAZy_short$combined_annotation <- apply(CAZy_short[, c("dbCAN_HMMER", "dbCAN_eCAMI", "dbCAN_DIAMOND")], 1, find_common_annotation)

# remove rows with NA in combined_annotation after manually checking their annotation through NR and kegg. Remove only if no annotation as CAZyme
CAZy_short <- CAZy_short[!is.na(CAZy_short$combined_annotation), ]


# subset dataframe to genome, combined_annotation and dbCAN_no.Tools
CAZy_subset <- CAZy_short[,c("genome", "combined_annotation","dbCAN_no.Tools") ]

# combine single CAZy families
CAZy_summed <- CAZy_subset %>%
  group_by(genome, combined_annotation) %>%
  summarize(count = n(), .groups = 'drop')

unique(CAZy_summed$genome)

CAZy_summed$genome <- factor(
  CAZy_summed$genome, 
  levels = c("E4_d458_spades_bin_53_orig-contigs",
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
ggplot(CAZy_summed, aes(x = genome, y = combined_annotation, fill = count)) +
  geom_tile()+
  geom_text(size = 2, aes(label = ifelse(count !=0, count, "")))+
  theme_bw() +
  xlab("") + ylab("CAZy family")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=7),
        axis.text.y =element_text(color='black',size=6),
        legend.position = 'right',
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(limits=rev)+
  scale_fill_gradientn("counts",colours = c("#e0f3db","#7fcdbb","#1d91c0","#253494","#081d58"),
                       na.value = "white",
                       breaks = c(0, 2, 4, 6, 8, 10, 12))

ggsave("CAZymes/CAZy_families_overview.png",  width = 18, height = 21, units = "cm", device = "png" )
ggsave("CAZymes/CAZy_families_overview.pdf",  width = 18, height = 21, units = "cm", device = "pdf" )
ggsave("CAZymes/CAZy_families_overview.tiff",  width = 18, height = 21, units = "cm", device = "tiff" )


# save for SI files
CAZy_wide <- CAZy_summed %>%
  pivot_wider(names_from = genome, values_from = count, values_fill = list(count = 0))

write.table(CAZy_wide, "CAZymes/CAZy_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)


