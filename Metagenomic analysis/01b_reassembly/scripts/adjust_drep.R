#!/usr/bin/env Rscript

# automatically adjust drep output prior to bin reassembly

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "tidyverse",
  "scales",
  "igraph",
  "reshape",
  "data.table"
)

# Function to check if packages are installed
is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

# If not all packages are available
if(any(!is.installed(package.list))) {
  cat("Not all required packages are available. They will now be installed.\n")
  
  # give the user the chance to abort manually
  Sys.sleep(20)
  
  # then install packages
  for(i in which(!is.installed(package.list))) {
    suppressMessages(install.packages(package.list[i], repos = "http://cran.us.r-project.org"))
  }
}

# Break the script if the package installation was unsuccessful
if(any(!is.installed(package.list))) {
  cat(
    paste0(
      "Unable to install all required packages.\nPlease install ",
      paste0(package.list[!is.installed(package.list)], collapse = ", "),
      " manually."
    )
  )
  break
}

# Load packages
cat("Loading libraries...")
silent <- suppressMessages(lapply(package.list, function(X) {require(X, character.only = TRUE)}))
rm(silent)
cat(" done\n")

# Some functions for message output
msg <- function(X){
  cat(crayon::white(paste0("[",format(Sys.time(), "%T"), "]")), X)
}
msg_sub <- function(X){
  cat(crayon::white(paste0("  [",format(Sys.time(), "%T"), "]")), X)
}


### Reading command line options ####

# define command line options
option_list <- list(
  make_option(
    c("-c", "--cdb"),
    type = "character",
    default = NULL,
    help = "drep Cdb output file",
    metavar = "character"
  ),
  make_option(
    c("-w", "--wdb"),
    type = "character",
    default = NULL,
    help = "drep Wdb output file",
    metavar = "character"
  ),
  make_option(
    c("-g", "--genome_info"),
    type = "character",
    default = NULL,
    help = "drep genome information output",
    metavar = "character"
  ),
  make_option(
    c("-f", "--fastani"),
    type = "character",
    default = NULL,
    help = "fastani check summary file",
    metavar = "character"
  ),
  make_option(
    c("-q", "--quant"),
    type = "character",
    default = NULL,
    help = "bin coverage file list",
    metavar = "character"
  ),
  make_option(
    c("-a", "--ani"),
    type = "double",
    default = 95,
    help = "ANI similarity cut-off (default: 95)",
    metavar = "number"
  ),
  make_option(
    c("-A", "--fraction"),
    type = "double",
    default = 0.5,
    help = "aligned fraction cut-off (default: 0.5)",
    metavar = "number"
  ),
  make_option(
    c("-d", "--diff_cov"),
    type = "double",
    default = 0.8,
    help = "correlation between bins (default: 0.8)",
    metavar = "number"
  ),
  make_option(
    c("-C", "--cdb_out"),
    type = "character",
    default = NULL,
    help = "updated non-redundant Cdb file",
    metavar = "character"
  ),
  make_option(
    c("-W", "--wdb_out"),
    type = "character",
    default = NULL,
    help = "updated non-redundant Wdb file",
    metavar = "character"
  ),
  make_option(
    c("-r", "--workspace"),
    type = "character",
    default = NULL,
    help = "R workspace to inspect intermediate results",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$cdb) | is.null(opt$wdb) | is.null(opt$fastani) | is.null(opt$quant) | 
    is.null(opt$cdb_out) | is.null(opt$wdb_out) | is.null(opt$workspace)) {
  print_help(opt_parser)
  stop("Provide all input and output file names.\n", call. = FALSE)
}


### automatically adjust drep output prior to bin reassembly

# opt$cdb = "Bin_dereplication/drep_out/data_tables/Cdb.csv"
# opt$wdb = "Bin_dereplication/drep_out/data_tables/Wdb.csv"
# opt$genome_info = "Bin_dereplication/drep_out/data_tables/genomeInfo.csv"
# opt$fastani = "Bin_dereplication/fastani_check/check_summary.txt"
# opt$quant = "Co_refinement/Co_bins_quant.txt Single_refinement/Mg_ESP_St14_12m/Mg_ESP_St14_12m_bins_quant.txt Single_refinement/Mg_ESP_St14_18m/Mg_ESP_St14_18m_bins_quant.txt Single_refinement/Mg_ESP_St14_22m/Mg_ESP_St14_22m_bins_quant.txt Single_refinement/Mg_ESP_St14_27m/Mg_ESP_St14_27m_bins_quant.txt Single_refinement/Mg_ESP_St14_37m/Mg_ESP_St14_37m_bins_quant.txt Single_refinement/Mg_ESP_St14_58m/Mg_ESP_St14_58m_bins_quant.txt Single_refinement/Mg_ESP_St18_15_20m/Mg_ESP_St18_15_20m_bins_quant.txt Single_refinement/Mg_ESP_St18_30_50m/Mg_ESP_St18_30_50m_bins_quant.txt Single_refinement/Mg_ESP_St18_70_88m/Mg_ESP_St18_70_88m_bins_quant.txt Single_refinement/Mg_ESP_St26_25m/Mg_ESP_St26_25m_bins_quant.txt Single_refinement/Mg_ESP_St26_30m/Mg_ESP_St26_30m_bins_quant.txt Single_refinement/Mg_ESP_St26_35_60m/Mg_ESP_St26_35_60m_bins_quant.txt Single_refinement/Mg_ESP_St26_80_120m/Mg_ESP_St26_80_120m_bins_quant.txt Single_refinement/Mg_ESP_St31_150m/Mg_ESP_St31_150m_bins_quant.txt Single_refinement/Mg_ESP_St31_200m/Mg_ESP_St31_200m_bins_quant.txt Single_refinement/Mg_ESP_St31_70_80m/Mg_ESP_St31_70_80m_bins_quant.txt Single_refinement/Mg_ESP_St31_90_120m/Mg_ESP_St31_90_120m_bins_quant.txt Single_refinement/Mg_ESP_St39_600_1000m/Mg_ESP_St39_600_1000m_bins_quant.txt"
# opt$cdb_out = "Bin_dereplication/drep_adj/Cdb_adjusted.csv"
# opt$wdb_out = "Bin_dereplication/drep_adj/Wdb_adjusted.csv"
# opt$workspace = "Bin_dereplication/drep_adj/drep_adjustment.Rdata"

# read coverage patterns and correlate bins within each primary cluster (differential coverage)
cov_list <- strsplit(opt$quant, " ")[[1]]
mean_cov <- map_dfr(
  cov_list,
  function(x) {
    tmp <- read.table(
      x,
      h = T, 
      sep = "\t",
      check.names = F
    )
    rownames(tmp) <- paste0(gsub("_bins_quant.txt", "", basename(x)), "_", tmp$Genome, ".fa")
    colnames(tmp) <- gsub(" Mean", "", colnames(tmp))
    tmp <- tmp[, gsub("_bins_quant.txt", "", basename(cov_list))[-1]]
    return(tmp)
  }
)

# read dRep output
wdb_drep <- read.table(opt$wdb, h = T, sep = ",")
cdb_drep <- read.table(opt$cdb, h = T, sep = ",")
genome_info <- read.table(opt$genome_info, h = T, sep = ",")

# read and filter fastani output
fastani_out <- read.table(
  opt$fastani,
  h = F, 
  sep = "\t",
  col.names = c("primary_cluster", "X1", "X2", "ani", "match", "all")
) %>% 
  mutate(
    af = match/all,
    X1 = basename(X1),
    X2 = basename(X2)
  ) %>% 
  # remove self matches
  filter(X1 != X2) %>% 
  # only look at aligned fraction relative to smaller genome
  mutate(
    X1_length = genome_info[match(.$X1, genome_info$genome), "length"],
    X2_length = genome_info[match(.$X2, genome_info$genome), "length"]
  ) %>% 
  filter(
    X1_length <= X2_length
  ) %>% 
  # filter to ANI 95% and aligned fraction 0.5
  filter(
    ani >= opt$ani,
    af >= opt$fraction
  ) %>% 
  # make sure that genomes within a primary cluster are not from the same assembly
  mutate(
    X1_assembly = gsub("_bin.*", "", X1),
    X2_assembly = gsub("_bin.*", "", X2),
  ) %>% 
  filter(X1_assembly != X2_assembly) %>% 
  # calculate correlation between coverage profiles for each comparison
  mutate(
    cov_cor = sapply(
      1:nrow(.),
      function(i) {
        cor(t(mean_cov[.$X1[i], ]), t(mean_cov[.$X2[i], ]))
      }
    )
  ) %>% 
  # filter to comparisons with >= 0.8 correlation
  filter(cov_cor >= opt$diff_cov)

# create graph for primary clusters
# in case there are sub clusters in the fastani output within a primary cluster, those are easier to detect this way
ani_edge <- fastani_out[, c("X1", "X2", "ani")]
colnames(ani_edge)[3] <- "weight"
ani_g <- graph_from_data_frame(ani_edge, directed = F)
length(E(ani_g))
ani_groups <- data.frame(nr_cluster = igraph::components(ani_g)$membership) %>% 
  rownames_to_column(var = "genome") %>% 
  mutate(
    primary_cluster = cdb_drep[match(.$genome, cdb_drep$genome), "primary_cluster"],
    secondary_cluster = cdb_drep[match(.$genome, cdb_drep$genome), "secondary_cluster"],
    score = wdb_drep[match(.$genome, wdb_drep$genome), "score"],
    completeness = genome_info[match(.$genome, genome_info$genome), "completeness"],
    contamination = genome_info[match(.$genome, genome_info$genome), "contamination"],
    length = genome_info[match(.$genome, genome_info$genome), "length"],
    N50 = genome_info[match(.$genome, genome_info$genome), "N50"]
  ) %>% 
  # very naive approach: just select the genome with the best score from each graph cluster
  # this may have to be adjusted later
  group_by(nr_cluster) %>% 
  arrange(desc(score)) %>% 
  mutate(
    winning = ifelse(row_number() == 1, "keep", "redundant"),
    secondary_cluster_adjusted = ifelse(winning == "keep", secondary_cluster, secondary_cluster[winning == "keep"])
  ) %>%
  ungroup() %>% 
  arrange(primary_cluster, nr_cluster) %>% 
  as.data.frame()

# write output
cdb_out <- rbind(
  cdb_drep[!cdb_drep$genome %in% ani_groups$genome, ],
  data.frame(
    genome = ani_groups$genome,
    secondary_cluster = ani_groups$secondary_cluster_adjusted,
    threshold = (100 - opt$ani)/100,
    cluster_method = "post_drep_adjustment",
    comparison_algorithm = "fastANI",
    primary_cluster = ani_groups$primary_cluster
  )
)
cdb_out <- cdb_out[match(cdb_drep$genome, cdb_out$genome), ]
# also adjust secondary cluster ID for members of redunant (previously) winning genome
for(i in which(ani_groups$winning == "redundant")) {
  if(sum(cdb_out$secondary_cluster == ani_groups$secondary_cluster[i]) > 0) {
    cdb_out$cluster_method[cdb_out$secondary_cluster == ani_groups$secondary_cluster[i]] <- "post_drep_adjustment"
    cdb_out$secondary_cluster[cdb_out$secondary_cluster == ani_groups$secondary_cluster[i]] <- ani_groups$secondary_cluster_adjusted[i]
  }
}
write.table(cdb_out, opt$cdb_out, quote = F, sep = ",", row.names = F)

wdb_out <- wdb_drep %>% filter(!genome %in% ani_groups$genome[ani_groups$winning == "redundant"])
write.table(wdb_out, opt$wdb_out, quote = F, sep = ",", row.names = F)

save.image(opt$workspace)
