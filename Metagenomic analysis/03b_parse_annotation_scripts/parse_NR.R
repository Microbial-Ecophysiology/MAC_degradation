#!/usr/bin/env Rscript

# parse assembly metadata for NCBI genomes

### setting up environment ####

package.list <- c(
  "crayon",
  "optparse",
  "config",
  "data.table",
  "tidyverse"
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
silent <- suppressMessages(lapply(package.list, function(X) {require(X, character.only = TRUE)}))
rm(silent)


### Reading command line options ####

# define command line options
option_list <- list(
  make_option(
    c("-b", "--blastout"),
    type = "character",
    default = NULL,
    help = "blast output",
    metavar = "character"
  ),
  make_option(
    c("-s", "--selfout"),
    type = "character",
    default = NULL,
    help = "blast output of self matches (only 1 hit per query)",
    metavar = "character"
  ),
  make_option(
    c("-c", "--cutoff"),
    type = "double",
    default = 0.4,
    help = "blast ratio score cut-off",
    metavar = "number"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "path and base of output files",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blastout) | is.null(opt$selfout) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("You need to provide the input and output file names.\n", call. = FALSE)
}


### parse blast output ####

# self matches for calculation of blast score ratio
selfout <- fread(
  opt$selfout,
  h = F,
  sep = "\t",
  col.names = c(
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
  )
)[, c("qseqid", "bitscore")]
setkey(selfout, "qseqid")

# kegg blast out
# append blast score ratio
# filter based on blast score ratio
blastout <- fread(
  opt$blastout,
  h = F,
  sep = "\t",
  col.names = c(
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "stitle",
    "staxids",
    "sscinames",
    "salltitles"
  )
) %>% 
  mutate(
    bsr = round(.$bitscore/selfout[.$qseqid, "bitscore"], 4) %>% pull(1)
  ) %>% 
  filter(
    bsr >= opt$cutoff
  )

# parse table with only best blast hit
blastout_best <- blastout %>%
  group_by(qseqid) %>%
  arrange(desc(bitscore), evalue) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(qseqid) %>%
  as.data.frame()

# write output
fwrite(
  blastout, 
  paste0(opt$output, ".txt"),
  quote = F, 
  sep = "\t",
  col.names = T,
  row.names = F
)

fwrite(
  blastout_best,
  paste0(opt$output, "_best.txt"),
  quote = F,
  sep = "\t",
  col.names = T,
  row.names = F
)
