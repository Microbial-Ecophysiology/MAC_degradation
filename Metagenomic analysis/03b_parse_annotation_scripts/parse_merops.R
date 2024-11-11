#!/usr/bin/env Rscript
# parse assembly metadata for NCBI genomes
### setting up environment ####
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "tidyverse",
  "data.table"
)
silent <- suppressMessages(lapply(package.list, function(X) {require(X, character.only = TRUE)}))
rm(silent)
### Reading command line options ####
# define command line options
option_list <- list(
  make_option(
    c("-f", "--filter"),
    type = "character",
    default = NULL,
    help = "merops scan blast output",
    metavar = "character"
  ),
  make_option(
    c("-p", "--pepunit"),
    type = "character",
    default = NULL,
    help = "merops pepunit blast output",
    metavar = "character"
  ),
  make_option(
    c("-m", "--merops"),
    type = "character",
    default = NULL,
    help = "merops metadata table",
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
    help = "path and filebase for output files",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$filter) | is.null(opt$pepunit) | is.null(opt$selfout) | is.null(opt$output) | is.null(opt$merops)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}
### filter merops blast results ####
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
# merops scan blast output
scanout <- read.table(
  opt$filter,
  h = F,
  sep = "\t",
  stringsAsFactors = F,
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
)
# merops pepunit blast output
# only considering hits present in scan
pepunitout <- read.table(
  opt$pepunit,
  h = F,
  sep = "\t",
  stringsAsFactors = F,
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
) %>% 
  filter(qseqid %in% scanout$qseqid) %>%
  mutate(
    bsr = round(.$bitscore/selfout[.$qseqid, "bitscore"], 4) %>% pull(1)
  ) %>% 
  filter(
    bsr >= opt$cutoff
  )
### parse merops metadata ####
merops <- read.table(
  opt$merops,
  h = T,
  sep = "\t",
  comment.char = "",
  quote = "",
  stringsAsFactors = F
)
output <- data.frame(
  pepunitout,
  merops[match(pepunitout$sseqid, merops$id), -1]
)
output_best <- output %>%
  group_by(qseqid) %>%
  arrange(desc(bitscore), evalue) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(qseqid) %>%
  as.data.frame()
write.table(
  output,
  paste0(opt$output, ".txt"),
  quote = F,
  sep = "\t",
  row.names = F
)
write.table(
  output_best,
  paste0(opt$output, "_best.txt"),
  quote = F,
  sep = "\t",
  row.names = F
)
