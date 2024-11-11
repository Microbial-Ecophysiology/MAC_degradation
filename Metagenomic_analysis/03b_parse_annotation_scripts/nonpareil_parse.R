#!/usr/bin/env Rscript

# parse assembly metadata for NCBI genomes

### setting up environment ####

# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
  "Nonpareil"
)

# Function to check if packages are installed
is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

# If not all packages are available
if(any(!is.installed(package.list))) {
  cat("Not all required packages are available.\n")
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
    c("-d", "--dir"),
    type = "character",
    default = NULL,
    help = "directory with npo files (all such files will be processed)",
    metavar = "character"
  ),
  make_option(
    c("-s", "--summary"),
    type = "character",
    default = NULL,
    help = "tabular summary output",
    metavar = "character"
  ),
  make_option(
    c("-p", "--plot"),
    type = "character",
    default = NULL,
    help = "nonpareil curves output",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dir) | is.null(opt$summary) | is.null(opt$plot)) {
  print_help(opt_parser)
  stop("All parameters are mandatory.\n", call. = FALSE)
}


### plot nonpareil curves ####

npo_files <- list.files(path = opt$dir, pattern = ".npo")
npo_cols <- rainbow(length(npo_files))

pdf(opt$plot, width = 7, height = 7)
nps <- Nonpareil.set(
  file.path(opt$dir, npo_files), 
  col = npo_cols,
  labels = gsub("_nonpareil.npo", "", basename(npo_files), fixed = T),
  plot.opts = list(plot.observed = F)
)
dev.off()


### summary statistics


nonpareil_summary <- data.frame(
  current_coverage = summary(nps)[, "C"] * 100,
  current_seq_effort_Gbp = sapply(nps@np.curves, function(x) x@LR)/1e9,
  near_complete_seq_effort_Gpb = summary(nps)[,"LRstar"]/1e9,
  coverage_10Gbp = sapply(nps$np.curves, predict, 10e9)
)
write.table(
  nonpareil_summary,
  opt$summary,
  quote = F,
  sep = "\t"
)
