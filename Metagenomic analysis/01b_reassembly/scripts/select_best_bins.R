#!/usr/bin/env Rscript

### setting up environment ####
# Check if packages are installed
package.list <- c(
  "crayon",
  "optparse",
  "config",
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
    c("-s", "--checkm"),
    type = "character",
    default = NULL,
    help = "checkM stats",
    metavar = "character"
  ),
  make_option(
    c("-a", "--sfa"),
    type = "numeric",
    default = 5,
    help = "scoring factor A",
    metavar = "number"
  ),
  make_option(
    c("-b", "--sfb"),
    type = "numeric",
    default = 1.2,
    help = "scoring factor B",
    metavar = "number"
  ),
  make_option(
    c("-c", "--completeness"),
    type = "numeric",
    default = 50,
    help = "minimum completeness",
    metavar = "number"
  ),
  make_option(
    c("-r", "--redundancy"),
    type = "numeric",
    default = 10,
    help = "maximum contamination",
    metavar = "number"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "checkM stats of best bins",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$checkm)) {
  print_help(opt_parser)
  stop("CheckM stats have to be provided.\n", call. = FALSE)
}


### Determine best bins based on scoring ####

A <- opt$sfa
B <- opt$sfb

checkm <- read.table(opt$checkm, h = T, sep = "\t") %>% 
  
  # remove any bins with too low completeness and too high contamination from further scoring and selection
  filter(completeness >= opt$completeness & contamination <= opt$redundancy) %>% 
  
  # define score weighted by A and B
  # unlike the default metawrap scoring method, contamination is weighted non-linearly
  # other than that the scale of the score for contamination less than 10 is comparable
  # for contamination cutoff larger than 10 a different scoring algorithm is advisable
  # for contamination cutoff smaller than 10, A and B should be adapted
  mutate(
    score = completeness + A * (100 - 1/A * (contamination/B)^2),
    group = gsub("\\.orig|\\.permissive|\\.strict|\\.indstr|\\.indper|_aug$|_ori$", "", bin)
  ) %>% 
  group_by(group) %>% 
  
  # for each bin, select the best bin by highest score and break ties by larger N50
  arrange(desc(score), desc(N50), .by_group = T) %>% 
  filter(row_number() == 1) %>% 
  
  # some maintenance for nice output
  ungroup() %>% 
  select(-group)

write.table(
  checkm,
  opt$output,
  quote = F,
  sep = "\t",
  row.names = F
)
