library(tools)
library(Biostrings)
library(parallel)
library(RCurl)

# In R sourcing other files is not trivial, unfortunately.
# WARNING:
# This method ONLY works for project files in depth one sub dirs!
project.file.path <- function(...) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dir <- dirname(file_path_as_absolute(script.name))
  project.dir <- sub(basename(script.dir),'',script.dir)
  normalizePath(file.path(project.dir,...))
}
src.project.file <- function(...) {
  source(project.file.path(...))
}
# Helper functions:
src.project.file( "src", "phyloFunTools.R" )

# Hail User:
print( "Usage: Rscript runPhyloFun.R path/2/query_proteins.fasta path/2/jackhmmer_results.tbl path/2/interproscan_results.tsv" )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[1]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[1]]) )

# Parse Jackhmmer results:
jr <- parseJackhmmerTable( 
  scan( file=trailing.args[[2]], what=character(), sep="\n" )
)

# for each protein, do
# 
# get jackhmmer results for that protein
# download their sequences
# us <- lapply( accs, function(a) uniprotkb.url( a, frmt='fasta' ) )
# writeLines( getURL( us ), file )

