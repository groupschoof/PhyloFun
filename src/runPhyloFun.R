library(tools)
library(Biostrings)
library(RCurl)
library(stringr)

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
print( "Usage: Rscript runPhyloFun.R path/2/query_proteins.fasta path/2/jackhmmer_results.tbl [ query-accession-regex='^\\s*(\\S+)\\s' ]")

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[1]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[1]]) )

# Parse Jackhmmer results:
jr <- parseJackhmmerTable( 
  scan( file=trailing.args[[2]], what=character(), sep="\n" )
)
print( paste( "Parsed JACKHMMER result table. Got", nrow(jr), "query-hit-pairs" ) )

# Parse accessions as String between '>' and the first blank character:
query.acc.regex <- if ( length(trailing.args) == 3 ) trailing.args[[3]] else '^\\s*(\\S+)\\s' 
accs <- unlist( lapply( names(aa.seqs), function(n) str_match( n, query.acc.regex )[[ 1, 2 ]] ) )
print( paste( "Parsed the query proteins' accessions as text matching",
  query.acc.regex, "after the '>' character in file", trailing.args[[ 1 ]],
  "If your query accessions are not matching, PhyloFun will fail to find their accessions in the JACKHMMER results!" )
)

# For each query protein, do:
for ( acc in accs ) {
  hits <- jr[ which( jr[ , 'query.name' ] == acc ), , drop=F ]
  if ( nrow( hits ) > 0 ) {
    hit.accs <- hits[ , 'hit.name' ]
    hit.uniprot.docs <- getURL( sapply( hit.accs, function( a ) uniprotkb.url( a ) ) )
    hit.seqs <- retrieveSequences( hit.uniprot.docs )
  }
}
# 
# get jackhmmer results for that protein
# download their sequences
# us <- lapply( accs, function(a) uniprotkb.url( a, frmt='fasta' ) )
# writeLines( getURL( us ), file )

