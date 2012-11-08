library(tools)
library(Biostrings)
library(parallel)

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
src.project.file('src','domainArchitectureSimilarity.R')
src.project.file('src','loadUniprotKBEntries.R')

# Usage:
print( "Usage: Rscript measureDistances.R path/2/proteins.fasta path/2/accession_per_line.txt path/2/protein_function_annotation_matrix.r_serialized domain_weights_table.tbl path/2/output_name cores.2.use[default=all]")
print( "WARNING: Make sure the complete sequence names in the FASTA file are exactly the same as in the accession_per_line file!" )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[1]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[1]]) )

# Read Accessions to process:
accs <- as.character( read.table( trailing.args[[2]] )$V1 )
print( paste("Read", length(accs), "accessions to compute partial distance matrices for. Accessions-File is:", trailing.args[[2]]) )

# Load Uniprot function (InterPro and GO) Annotations:
f <- file( trailing.args[[3]], "r" )
annos <- unserialize( f )
close( f )
print( paste("Read", ncol(annos), "function annotations from", trailing.args[[3]]) )

# Load domain weights table:
dom.weights <- read.table( trailing.args[[4]] )
print( paste("Read", nrow(dom.weights), "domain weights from table", trailing.args[[4]]) )

# Write output to files beginning with..
out.head <- trailing.args[[ 5 ]]
print( paste("Writing output to files whose names will start with", out.head) )

# How many cores to use:
if( is.null(trailing.args[[ 6 ]]) ) {
  # All available cores:
  options( 'mc.cores'=detectCores() )
} else {
  # As specified by the sixth command line argument:
  options( 'mc.cores'=trailing.args[[ 6 ]] )
}

print( "Starting computation" )
print( paste( "Will be using", options('mc.cores'), "cores" ) )

# Partial Sequence Distances:
part.seq.dists.file <- paste( out.head, "_seq_dists.tbl", sep="" )
if( ! file.exists( part.seq.dists.file ) ) {
  write.table( partialSequenceDistances( aa.seqs, accs ),
    part.seq.dists.file
  )
  print( "computed sequence distances" )
} else {
  print( paste( "Partial sequence distances file already exists. Did not compute them again.", 
    part.seq.dists.file)
  )
}

# Partial Domain Architecture Distances:
part.das.dists.file <- paste( out.head, "_das_dists.tbl", sep="" )
if( ! file.exists( part.das.dists.file ) ) {
  write.table( partialDomainArchitectureDistances( annos, dom.weights, accs ),
    part.das.dists.file
  )
  print( "computed das distances" )
} else {
  print( paste( "Partial domain architecture sequence distances file already exists. Did not compute them again.",
      part.das.dists.file)
  )
}

# Partial boolean 'shared function' matrix:
part.shrd.funct.file <- paste( out.head, "_shrd_func.tbl", sep="" )
if( file.exists( part.shrd.funct.file ) ) {
  write.table( sharedFunction( annos, accs ), 
    part.shrd.funct.file
  )
  print( "computed shared functions" )
} else {
  print( paste("Partial shared function file already exists. Did not compute them again.",
      part.shrd.funct.file)
  )
}

print( "DONE" )
