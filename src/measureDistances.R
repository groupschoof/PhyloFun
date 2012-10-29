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
src.project.file('src','domainArchitectureDistances.R')

# Usage:
print( "Usage: Rscript measureDistances.R path/2/proteins.fasta path/2/accession_per_line.txt path/2/protein_function_annotation_matrix.r_serialized domain_weights_table.tbl path/2/output_name")
print( "WARNING: Make sure the complete sequence names in the FASTA file are exactly the same as in the accession_per_line file!" )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[1]] ), function(s) toString(s) )

# Read Accessions to process:
accs <- as.character( read.table( trailing.args[[2]] )$V1 )

# Load Uniprot function (InterPro and GO) Annotations:
f <- file( trailing.args[[3]], "r" )
annos <- unserialize( f )
close( f )

# Load domain weights table:
dom.weights <- read.table( trailing.args[[4]] )

# Write output to files beginning with..
out.head <- trailing.args[[ 5 ]]

# Partial Sequence Distances:
write.table( partialSequenceDistances( aa.seqs, accs ),
  paste( out.head, "_seq_dists.tbl", sep="" )
)

# Partial Domain Architecture Distances:
write.table( partialDomainArchitectureDistances( annos, dom.weights, accs ),
  paste( out.head, "_das_dists.tbl", sep="" )
)

# Partial boolean 'shared function' matrix:
write.table( sharedFunction( annos, accs ), 
  paste( out.head, "_shrd_func.tbl", sep="" )
)

print( "DONE" )
