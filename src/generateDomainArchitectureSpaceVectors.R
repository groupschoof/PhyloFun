library(tools)

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
print( "Usage: Rscript generateDomainArchitectureSpaceVectors.R path/2/domain_weights_table.tbl path/2/function_annotations.tbl path/2/domain_architecture_space_vectors_list.r_serialized")
print( "WARNING: Make sure the complete sequence names in the FASTA file are exactly the same as in the accession_per_line file!" )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read domain weights database:
dom.weights <- read.table( trailing.args[[1]] )

# Read annotations:
f <- file( trailing.args[[2]], "r")
annos <- unserialize(f)
close(f)

# Generate vector space model:
vsm <- constructVectorSpaceModel( annos )

# Vectors in VSM for each Protein:
dom.vects <- generateDomainArchitectureSpaceVectors( vsm, annos, dom.weights )

# Write out result:
f <- file( trailing.args[[3]], "w" )
serialize( dom.vects, file=f )
close(f)
