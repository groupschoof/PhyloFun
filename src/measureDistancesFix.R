library(tools)
library(Biostrings)
library(phangorn)
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
src.project.file('src','domainArchitectureSimilarityRedis.R')
src.project.file('src','loadUniprotKBEntries.R')
src.project.file('src','measureDistanceFunctions.R')

# Usage:
print( "Usage: Rscript measureDistancesFix.R path/2/distances.tbl go_term_id path/2/output_dir")

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read Distances Tbl:
go.dists <- read.table( trailing.args[[ 1 ]] )

# GO term:
go.term <- as.character( trailing.args[[ 2 ]] )

# Path to output directory:
path.2.output.dir <- as.character( trailing.args[[ 3 ]] )
print( paste( "Will write output to files in directory", path.2.output.dir ) )

go.term.out.path <- paste( path.2.output.dir, "/", go.term, "_", sep="" )

# Begin
print( "Starting computation" )

write.table( 
  pMutationMinMaxParentValues(
    mutationProbabilityDistribution( go.dists, "Sequence.Distance" )
    , "p.mutation|Sequence.Distance"
  ),
  file=paste( go.term.out.path, "p_mut_seq_dist.tbl", sep="" )
)

# mutation probability depending on domain architecture distance
write.table(
  pMutationMinMaxParentValues(
    mutationProbabilityDistribution( go.dists,
      "Domain.Architecture.Distance"
    ),
    "p.mutation|Domain.Architecture.Distance"
  ), 
  file=paste( go.term.out.path, "p_mut_das_dist.tbl", sep="" )
)

# mutation probability depending on euclidean distance to origin
write.table(
  pMutationMinMaxParentValues(
    mutationProbabilityDistribution( go.dists,
      "Euclidean.Distance.To.Origin"
    ),
    "p.mutation|Euclidean.Distance.To.Origin"
  ),
  file=paste( go.term.out.path, "p_mut_seq_das_dist.tbl", sep="" )
)

print( "DONE" )
