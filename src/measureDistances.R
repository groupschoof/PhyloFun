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
print( "Usage: Rscript measureDistances.R path/2/proteins.fasta path/2/protein_function_annotation_matrix.r_serialized domain_weights_table.tbl path/2/blast_out.tbl go_term_id_per_line.tbl cores.2.use path/2/output_dir")
print( "WARNING: Make sure the _complete_ sequence names in the FASTA file are _exactly the same_ as in the function annotation file!" )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[ 1 ]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[ 1 ]]) )

# Load Uniprot function (InterPro and GO) Annotations:
f <- file( trailing.args[[ 2 ]], "r" )
annos <- unserialize( f )
close( f )
print( paste("Read", ncol(annos), "function annotations from", trailing.args[[ 2 ]]) )

# Load domain weights table:
dom.weights <- read.table( trailing.args[[ 3 ]] )
print( paste("Read", nrow(dom.weights), "domain weights from table", trailing.args[[ 3 ]]) )

# Read tabular Blast Output
blast.out <- read.table( trailing.args[[ 4 ]] )
print( paste("Read", nrow(blast.out), "blast results (query-hit-pairs) from table", trailing.args[[ 4 ]]) )

# Read GO terms to measure distances for:
go.terms <- as.character( read.table( trailing.args[[ 5 ]] )[ , 1 ] )
print( paste("Read", length( go.terms ),
    "GO terms to measure mutation probabilities",
    "depending on sequence and domain architecture distances for from file",
    trailing.args[[ 5 ]])
)

# How many cores to use:
options( 'mc.cores'=trailing.args[[ 6 ]] )
print( paste( "Will be using", options('mc.cores'), "parallel threads." ) )

# Path to output directory:
path.2.output.dir <- as.character( trailing.args[[ 7 ]] )
print( paste( "Will write output to files in directory", path.2.output.dir ) )

# Begin
print( "Starting computation" )

for ( go.term in go.terms ) {
  # Output path
  go.term.out.path <- paste( path.2.output.dir, "/", go.term, "_", sep="" )

  # Measure distances
  go.dists <- measureDistances( go.term, annos, blast.out,
    aa.seqs, dom.weights, lapply.funk=mclapply
  )
  write.table( go.dists,
    file=paste( go.term.out.path, "distances.tbl", sep="" )
  )
  
  # mutation probability depending on sequence distance
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
}

print( "DONE" )
