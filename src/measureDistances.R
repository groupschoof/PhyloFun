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
print( "Usage: Rscript measureDistances.R path/2/proteins.fasta path/2/protein_function_annotation_matrix.r_serialized domain_weights_table.tbl batch.no batch.size cores.2.use REDIS-URL REDIS-Port")
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

# Read Batch Number and Batch Size:
batch.no <- as.integer( trailing.args[[ 4 ]] )
batch.size <- as.integer( trailing.args[[ 5 ]] )

# How many cores to use:
options( 'mc.cores'=trailing.args[[ 6 ]] )

# Redis URL and port:
redis.host <- as.character( trailing.args[[ 7 ]] )
redis.port <- as.integer( trailing.args[[ 8 ]] )

# Begin
print( "Starting computation" )
print( paste( "Will be using", options('mc.cores'), "parallel threads." ) )

# Connect to REDIS:
tryCatch(
  redisConnect( host=redis.host, port=redis.port ),
  error=function( e ) {
    err.msg <- paste( "Could not connect to REDIS using URL ",
      redis.host, ":", redis.port, "\n", as.character( e )
    )
    stop( err.msg )
  }
)

# Generate pairs of protein accessions to compute distances for:
accessions <- colnames( annos )
dist.pairs <- distanceIndices( batch.no, batch.size, accessions )

# Partial Sequence Distances:
redisMSet(
  partialSequenceDistances( aa.seqs, dist.pairs, lapply.funk=mclapply )
)

# Partial Domain Architecture Distances:
redisMSet(
  partialDomainArchitectureDistances( annos, dom.weights, dist.pairs,
    lapply.funk=mclapply
  )
)

print( "DONE" )
