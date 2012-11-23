library(tools)
library(Biostrings)
library(parallel)
library(rredis)

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

# Usage:
print( "Usage: Rscript measureDistancesRedis.R path/2/proteins.fasta path/2/accession_per_line.txt redis.host redis.port")
print( "WARNING: Make sure the complete sequence names in the FASTA file are exactly the same as in the accession_per_line file!" )
print( "Remember, that it is required to have setup this computation beforehand. See R script 'initalizeMeasureDistances.R' for details." )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Connect to redis:
redis.host <- trailing.args[[ 3 ]]
redis.port <- trailing.args[[ 4 ]]
# Test, if REDIS is available:
tryCatch(
  redisConnect( host=redis.host, port=redis.port ),
  error=function( e ) {
    err.msg <- paste( "Could not connect to REDIS using URL ",
      redis.host, ":", redis.port, "\n", as.character( e )
    )
    stop( err.msg )
  },
  finaly=function() {
    redisClose()
  }
)

# Init each loop iteration with this function:
init.thread.funk <- function() {
  redisConnect( host=redis.host, port=redis.port )
}

# Close each loop iteration with this function:
close.thread.funk <- function() {
  redisClose()
}

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[1]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[1]]) )

# Read Accessions to process:
accs <- as.character( read.table( trailing.args[[2]] )$V1 )
print( paste("Read", length(accs), "accessions to compute partial distances for. Accessions-File is:", trailing.args[[2]]) )

print( "Starting computation" )

# Partial Sequence Distances:
no.res <- partialSequenceDistancesRedis( aa.seqs, accs )
print( "Computed sequence distances" )

# Partial Domain Architecture Distances:
no.res <- partialDomainArchitectureDistancesRedis( names( aa.seqs ), accs )

print( "DONE" )
