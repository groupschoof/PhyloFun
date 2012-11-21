library(tools)
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
src.project.file('src','domainArchitectureSimilarityRedis.R')

# Usage:
print( "Usage: Rscript loadPartialDistanceMatricesIntoRedis.R path/2/partialSequenceDistanceMatrix.tbl path/2/partialDomainArchitectureDistanceMatrix.tbl redis.host redis.port")

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Connect to redis:
redis.host <- trailing.args[[ 3 ]]
redis.port <- trailing.args[[ 4 ]]
tryCatch(
  redisConnect( host=redis.host, port=redis.port ),
  error=function( e ) {
    err.msg <- paste( "Could not connect to REDIS using URL ",
      redis.host, ":", redis.port, "\n", as.character( e )
    )
    stop( err.msg )
  }
)

# Partial Sequence Distances:
part.seq.dists.file <- trailing.args[[ 1 ]]
if ( file.exists( part.seq.dists.file ) ) {
  seq.dists <- read.table( part.seq.dists.file )
  for ( acc.1 in rownames( seq.dists ) ) {
    lapply( colnames( seq.dists ), function( acc.2 ) {
      k <- pairwiseDistanceKey( acc.1, acc.2 )
      if ( is.null( redisGet(k) ) ) {
        redisSet( k, seq.dists[[ acc.1, acc.2 ]] )
      }
    })
  }
  print( "uploaded sequence distances" )
} else {
  print( paste( "Partial sequence distances file does not exists. Did not upload them.", 
    part.seq.dists.file)
  )
}

# Partial Domain Architecture Distances:
part.das.dists.file <- trailing.args[[ 2 ]]
if ( file.exists( part.das.dists.file ) ) {
  das.dists <- read.table( part.das.dists.file )
  for ( acc.pattern in rownames( das.dists ) ) {
    lapply( colnames( das.dists ), function( acc.subject ) {
      k <- pairwiseDistanceKey( acc.pattern, acc.subject, distance.type="seq_dist" )
      if ( is.null( redisGet(k) ) ) {
        redisSet( k, das.dists[[ acc.pattern, acc.subject ]] )
      }
    })
  }
  print( "uploaded das distances" )
} else {
  print( paste( "Partial domain architecture sequence distances file does not exists. Did not upload them.",
      part.das.dists.file)
  )
}

print( "DONE" )
