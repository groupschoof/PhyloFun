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
print( "Usage: Rscript initializeMeasureDistances.R path/2/domain_weights_database.tbl path/2/domain_annotations_R_serialized.txt redis.host redis.port")

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

# Load domain weights table:
dom.weights <- read.table( trailing.args[[ 1 ]] )
initializeDomainWeights( dom.weights )
print( paste("Initialized", nrow(dom.weights), "domain weights from table", trailing.args[[ 1 ]]) )

# Load Uniprot function (InterPro and GO) Annotations:
f <- file( trailing.args[[ 2 ]], "r" )
annos <- unserialize( f )
close( f )
initializeDomainAnnotations( annos )
print( paste("Initialized", ncol(annos), "function annotations from", trailing.args[[ 2 ]]) )

print( "DONE" )
