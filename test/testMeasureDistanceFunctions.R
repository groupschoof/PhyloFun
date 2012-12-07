library(RUnit)
library(tools)
library(RCurl)
library(Biostrings)
library(phangorn)

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
src.project.file( 'src','measureDistanceFunctions.R' )
src.project.file( 'src','domainArchitectureSimilarity.R' )
src.project.file( 'src','domainArchitectureSimilarityRedis.R' )
src.project.file( 'src','loadUniprotKBEntries.R' )

# Test measureDistances
print("Testing measureDistances(...)")
f <- file( project.file.path("test", "test_annotations.tbl"), "r" )
annos <- unserialize(f)
close(f)
aa.seqs <- sapply(
  read.AAStringSet( project.file.path("test", "measureDistanceFunctions.fasta") ),
  function(s) toString(s)
)
blast.rslt.tbl <- matrix(
  c( 
    "A0RLX8", "A0PKB2",
    "A0RLX8", "A0K2M8",
    "A0K2M8", "A0KR35",
    "A0Q3U6", "A0KEC3"
  ),
  ncol=2, byrow=T
)
domain.weights.table <- read.table( project.file.path( "test",
    "domain_weights_database.tbl" )
)
go.term <- "GO:0017111"
dists <- measureDistances( go.term, annos, blast.rslt.tbl, aa.seqs,
  domain.weights.table )
print( dists )
checkTrue( ! is.null( dists ) )
checkEquals( class(dists), 'matrix' )
checkEquals( nrow(dists), 4 )
exp.dists <- matrix(
  c(
    1.13, 0.56, T,
    1.13, 0.00, T,
    0.67, 0.32, 0,
    0.79, 1.00, 0
  ),
  byrow=T, nrow=4,
  dimnames=list(
    c( "A0PKB2_A0RLX8", "A0K2M8_A0RLX8", "A0K2M8_A0KR35", "A0KEC3_A0Q3U6" ),
    c( "Sequence.Distance", "Domain.Architecture.Distance", "Share.GO:0017111")
  )
)
checkEquals( dists, exp.dists )

# Test mutationProbabilityDistribution
print("Testing mutationProbabilityDistribution(...)")
dists.test <- matrix( c(0.1, 0.3, 0.6, 1.0, 0.1, 0.3, 0.6, 1.0, 1, 1, 0, 1),
  nrow=4, ncol=3, dimnames=list(
    c("A_B", "A_C", "B_C", "C_D"),
    c("Sequence.Distance", "Domain.Architecture.Distance", "Share.GO:7272727")
  )
)
p.mut.seq <- mutationProbabilityDistribution( dists.test, "Sequence.Distance", "Domain.Architecture.Distance" )
f <- file( project.file.path( "test", "p_mutation_sequence_distance_R_serialized.txt" ), "r" )
exp.p.mut.seq <- unserialize( f )
close( f )
checkEquals( p.mut.seq, exp.p.mut.seq )

p.mut.das <- mutationProbabilityDistribution( dists.test, "Domain.Architecture.Distance", "Sequence.Distance" )
f <- file( project.file.path( "test", "p_mutation_das_distance_R_serialized.txt" ), "r" )
exp.p.mut.das <- unserialize( f )
close( f )
checkEquals( p.mut.das, exp.p.mut.das )
