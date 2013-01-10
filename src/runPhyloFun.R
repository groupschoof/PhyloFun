library(tools)
library(Biostrings)
library(RCurl)
library(stringr)
library(gRain)
library(XML)

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
# Helper functions:
src.project.file( "src", "phyloFunTools.R" )

# Hail User:
print( "Usage: Rscript runPhyloFun.R path/2/query_proteins.fasta path/2/jackhmmer_results.tbl [ query-accession-regex='^\\s*(\\S+)\\s' ]")

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( trailing.args[[1]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[1]]) )

# Parse Jackhmmer results:
jr <- parseJackhmmerTable( 
  scan( file=trailing.args[[2]], what=character(), sep="\n" )
)
print( paste( "Parsed JACKHMMER result table. Got", nrow(jr), "query-hit-pairs" ) )

# Parse accessions as String between '>' and the first blank character:
query.acc.regex <- if ( length(trailing.args) == 3 ) trailing.args[[3]] else '^\\s*(\\S+)\\s' 
accs <- unlist( setNames( lapply( names(aa.seqs), function(n) str_match( n, query.acc.regex )[[ 1, 2 ]] ), names(aa.seqs) ) )
print( paste( "Parsed the query proteins' accessions as text matching",
  query.acc.regex, "after the '>' character in file", trailing.args[[ 1 ]],
  "If your query accessions are not matching, PhyloFun will fail to find their accessions in the JACKHMMER results, nor will Gblocks accept such sequence names!" )
)

# Will need DB access to gene ontology:
go.con <- connectToGeneOntology()

# For each query protein, do:
for ( acc in accs ) {
  homologs <- jr[ which( jr[ , 'query.name' ] == acc ), , drop=F ]
  if ( nrow( homologs ) > 0 ) {
    orig.acc <- names( accs[ accs[] == acc ] )
    dir.create( acc )
    hit.accs <- homologs[ , 'hit.name' ]
    hit.uniprot.docs <- downloadUniprotDocuments( hit.accs )
    hit.seqs <- unlist( lapply( hit.uniprot.docs, retrieveSequence, return.error=F ) )
    
    # Generate multiple sequence alignment ( MSA ) using MAFFT:
    acc.hmlgs <- setNames( c( hit.seqs, aa.seqs[ orig.acc ] ), c( names(hit.seqs), acc ) )
    acc.hmlgs.file <- gsub( '|', '\\|', paste( acc, '/homologs.fasta', sep='' ) )
    write.XStringSet( AAStringSet( acc.hmlgs ), file=acc.hmlgs.file )
    acc.msa.file <- gsub( '|', '\\|', paste( acc, "/msa.fasta", sep="" ) )
    system( paste( "mafft --auto", acc.hmlgs.file, ">", acc.msa.file ) )

    # Filter the MSA for highly conserved regions using GBlocks:
    acc.filtered.msa.file <- paste( acc.hmlgs.file, '-gb', sep='' )
    system( paste( 'Gblocks', acc.msa.file, '-b5=h -t=p -p=n' ) )

    # Construct the phylogenetic max likelihoo tree of the above alignment
    # using FastTree[MP]:
    acc.phyl.tree.file <- gsub( '|', '\\|', paste( acc, '/ml_tree.newick', sep='' ), fixed=T )
    system( paste( 'FastTreeMP <', acc.filtered.msa.file, '>', acc.phyl.tree.file ) ) 

    # Compute probability distributions for GO terms of the three different
    # types 'biological_process', 'cellular_component', and
    # 'molecular_function':
    acc.phyl.tree <- read.tree( acc.phyl.tree.file )
    acc.hmlgs.annos <- retrieveExperimentallyVerifiedGOAnnotations( acc.phyl.tree$tip.label )
    acc.go.predictions <- lapply( c( 'biological_process', 'cellular_component',
      'molecular_function' ), function( go.type ) {

    })

  }
}
