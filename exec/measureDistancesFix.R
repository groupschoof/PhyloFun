require( 'PhyloFun' )

# Usage:
print( "Usage: Rscript measureDistancesFix.R path/2/protein_pairs_with_distances.tbl path/2/gene_ontology_annotations.tbl path/2/output [path/2/go_term_sublist.tbl]")

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read protein pairs with their respective distances:
prot.pairs <- read.table( trailing.args[[ 1 ]], stringsAsFactors=FALSE,
  comment.char='', quote='',
  colClasses=c( 'character', 'character', 'numeric' )
)

# Read Gene Ontology (GO) annotations:
go.annos <- read.table( trailing.args[[ 2 ]], stringsAsFactors=FALSE,
  comment.char='', quote='', colClasses=c( 'character' ) )

# Write Mutation Probability Distributions measured for each argument GO term
# into output:
output <- trailing.args[[ 3 ]]

# For which GO terms shall mutation probability distributions be calculated? -
# Default is ALL!
go.terms <- if ( length( trailing.args ) > 3 ) {
  read.table( trailing.args[[ 4 ]], stringsAsFactors=FALSE, comment.char='',
    quote='', colClasses=c( 'character' )
  )[ , 1 ]
} else {
  unique( go.annos[ , 1 ] )
}

# Begin
print( "Starting computation" )

p.mut.dists <- setNames(
  lapply( go.terms, function( go.term ) {
    pMutationMinMaxParentValues(
      mutationProbabilityDistribution( prot.pairs, go.annos, go.term ),
      p.column='p.mutation|Sequence.Distance',
      parent.columns='Sequence.Distance'
    )
  } ),
  go.terms
)
# Write named list of each GO terms' p( mut | seq.dist ) into a binary RData
# file:
pmds.no.null <- p.mut.dists[
  which( ! as.logical( lapply( p.mut.dists, is.null ) ) )
]
save( pmds.no.null, file=output )

# End
print( "DONE" )
