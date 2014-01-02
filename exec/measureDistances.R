require( PhyloFun )

# Usage:
print( "Usage: Rscript measureDistances.R path/2/proteins.fasta path/2/blast_out.tbl output_file.RData [round_sequence_distances_to_x_digits=2]")
print( "WARNING: Make sure the _complete_ sequence names in the FASTA file are _exactly the same_ as in the Blast output file!" )

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read fasta:
aa.seqs <- readAAStringSet( trailing.args[[ 1 ]] )
print( paste("Read", length(aa.seqs), "sequences from", trailing.args[[ 1 ]]) )

# Read (preprocessed) Blast 'all vs all' result table:
prot.pairs <- uniqueProteinPairs( read.table( trailing.args[[ 2 ]],
  stringsAsFactors=FALSE, comment.char='', quote='',
  colClasses=c( 'character' ) )
)
print( paste( "Extracted", nrow( prot.pairs ), "UNIQUE protein pairs from",
  trailing.args[[ 2 ]] ) )

# Path to output:
path.2.output <- as.character( trailing.args[[ 3 ]] )
print( paste( "Will write output as binary R table into file", path.2.output ) )

# Round sequence-distances to how many digits:
round.2.digits <- if ( length( trailing.args ) == 4 ) {
  as.numeric( trailing.args[[ 4 ]] )
} else {
  2
}

# Begin
print( "Starting computation" )

prot.pairs <- cbind( as.data.frame( prot.pairs, stringsAsFactors=FALSE ),
  'sequence_distance'=numeric( nrow( prot.pairs ) ) )

for ( i in 1:nrow( prot.pairs ) ) {
  p1 <- toString( aa.seqs[[ prot.pairs[[ i, 1 ]] ]] )
  p2 <- toString( aa.seqs[[ prot.pairs[[ i, 2 ]] ]] )
  prot.pairs[[ i, 3 ]] <- round(
    pairwiseSequenceDistance( p1, p2 ),
    round.2.digits
  )
}

# Save results as binary RData:
save( prot.pairs, file=path.2.output )

print( "DONE" )
