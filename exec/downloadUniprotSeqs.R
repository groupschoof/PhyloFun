require( PhyloFun )

input.args <- commandArgs( trailingOnly = TRUE )

print("Usage: Rscript downloadUniprotSeqs.R path/2/output.fasta [ path/2/accession_per_line.txt ] ")

# Read input
accs <- if ( length( input.args ) > 1 ) {
  scan( input.args[[ 2 ]], what=character(), sep="\n" ) 
} else {
  data( "ukb_proteins_with_trusted_go_annos", package = "PhyloFun" )
  UKB_PROTEINS_WITH_TRUSTED_GO_ANNOS
}

# Download:
downloadSequences( accs, input.args[[ 1 ]] )
