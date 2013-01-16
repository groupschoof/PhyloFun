library(tools)
library(Biostrings)
library(RCurl)
library(stringr)
library(gRain)
library(RMySQL)
library(XML)
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
# Helper functions:
src.project.file( "src", "phyloFunTools.R" )
src.project.file( "src", "phyloFun.R" )
src.project.file( "src", "loadUniprotKBEntries.R" )
src.project.file( "src", "geneOntologySQL.R" )
load( project.file.path( "data", "p_mutation_tables_R_image.bin" ) )

# Hail User:
print( "Usage: Rscript runPhyloFun.R -q path/2/query_proteins.fasta -j path/2/jackhmmer_results.tbl [ -f path/2/FastTree[MP] ] [ -g path/2/GBlocks ] [ -m path/2/MAFFT ]")
print( '' )
print(
  paste( "WARNING: The PhyloFun pipeline uses other programs to generate multiple sequence alignments (MAFFT),",
    "filter them for conserved regions (GBlocks), and generate a phylogenetic tree of the MSA (FastTree[MP]).",
    "These programs need to be in your $PATH and require protein accessions of your JACKHMMER homolgy searches to be Uniprot accessions.",
    "Finally the accessions of your query proteins should consist only of the following character class [a-zA-Z0-9_-]"
  )
)

# Input
phylo.fun.args <- commandLineArguments( commandArgs(trailingOnly = TRUE), list( 'f'='FastTreeMP', 'g'='Gblocks', 'm'='mafft' ) )

# Read fasta:
aa.seqs <- sapply( read.AAStringSet( phylo.fun.args[[ 'q' ]] ), function(s) toString(s) )
print( paste("Read", length(aa.seqs), "sequences from", phylo.fun.args[[ 'q' ]] ) )

# Parse Jackhmmer results:
jr <- parseJackhmmerTable( 
  scan( file=phylo.fun.args[[ 'j' ]], what=character(), sep="\n" )
)
print( paste( "Parsed JACKHMMER result table. Got", nrow(jr), "query-hitirs" ) )

# Sanitize protein accessions:
accs <- unlist( setNames( lapply( names(aa.seqs), sanitizeUniprotAccession ), names(aa.seqs) ) )
print( "Parsed the query proteins' accessions using function sanitizeUniprotAccession(...) ." )
print( 
  paste( "WARNING: If your query accessions are not matching, PhyloFun will fail to find their",
  "accessions in the JACKHMMER results, nor will Gblocks accept such sequence names!",
  "See function sanitizeUniprotAccession for details." )
)
print( accs )

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

    # Construct the phylogenetic max likelihood tree of the above alignment
    # using FastTree[MP]:
    acc.phyl.tree.file <- gsub( '|', '\\|', paste( acc, '/ml_tree.newick', sep='' ), fixed=T )
    system( paste( 'FastTreeMP <', acc.filtered.msa.file, '>', acc.phyl.tree.file ) ) 

    # Compute probability distributions for GO terms of the three different
    # types 'biological_process', 'cellular_component', and
    # 'molecular_function':
    acc.phyl.tree      <- read.tree( acc.phyl.tree.file )
    acc.hmlgs.annos    <- retrieveExperimentallyVerifiedGOAnnotations( acc.phyl.tree$tip.label )
    acc.go.type.annos  <- goTypeAnnotationMatrices( acc.hmlgs.annos, go.con=go.con )
    acc.go.anno.spaces <- goAnnotationSpaceList( acc.go.type.annos )

    acc.go.predictions <- lapply( c( 'biological_process', 'cellular_component',
      'molecular_function' ), function( go.type ) {
      acc.bayes.evdnc <- annotationMatrixForBayesNetwork( acc.go.type.annos[[ go.type ]] )
      acc.bayes.netw <- grain( compileCPT(
        bayesNodes( acc.phyl.tree, acc.go.type.annos[[ go.type ]], acc.go.anno.spaces[[ go.type ]], lapply.funk=mclapply )
      ) )
      predict.grain( acc.bayes.netw, response=surroundEachWithQuotes( acc ),
        newdata=acc.bayes.evdnc, type='dist'
      )
    })
    
    f <- file( paste( acc, "/phyloFun_R_serialized.txt", sep="" ), "w" )
    serialize( acc.go.predictions, f )
    close( f )
  }
}
