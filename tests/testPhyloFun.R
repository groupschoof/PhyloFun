library(RUnit)
library(tools)
library(ape)
library(RMySQL)
# In R sourcing other files is not trivial, unfortunately.
# WARNING:
# This method ONLY works for project files in depth one sub dirs!
project.file.path <- function(...) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dir <- dirname(file_path_as_absolute(script.name))
  project.dir <- sub(basename(script.dir), '', script.dir)
  normalizePath(file.path(project.dir, ...))
}
src.project.file <- function(...) {
  source(project.file.path(...))
}
# We set-up required libraries in the test case, not in the R file, as path
# problems will be resolved, as soon as this R package is loaded as such.
src.project.file('src', 'loadUniprotKBEntries.R')
src.project.file('src', 'phyloFun.R')
src.project.file('src', 'geneOntologySQL.R')

# Initialize test data:
#######################
# Load list of mutation probability tables for all measured GO terms:
load( project.file.path( "data", "p_mutation_tables_R_image.bin" ) )

# Test tree is midpoint rooted!
phylo.tree <- read.tree(project.file.path('test', 'test_tree.newick'))
fl <- file(project.file.path('test','test_annotations_2.tbl'),"r")
annotation.matrix <- unserialize(fl)
close(fl)

# Test roundBranchLengths
print("Testing roundBranchLengths(...)")
res.roundBranchLengths <- roundBranchLengths( phylo.tree )
exp.roundBranchLengths <- c( 0.01, 0.90, 1.17, 0.13, 0.13, 0.11, 0.62, 0.75,
                            0.16, 0.27, 0.31, 0.16, 0.10, 0.20, 0.59, 0.23,
                            0.07, 0.14, 0.67, 0.57 )
# print( phylo.tree$edge.length )
# print( res.roundBranchLengths$edge.length )
checkEquals( res.roundBranchLengths$edge.length, exp.roundBranchLengths )

# Test cumulativeBranchLengthsToRoot
print("Testing cumulativeBranchLengthsToRoot(...)")
res.cumulativeBranchLengthsToRoot <- cumulativeBranchLengthsToRoot( 2, phylo.tree )
exp.cumulativeBranchLengthsToRoot <- 1.17788
checkEquals( res.cumulativeBranchLengthsToRoot, exp.cumulativeBranchLengthsToRoot ) 
res.cumulativeBranchLengthsToRoot <- cumulativeBranchLengthsToRoot( 7, phylo.tree )
checkEquals( res.cumulativeBranchLengthsToRoot, exp.cumulativeBranchLengthsToRoot ) 

# Test maxDistanceToRoot
print("Testing maxDistanceToRoot(...)")
res.maxDistanceToRoot <- maxDistanceToRoot( phylo.tree )
exp.maxDistanceToRoot <- 1.17788
checkEquals( res.maxDistanceToRoot, exp.maxDistanceToRoot ) 

# Test mostAppropriateAnnotation
print("Testing mostAppropriateAnnotation(...)")
preds <- setNames( c( 0.5, 0.5, 0.5 ),
  c( 'GO:0009267 & GO:0036170 & GO:0036180 & GO:0071216', 'unknown',
    'GO:0009267 & GO:0036170' )
)
res.mostAppropriateAnnotation <- mostAppropriateAnnotation( preds )
exp.mostAppropriateAnnotation <- c( 'GO:0009267', 'GO:0036170', 'GO:0036180', 'GO:0071216' )
checkEquals( res.mostAppropriateAnnotation, exp.mostAppropriateAnnotation ) 

# Some PhyloFun results to be arguments for tesed methods:
f <- file( project.file.path( 'test', 'test_phyloFun_serialized_result.bin' ), 'r' )
phylo.fun.rslt <- unserialize( f )
close( f )

# Test highScoringAnnotations
print("Testing highScoringAnnotations(...)")

res.highScoringAnnotations <- highScoringAnnotations( phylo.fun.rslt,
  '"Query_Q9XWC3"' )
checkEquals( res.highScoringAnnotations$biological_process,
  c( "GO:0009267", "GO:0036170", "GO:0036180", "GO:0071216" ) ) 
checkTrue( is.na( res.highScoringAnnotations$cellular_component ) )
checkTrue( is.na( res.highScoringAnnotations$molecular_function ) )

f <- file( project.file.path( 'test',
  'Query_Q9NTK1_phyloFun_annotations_serialized.txt'), "r" )
q9ntk1.pf.res <- unserialize( f )
close( f )

res.highScoringAnnotations <- highScoringAnnotations( q9ntk1.pf.res,
  '"Query_Q9NTK1"' )
checkEquals( res.highScoringAnnotations$biological_process,
  c( "GO:0045892", "GO:0046677", "GO:0051260", "GO:0001514", "GO:0007420",
    "GO:0007426", "GO:0051960" ) ) 
checkEquals( res.highScoringAnnotations$cellular_component,
  c( "GO:0005739", "GO:0005634", "GO:0005886", "GO:0019005" ) ) 
checkEquals( res.highScoringAnnotations$molecular_function,
  c( "GO:0003730", "GO:0008135" ) ) 

# Test tree whose Baysian network representation has nodes with unreachable
# states:
phylo.tree.unreachbl.stts <- read.tree(
  project.file.path( 'test', 'test_tree_unreachbl_stts.newick')
)

# Test predictionsToCharacterVector
print("Testing predictionsToCharacterVector(...)")
res.predictionsToCharacterVector <- predictionsToCharacterVector( phylo.fun.rslt, '"Query_Q9XWC3"' )
exp.predictionsToCharacterVector <- c( 'GO:0009267', 'GO:0036170',
  'GO:0036180', 'GO:0071216' )
checkEquals( res.predictionsToCharacterVector, exp.predictionsToCharacterVector ) 
res.predictionsToCharacterVector <- predictionsToCharacterVector( phylo.fun.rslt, 'Query_Q9XWC3' )
checkEquals( res.predictionsToCharacterVector, exp.predictionsToCharacterVector ) 
checkEquals( predictionsToCharacterVector( NULL, 'Query_Q9XWC3' ), 'unknown' )
checkEquals( predictionsToCharacterVector( list(), 'Query_Q9XWC3' ), 'unknown' )
phylo.fun.rslt.unknown <- phylo.fun.rslt
phylo.fun.rslt.unknown$biological_process$pred[ '"Query_Q9XWC3"' ][[1]][[ 1, 'unknown' ]] <- 0.6
checkEquals( predictionsToCharacterVector( phylo.fun.rslt.unknown, 'Query_Q9XWC3' ), 'unknown')

# Test goTermPredictionTable
print("Testing goTermPredictionTable(...)")
res.goTermPredictionTable <- goTermPredictionTable( phylo.fun.rslt, 'Query_Q9XWC3' )
checkEquals( class( res.goTermPredictionTable ), 'data.frame' ) 
checkEquals( ncol( res.goTermPredictionTable ), 8 ) 
checkTrue( nrow( res.goTermPredictionTable ) > 0 )
res.goTermPredictionTable <- goTermPredictionTable( phylo.fun.rslt.unknown, 'Query_Q9XWC3' )
checkEquals( res.goTermPredictionTable,
  matrix(
    c( NA, 'unknown', NA, NA, NA, NA, NA, NA ), byrow=T, ncol=8,
    dimnames=list( c(), c( "id", "name", "term_type", "acc", "is_obsolete",
      "is_root", "is_relation", "relation_distance" ) )
  )
)

# Initialize a database connection to the Gene Ontology
go.con <- connectToGeneOntology()

# Test getDescendantNodes
print("Testing getDescendantNodes(...)")
res.getDescendantNodes <- getDescendantNodes( phylo.tree, 12 )
exp.getDescendantNodes <- c( 13, 14 )
checkEquals( res.getDescendantNodes, exp.getDescendantNodes ) 
checkTrue( is.null( getDescendantNodes( phylo.tree, 8 ) ) )

# Test annotationMatrixForBayesNetwork
print("Testing annotationMatrixForBayesNetwork(...)")
res.annotationMatrixForBayesNetwork <- annotationMatrixForBayesNetwork( annotation.matrix )
# print( res.annotationMatrixForBayesNetwork )
checkTrue( ! is.null( res.annotationMatrixForBayesNetwork ) )
checkEquals( class( res.annotationMatrixForBayesNetwork ), 'matrix' )
checkEquals( ncol( res.annotationMatrixForBayesNetwork ),
  ncol( annotation.matrix ) )
checkEquals( rownames( res.annotationMatrixForBayesNetwork ), 'GO' )
checkEquals( colnames( res.annotationMatrixForBayesNetwork ),
  surroundEachWithQuotes( colnames( annotation.matrix ) ) )
checkEquals( res.annotationMatrixForBayesNetwork[[ 'GO', '"A0RLX8"' ]],
  "GO:0003688 & GO:0005524 & GO:0005737 & GO:0006270 & GO:0006275 & GO:0017111"
)
checkTrue( is.null( annotationMatrixForBayesNetwork( NULL ) ) )
# Check with homologs missing experimentally verified function annotations:
res.annotationMatrixForBayesNetwork <- annotationMatrixForBayesNetwork(
  annotation.matrix, all.accessions=c( colnames( annotation.matrix ),
    'Protein_A', 'Protein_B' )
)
# print( res.annotationMatrixForBayesNetwork )
checkTrue( ! is.null( res.annotationMatrixForBayesNetwork ) )
checkEquals( ncol( res.annotationMatrixForBayesNetwork ),
  ncol( annotation.matrix ) + 2 )
checkEquals( res.annotationMatrixForBayesNetwork[[ 'GO', '"Protein_A"' ]],
  'unknown' )
checkEquals( res.annotationMatrixForBayesNetwork[[ 'GO', '"Protein_B"' ]],
  'unknown' )

# Test annotationToString
print("Testing annotationToString(...)")
res.annotationToString <- annotationToString(  c( "GO_A", "GO_B", "GO_C" ) )
exp.annotationToString <- "GO_A & GO_B & GO_C"
checkEquals( res.annotationToString, exp.annotationToString ) 

# Test goTypeAnnotationMatrices
print("Testing goTypeAnnotationMatrices(...)")
go.type.annos.no.restriction <- goTypeAnnotationMatrices( annotation.matrix, NULL, go.con=go.con )
# print( go.type.annos.no.restriction )
checkEquals( names( go.type.annos.no.restriction ), c( 'biological_process', 'cellular_component', 'molecular_function' ) )
checkEquals( go.type.annos.no.restriction$biological_process[[ 'GO', 'A0K2M8' ]],
  'GO:0006275' )
checkEquals( go.type.annos.no.restriction$molecular_function[[ 'GO', 'A0K2M8' ]],
  'GO:0003688' )

# Test goAnnotationSpaceList
print("Testing goAnnotationSpaceList(...)")
res.annotationSpace <- goAnnotationSpaceList( go.type.annos.no.restriction, unknown.annot=NULL )
exp.annotationSpace <- list(
  'biological_process'=list( c("GO:0006270","GO:0006275"), "GO:0006275", "GO:0006270" ),
  'cellular_component'=list( "GO:0005737" ),
  'molecular_function'=list( c("GO:0003688","GO:0005524","GO:0017111"), "GO:0005524", "GO:0003688" )
)
# print( res.annotationSpace )
checkEquals( res.annotationSpace, exp.annotationSpace ) 
# With UNKOWN annotation
res.annotationSpace <- goAnnotationSpaceList( go.type.annos.no.restriction, unknown.annot='unknown' )
exp.annotationSpace <- list(
  'biological_process'=list( c("GO:0006270","GO:0006275"), "GO:0006275", "GO:0006270", "unknown" ),
  'cellular_component'=list( "GO:0005737", "unknown" ),
  'molecular_function'=list( c("GO:0003688","GO:0005524","GO:0017111"), "GO:0005524", "GO:0003688", "unknown" )
)
# print( res.annotationSpace )
checkEquals( res.annotationSpace, exp.annotationSpace ) 

# Test findMatchingRow
print("Testing findMatchingRow(...)")
p.mut.tbl <- as.matrix( read.table( header=T, text=
"   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
1                          0.00                  0.00                            0.00                             0.00                  0.00                                0                             0.00
2                          0.50                  0.03                            0.00                             0.03                  0.03                                0                             0.03
3                          0.67                  0.06                            0.00                             0.06                  0.06                                0                             0.06
4                          0.75                  0.07                            0.00                             0.07                  0.07                                0                             0.07
5                          0.83                  0.11                            0.00                             0.11                  0.11                                0                             0.11
6                          0.86                  0.12                            0.00                             0.12                  1.50                                1                             1.98
7                          0.87                  1.51                            0.00                             1.51                  1.52                                1                             1.82
8                          0.88                  1.53                            0.00                             1.53                  1.54                                1                             1.84
9                          0.89                  1.55                            0.00                             1.56                  1.58                                1                             1.87
10                         0.90                  1.59                            0.16                             1.61                  1.60                                1                             1.89
11                         0.91                  1.61                            0.00                             1.62                  1.66                                1                             2.29
12                         0.92                  1.67                            0.06                             1.71                  1.92                                1                             2.14" )
)
mtch.col.1 <- findMatchingRow( p.mut.tbl, 0.8, 5 )
# print( mtch.col.1 )
exp.mtch.col.1 <- as.matrix( read.table( header=T, text="
   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
6                          0.86                  0.12                            0.00                             0.12                  1.50                                1                             1.98" )
)
checkEquals( exp.mtch.col.1, mtch.col.1 )

mtch.col.2 <- findMatchingRow( p.mut.tbl, 1.54, 5 )
# print( mtch.col.2 )
exp.mtch.col.2 <- as.matrix( read.table( header=T, text="
   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
8                          0.88                  1.53                            0.00                             1.53                  1.54                                1                             1.84" )
)
checkEquals( exp.mtch.col.2, mtch.col.2 )

mtch.col.3 <- findMatchingRow( p.mut.tbl, 7272, 5 )
# print( mtch.col.3 )
exp.mtch.col.3 <- as.matrix( read.table( header=T, text="
   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
12                         0.92                  1.67                            0.06                             1.71                  1.92                                1                             2.14" )
)
checkEquals( exp.mtch.col.3, mtch.col.3 )

mtch.col.4 <- findMatchingRow(
  matrix( c(0.33, 0.66, 1.0, 0.5, 1.0, 1.5), ncol=2 ), 0.9, 2
)
checkEquals( matrix( c(0.66, 1.0), ncol=2 ), mtch.col.4 )

# Test eliminateUnreachableStates
print("Testing eliminateUnreachableStates(...)")
cpt.states <- c( 'GO_1', 'GO_2', 'unknown' )
cpt <- matrix( c( 0.5, 0.5, 0, 1, 0, 0, 0.5, 0.5, 0 ),
  nrow=3, ncol=3, byrow=T,
  dimnames=list( cpt.states, cpt.states )
)
res.eliminateUnreachableStates <- eliminateUnreachableStates( cpt )
exp.eliminateUnreachableStates <- matrix( c( 0.5, 0.5, 1, 0, 0.5, 0.5 ),
  nrow=3, ncol=2, byrow=T,
  dimnames=list( cpt.states, c( 'GO_1', 'GO_2' ) )
)
checkEquals( res.eliminateUnreachableStates, exp.eliminateUnreachableStates ) 
# No state is unreachable -> check equality:
cpt <- matrix( c( 0.5, 0.3, 0.2, 1, 0, 0, 0.5, 0.5, 0 ),
  nrow=3, ncol=3, byrow=T,
  dimnames=list( cpt.states, cpt.states )
)
res.eliminateUnreachableStates <- eliminateUnreachableStates( cpt )
checkEquals( res.eliminateUnreachableStates, cpt ) 

# Test conditionalProbabilityTable
print("Testing conditionalProbabilityTable(...)")
ua <- list( c( "GO_1", "GO_2", "GO_3" ), c( "GO_1", "GO_2" ), "GO_3" )
p.mut.tbl.lst <- list()
p.mut.tbl.lst[[ "GO_1" ]] <- matrix( c(0.33, 0.66, 1.0, 0.5, 1.0, 1.5), ncol=2 )
p.mut.tbl.lst[[ "GO_2" ]] <- matrix( c(0.25, 0.5, 0.75, 0.5, 1.0, 1.5), ncol=2 )
p.mut.tbl.lst[[ "GO_3" ]] <- matrix( c(0.45, 0.75, 0.98, 0.5, 1.0, 1.5), ncol=2 )
# print( p.mut.tbl.lst )
con.prbs.tbl <- conditionalProbabilityTable( 0.9, c( ua, 'unknown' ), p.mut.tbl.lst, 2 )
# print( con.prbs.tbl )
checkEquals( 1.0, sum( con.prbs.tbl[ , 1 ] ) )
# print( 1 - p.mut.tbl.lst[[ 1 ]][[ 2, 1 ]] )
checkEquals( 1 - p.mut.tbl.lst[[ 3 ]][[ 2, 1 ]], con.prbs.tbl[[ 1, 1 ]] ) 
checkEquals( 1.0, sum( con.prbs.tbl[ , 1 ] ) ) 
checkEquals( 1 - p.mut.tbl.lst[[ 1 ]][[ 2, 1 ]], con.prbs.tbl[[ 2, 2 ]] ) 
checkEquals( 1.0, sum( con.prbs.tbl[ , 3 ] ) ) 
checkEquals( 1 - p.mut.tbl.lst[[ 3 ]][[ 2, 1 ]], con.prbs.tbl[[ 3, 3 ]] ) 
checkEquals( 0, con.prbs.tbl[[ 'unknown', 'unknown' ]] )
checkEquals( 1.0, sum( con.prbs.tbl[ , 'unknown' ] ) )

# Test conditionalProbabilityTables
print("Testing conditionalProbabilityTables(...)")
phylo.tree.4.brnch.lngths <- read.tree( project.file.path( 'test',
  'test_tree_4_branch_lengths.newick' ) )
res.conditionalProbabilityTables <- conditionalProbabilityTables( phylo.tree.4.brnch.lngths,
  c( ua, 'unknown' ), p.mut.tbl.lst, mut.tbl.length.col.indx=2, lapply.funk=lapply )
# print( res.conditionalProbabilityTables )
checkEquals( names( res.conditionalProbabilityTables ),
  as.character( unique( phylo.tree.4.brnch.lngths$edge.length ) ) )
# print(  res.conditionalProbabilityTables[[ '0.59211' ]][ , 1 ])
checkEquals( setNames( c( 0.25, 0.25, 0.25, 0.25 ),
    c( 'GO_1 & GO_2 & GO_3', 'GO_1 & GO_2', 'GO_3', 'unknown' ) ),
  res.conditionalProbabilityTables[[ '0.59211' ]][ , 1 ]
)

# Test mutationProbability
print("Testing mutationProbability(...)")
anno <- c( "A", "B", "C" )
branch.length <- 0.5
mut.prob.tbls <- list( "A"=matrix( c(0.2, 0.5), nrow=1 ),
  "B"=matrix( c(0.3, 0.5), nrow=1 ), "C"=matrix( c(0.4, 0.5), nrow=1 )
)
res.mutationProbability <- mutationProbability( anno, branch.length, mut.prob.tbls, 2 )
exp.mutationProbability <- 0.4
checkEquals( res.mutationProbability, exp.mutationProbability ) 

# Test get.node.label
print("Testing get.node.label(...)")
checkEquals(get.node.label(phylo.tree, 21), "21")
checkTrue(identical(class(get.node.label(phylo.tree, 11)), "character"))

# Test edge.to.formula
print("Testing edge.to.formula(...)")
# Test a tip's formula
indx <- which(phylo.tree$edge[, 2] == 1)
frml <- edge.to.formula(phylo.tree, indx)
checkTrue(identical(class(frml), 'formula'))
checkTrue(length(as.character(frml)) == 2)
checkTrue(grepl('[a-zA-Z]+', as.character(frml)[[2]], perl=T))
indx <- which(phylo.tree$edge[, 1] == (length(phylo.tree$tip.label)+1))[1]
frml <- edge.to.formula(phylo.tree, indx)
checkTrue(! grepl('[a-zA-Z]+', as.character(frml)[[2]], perl=T))

# Test bayesNodes
# Tree in which the state 'unknown' is unreachable
print("Testing bayesNodes(...)")
res.bayesNodes <- bayesNodes( 
  phylo.tree.unreachbl.stts, c( 'GO:0043047', 'unknown' ), lapply.funk=lapply
)
exp.anno.space <- c( 'GO:0043047', 'unknown' )
# print( res.bayesNodes )
checkEquals( length( res.bayesNodes ), 38 ) 
root.cpt <- res.bayesNodes[[ 1 ]]
checkEquals( root.cpt[[ 'values' ]],
   matrix( c( 1, 1 ), ncol=1, dimnames=list( exp.anno.space, c() ) )
)
# print( exp.anno.space )
# print( root.cpt[[ 'levels' ]] )
checkEquals( exp.anno.space, root.cpt[[ 'levels' ]] )
for( i in 2:length(res.bayesNodes) ) {
  desc.cpt <- res.bayesNodes[[ i ]]
  print( 
    checkEquals( desc.cpt[[ 'values' ]],
      matrix( c( 1, 0, 1, 0 ), nrow=2, ncol=2,
        dimnames=list( exp.anno.space, exp.anno.space )
      )
    )
  )
  print( checkEquals( exp.anno.space, desc.cpt[[ 'levels' ]] ) )
}

# Test tree without unreachable 'unkown' state:
go.type.annos <- goTypeAnnotationMatrices( annotation.matrix, go.con=go.con )
anno.space.lst <- goAnnotationSpaceList( go.type.annos )
bys.nds <- bayesNodes( phylo.tree, anno.space.lst$molecular_function,
  lapply.funk=lapply )
checkTrue( length( bys.nds ) == 21 )
root.bys.nd <- bys.nds[[ 1 ]]
# print( root.bys.nd )
# print( uniq.annos )
# print( root.bys.nd$levels )
# print( root.bys.nd$values )
# print( anno.space.lst$molecular_function )
checkTrue( length( root.bys.nd$values ) == length( anno.space.lst$molecular_function ) )
checkEquals( root.bys.nd$levels, as.character( lapply( anno.space.lst$molecular_function, annotationToString ) ) )
# print(bys.nds[[2]]$values)

# Test create bayesian Network
print( "Testing create bayesian Network" )
plist <- compileCPT( bys.nds )
grain.res <- try( grain( plist ), silent=T )
checkTrue( ! identical( 'try-error', class( grain.res ) ) )

# Check conditional probability table for leaf '"A0KEC3"':
print( "Testing conditional probability tables of compiled Bayesian network" )
res.cpt <- as.matrix( plist$'"A0KEC3"' )
exp.states <- c( "GO:0005524 & GO:0017111", "GO:0005524", "unknown" )
checkTrue( ! is.null( res.cpt ) )
checkEquals( class( res.cpt ), c( "parray", "array" ) )
checkEquals( rownames( res.cpt ), exp.states )
checkEquals( colnames( res.cpt ), exp.states )
checkEquals( as.numeric( res.cpt[ , 1 ] ), c( 0.46, 0.27, 0.27 ) )
checkEquals( as.numeric( res.cpt[ , 2 ] ), c( 0.12, 0.76, 0.12 ) )
checkEquals( as.numeric( res.cpt[ , 3 ] ), c( 0.5, 0.5, 0 ) )

# Test getTipsWithNaAnnotation
print("Testing getTipsWithNaAnnotation(...)")
checkEquals( getTipsWithNaAnnotation( phylo.tree, annotation.matrix ),
  c("\"Protein_1\"") )
# Test multiple NA annotations
am <- annotation.matrix
am[1, 6] <- NA
checkEquals(getTipsWithNaAnnotation(phylo.tree, am),
  c("\"Protein_1\"", "\"A0AEI7\""))
checkEquals(
  getTipsWithNaAnnotation(phylo.tree, am, negate=T),
  c( "\"A0RLX8\"", "\"A0LE53\"", "\"A0PKB2\"", "\"A0Q3U6\"", "\"A0K2M8\"",
    "\"A0KR35\"", "\"A0KEC3\"", "\"A0Q3U7\"", "\"A0L3I7\"")
  )

# Clean Up:
dbDisconnect( go.con )
