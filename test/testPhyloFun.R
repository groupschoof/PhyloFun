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

go.con <- connectToGeneOntology()

# Test annotationToString
print("Testing annotationToString(...)")
res.annotationToString <- annotationToString(  c( "GO_A", "GO_B", "GO_C" ) )
exp.annotationToString <- "GO_A & GO_B & GO_C"
checkEquals( res.annotationToString, exp.annotationToString ) 

# Test goTypeAnnotationMatrices
print("Testing goTypeAnnotationMatrices(...)")
go.type.annos.no.restriction <- goTypeAnnotationMatrices( annotation.matrix, NULL, go.con )
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

# Test findMatchingColumn
print("Testing findMatchingColumn(...)")
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
mtch.col.1 <- findMatchingColumn( p.mut.tbl, 0.8, 5 )
# print( mtch.col.1 )
exp.mtch.col.1 <- as.matrix( read.table( header=T, text="
   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
6                          0.86                  0.12                            0.00                             0.12                  1.50                                1                             1.98" )
)
checkEquals( exp.mtch.col.1, mtch.col.1 )

mtch.col.2 <- findMatchingColumn( p.mut.tbl, 1.54, 5 )
# print( mtch.col.2 )
exp.mtch.col.2 <- as.matrix( read.table( header=T, text="
   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
8                          0.88                  1.53                            0.00                             1.53                  1.54                                1                             1.84" )
)
checkEquals( exp.mtch.col.2, mtch.col.2 )

mtch.col.3 <- findMatchingColumn( p.mut.tbl, 7272, 5 )
# print( mtch.col.3 )
exp.mtch.col.3 <- as.matrix( read.table( header=T, text="
   p.mutation.Sequence.Distance min.Sequence.Distance min.Domain.Architecture.Distance min.Euclidean.Distance.To.Origin max.Sequence.Distance max.Domain.Architecture.Distance max.Euclidean.Distance.To.Origin
12                         0.92                  1.67                            0.06                             1.71                  1.92                                1                             2.14" )
)
checkEquals( exp.mtch.col.3, mtch.col.3 )

mtch.col.4 <- findMatchingColumn(
  matrix( c(0.33, 0.66, 1.0, 0.5, 1.0, 1.5), ncol=2 ), 0.9, 2
)
checkEquals( matrix( c(0.66, 1.0), ncol=2 ), mtch.col.4 )

# Test conditional.probs.tbl
print("Testing conditional.probs.tbl(...)")
ua <- list( c( "GO_1", "GO_2", "GO_3" ), c( "GO_1", "GO_2" ), "GO_3" )
p.mut.tbl.lst <- list()
p.mut.tbl.lst[[ "GO_1" ]] <- matrix( c(0.33, 0.66, 1.0, 0.5, 1.0, 1.5), ncol=2 )
p.mut.tbl.lst[[ "GO_2" ]] <- matrix( c(0.25, 0.5, 0.75, 0.5, 1.0, 1.5), ncol=2 )
p.mut.tbl.lst[[ "GO_3" ]] <- matrix( c(0.45, 0.75, 0.98, 0.5, 1.0, 1.5), ncol=2 )
# print( p.mut.tbl.lst )
con.prbs.tbl <- conditional.probs.tbl( 0.9, c( ua, 'unknown' ), p.mut.tbl.lst, 2 )
# print( con.prbs.tbl )
checkEquals( 1.0, sum( con.prbs.tbl[ 1, ] ) )
# print( 1 - p.mut.tbl.lst[[ 1 ]][[ 2, 1 ]] )
checkEquals( 1 - p.mut.tbl.lst[[ 3 ]][[ 2, 1 ]], con.prbs.tbl[[ 1, 1 ]] ) 
checkEquals( 1.0, sum( con.prbs.tbl[ 1, ] ) ) 
checkEquals( 1 - p.mut.tbl.lst[[ 1 ]][[ 2, 1 ]], con.prbs.tbl[[ 2, 2 ]] ) 
checkEquals( 1.0, sum( con.prbs.tbl[ 3, ] ) ) 
checkEquals( 1 - p.mut.tbl.lst[[ 3 ]][[ 2, 1 ]], con.prbs.tbl[[ 3, 3 ]] ) 
checkEquals( 0, con.prbs.tbl[[ 'unknown', 'unknown' ]] )
checkEquals( 1.0, sum( con.prbs.tbl[ 'unknown', ] ) )

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
print("Testing bayesNodes(...)")
go.type.annos <- goTypeAnnotationMatrices( annotation.matrix, go.con=go.con )
anno.space.lst <- goAnnotationSpaceList( go.type.annos )
bys.nds <- bayesNodes( phylo.tree, go.type.annos$molecular_function, anno.space.lst$molecular_function )
# print( bys.nds )
checkTrue( length( bys.nds ) == nrow( phylo.tree$edge ) + 1 )
root.bys.nd <- bys.nds[[ 1 ]]
# print( root.bys.nd )
# print( uniq.annos )
# print( root.bys.nd$levels )
# print( root.bys.nd$values )
# print( anno.space.lst$molecular_function )
checkTrue( length( root.bys.nd$values ) == length( anno.space.lst$molecular_function ) )
checkEquals( root.bys.nd$levels, lapply( anno.space.lst$molecular_function, annotationToString ) )
# print(bys.nds[[2]]$values)

# Test create bayesian Network
print( "Testing create bayesian Network" )
grain.res <- try( grain( compileCPT( bys.nds ) ), silent=T )
checkTrue( ! identical( 'try-error', class( grain.res ) ) )

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

# Test queryPhylBayesNetwork
print("Testing queryPhylBayesNetwork(...)")
# Test with NO experimentally verified evidence in any tip:
phyl.tree.1 <- read.tree( project.file.path( "test", "test_tree.newick" ) )
annos.1 <- retrieveAnnotationsBiomart( phyl.tree.1$tip.label )
prediction.result <- try(
  queryPhylBayesNetwork( phyl.tree.1, annos.1, "\"Protein_1\"" ),
  silent=F
)
# print( prediction.result )
checkTrue( ! identical( class( prediction.result ), 'try-error' ) )
checkTrue( identical( class( prediction.result ), "list" ) )
# Test with a large tree ( 646 leaves ) with evidence of two Proteins with
# experimentally verified functions:
phyl.tree.2 <- read.tree( project.file.path( "test", "test_tree_large.newick" ) )
annos.2 <- retrieveAnnotationsBiomart( phyl.tree.2$tip.label )
prediction.result <- try(
  queryPhylBayesNetwork( phyl.tree.2, annos.2, "\"Protein_1\"" ),
  silent=F
)
# print( prediction.result )
checkTrue( ! identical( class( prediction.result ), 'try-error' ) )
checkTrue( identical( class( prediction.result ), "list" ) )

# Clean Up:
dbDisconnect( go.con )
