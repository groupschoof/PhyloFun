library(RUnit)
library(tools)
library(stringr)
library( RMySQL )

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
src.project.file( 'src', 'evaluation.R' )
src.project.file( 'src', 'phyloFunTools.R' )
src.project.file( 'src', 'geneOntologySQL.R' )

# Setup:
go.con <- connectToGeneOntology()

# Test precision
print("Testing precision(...)")
true.gos <- c( 'A', 'B', 'C' )
res.precision <- precision( true.gos, true.gos, go.con=go.con )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( true.gos, 'D' ), true.gos, go.con=go.con )
exp.precision <- 0.75
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( 'A', 'B' ), true.gos, go.con=go.con )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( 'D' ), true.gos, go.con=go.con )
exp.precision <- 0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( 'A' ), c(), go.con=go.con )
exp.precision <- 0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( ), true.gos, go.con=go.con )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( ), true.gos, go.con=go.con )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

# pred.gos have just a single false positive:
pred.gos <- c( "GO:0004252", "GO:0005634", "GO:0004175", "GO:0044699",
  "GO:0000122" )
ref.gos <- c( "GO:0004252", "GO:0005634", "GO:0021551" )
res.precision <- precision( pred.gos, ref.gos, go.con=go.con )
exp.precision <- 2 / 3
# print( res.precision )
checkEquals( res.precision, exp.precision ) 

# Test recall
print("Testing recall(...)")
true.gos <- c( 'A', 'B', 'C' )
res.recall <- recall( true.gos, true.gos )
exp.recall <- 1.0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( true.gos, 'D' ), true.gos )
exp.recall <- 1.0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( 'A', 'B' ), true.gos )
exp.recall <- 2 / 3
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( 'D' ), true.gos )
exp.recall <- 0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( 'A' ), c() )
exp.recall <- 1.0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( ), true.gos )
exp.recall <- 0
checkEquals( res.recall, exp.recall ) 

# Test fScore
print("Testing fScore(...)")
res.fScore <- fScore( true.gos, true.gos )
exp.fScore <- 1.0
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c( 'A' ), true.gos )
exp.fScore <- 0.5
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c(), c() )
exp.fScore <- 1.0
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c(), true.gos )
exp.fScore <- 0
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c( 'A' ), c() )
exp.fScore <- 0
checkEquals( res.fScore, exp.fScore ) 

# Test parseBlast2GOresults
print("Testing parseBlast2GOresults(...)")
res.parseBlast2GOresults <- parseBlast2GOresults( readLines( project.file.path( 'test', 'blast2GO.annot' ) ) )
checkTrue( ! is.null( res.parseBlast2GOresults ) )
checkEquals( class( res.parseBlast2GOresults ), 'matrix' )
checkEquals( ncol( res.parseBlast2GOresults ), 97 )
checkEquals( res.parseBlast2GOresults[[ 'GO', 'G4MTK6' ]], 'GO:0043581' )
checkEquals( res.parseBlast2GOresults[[ 'GO', 'Q9USJ5' ]], c( "GO:0006779", "GO:0004853" ) )

# Test parseInterProScan2GOresults
print("Testing parseInterProScan2GOresults(...)")
res.parseInterProScan2GOresults <- parseInterProScan2GOresults( readLines( project.file.path( 'test', 'interproscan_out.tsv' ) ) )
checkTrue( ! is.null( res.parseInterProScan2GOresults ) )
checkEquals( class( res.parseInterProScan2GOresults ), 'matrix' )
checkEquals( ncol( res.parseInterProScan2GOresults ), 7 )
checkEquals( res.parseInterProScan2GOresults[[ 'GO', 'Protein_4' ]], 'GO:0016787' )
checkEquals( res.parseInterProScan2GOresults[[ 'GO', 'Protein_4' ]], 'GO:0016787' )

# Test falsePositiveRate
print("Testing falsePositiveRate(...)")
res.falsePositiveRate <- falsePositiveRate( c( 'A', 'B', 'B', 'C' ), c( 'A', 'B' ) )
checkEquals( 1/3, res.falsePositiveRate ) 
res.falsePositiveRate <- falsePositiveRate( c( 'A', 'B', 'B' ), c( 'A', 'B' ) )
checkEquals( 0.0, res.falsePositiveRate ) 
res.falsePositiveRate <- falsePositiveRate( c( 'A' ), c( 'B' ) )
checkEquals( 1.0, res.falsePositiveRate ) 

##############
# Test rates #
##############
# Predicted annotations
pred.annos <- matrix( list(), nrow=1, ncol=2,
  dimnames=list( 'GO', c( 'Prot_A', 'Prot_B' ) )
)
pred.annos[[ 1, 1 ]] <- c( 'GO_1' )
pred.annos[[ 1, 2 ]] <- c( 'GO_3' )
# Reference annotations
ref.annos <- matrix( list(), nrow=1, ncol=2,
  dimnames=list( 'GO', c( 'Prot_A', 'Prot_B' ) )
)
ref.annos[[ 1, 1 ]] <- c( 'GO_1', 'GO_2' )
ref.annos[[ 1, 2 ]] <- c( 'GO_3', 'GO_4' )
# Protein accessions
prot.accs <- colnames( ref.annos )

# Test fScores
print("Testing fScores(...)")
checkEquals( c( 2/3, 2/3 ),
  as.numeric( fScores( prot.accs, pred.annos, reference.annotations=ref.annos ) )
)

# Test falsePositiveRates
print("Testing falsePositiveRates(...)")
res.falsePositiveRates <- as.numeric( falsePositiveRates( prot.accs, pred.annos, reference.annotations=ref.annos ) )
exp.falsePositiveRates <- c( 0, 0 )
checkEquals( res.falsePositiveRates, exp.falsePositiveRates ) 

# Test recallRates
print("Testing recallRates(...)")
res.recallRates <- as.numeric( recallRates( prot.accs, pred.annos, reference.annotations=ref.annos ) )
exp.recallRates <- c( 0.5, 0.5 )
# print( res.recallRates )
checkEquals( res.recallRates, exp.recallRates ) 

# Test falsePositives
print("Testing falsePositives(...)")
pred.gos <- c( "GO:0004175", "GO:0044699", "GO:0000122" )
ref.gos <- c( "GO:0004252", "GO:0005634", "GO:0021551" )
res.falseNegatives <- falsePositives( pred.gos, ref.gos )
exp.falseNegatives <- "GO:0000122"
# print( res.falseNegatives )
checkEquals( res.falseNegatives, exp.falseNegatives ) 

# Test falsePositivesUpperBound
print("Testing falsePositivesUpperBound(...)")
# Child, Parent and a false positive:
pred.gos <- c( 'GO:0003824', 'GO:0070011', 'GO:0001906' )
ref.go <- 'GO:0016787'
res.falsePositivesUpperBound <- falsePositivesUpperBound( pred.gos, ref.go,
  go.con=go.con )
exp.falsePositivesUpperBound <- c( 'GO:0001906' )
checkEquals( res.falsePositivesUpperBound, exp.falsePositivesUpperBound ) 

# Test truePositivesUpperBound
print("Testing truePositivesUpperBound(...)")
res.truePositivesUpperBound <- truePositivesUpperBound( pred.gos, ref.go,
  go.con=go.con )
exp.truePositivesUpperBound <- c( 'GO:0003824', 'GO:0070011' ) 
checkEquals( res.truePositivesUpperBound, exp.truePositivesUpperBound ) 

# Clean up:
dbDisconnect( go.con )
