library(RUnit)
library(tools)
library(stringr)
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
src.project.file('src', 'phyloFunTools.R')

# Test parseJackhmmerTable
print("Testing parseJackhmmerTable(...)")
jr <- scan( file=project.file.path( 'test', 'jackhmmer_out_10_prots.tbl' ),
  what=character(), sep="\n" )
jack.mat <- parseJackhmmerTable( jr )
print( jack.mat[1,] )
checkTrue( ! is.null( jack.mat ) )
checkEquals( class(jack.mat) , 'matrix' )
checkEquals( nrow(jack.mat), length(jr) - 3 )
checkEquals( colnames(jack.mat), c( "hit.name", "query.name" ) )
checkEquals( jack.mat[ 1, ],
  list( "hit.name"="sp|B5YXA4|DNAA_ECO5E", "query.name"="Protein_1" )
)

# Test extractUniprotAccessionFromUniprotName
print("Testing extractUniprotAccessionFromUniprotName(...)")
uni.acc <- extractUniprotAccessionFromUniprotName( "sp|B5YXA4|DNAA_ECO5E" )
checkEquals( uni.acc, "B5YXA4" )
checkTrue( is.na( extractUniprotAccessionFromUniprotName( NA ) ) )

# Test parseInterProScanTable
print("Testing parseInterProScanTable(...)")
ipr <- scan( file=project.file.path( "test", "interproscan_out.tsv" ),
  what=character(), sep="\n" )
ipr.res <- parseInterProScanTable( ipr )
# print( ipr.res )
checkTrue( ! is.null( ipr.res ) )
checkEquals( class(ipr.res), "matrix" )
checkEquals( ncol(ipr.res), 9 )
checkEquals( ipr.res[[ "InterPro", "Protein_5" ]], "IPR011990" )
checkEquals( ipr.res[[ "InterPro", "Protein_4" ]], 
  c("IPR011650", "IPR017439", "IPR002933")
)

# Test commandLineArguments
print("Testing commandLineArguments(...)")
res.commandLineArguments <- commandLineArguments(
  c( '-a', '1', '-b', '2' ), list('b'='0', 'c'='3')
)
exp.commandLineArguments <- list( 'a'='1', 'b'='2', 'c'='3' )
# print( res.commandLineArguments )
checkEquals( res.commandLineArguments, exp.commandLineArguments ) 

# Test sanitizeUniprotAccession
print("Testing sanitizeUniprotAccession(...)")
res.sanitizeUniprotAccession <- sanitizeUniprotAccession( '  sp|MyAccession|MOUSE_SHEEP   ' )
exp.sanitizeUniprotAccession <- 'MyAccession'
checkEquals( res.sanitizeUniprotAccession, exp.sanitizeUniprotAccession ) 
# without pipes
res.sanitizeUniprotAccession <- sanitizeUniprotAccession( '  MyAccession   ' )
exp.sanitizeUniprotAccession <- 'MyAccession'
checkEquals( res.sanitizeUniprotAccession, exp.sanitizeUniprotAccession ) 
# test NULL:
checkTrue( is.null( sanitizeUniprotAccession( NULL ) ) )
