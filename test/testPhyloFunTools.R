library(RUnit)
library(tools)
library(stringr)
library(Biostrings)
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

# Test parsePhmmerTable
print("Testing parsePhmmerTable(...)")
jr <- scan( file=project.file.path( 'test', 'jackhmmer_out_10_prots.tbl' ),
  what=character(), sep="\n" )
jack.mat <- parsePhmmerTable( jr )
# print( jack.mat[1,] )
checkTrue( ! is.null( jack.mat ) )
checkEquals( class(jack.mat) , 'matrix' )
checkEquals( nrow(jack.mat), length(jr) - 3 )
checkEquals( colnames(jack.mat), c( "hit.name", "query.name", "bit.score" ) )
checkEquals( jack.mat[ 1, ],
  list( "hit.name"="sp|B5YXA4|DNAA_ECO5E", "query.name"="Protein_1", "bit.score"="697.8" )
)

# Test parseBlastTable
print("Testing parseBlastTable(...)")
res.parseBlastTable <- parseBlastTable( read.table(
  project.file.path( 'test', 'test_blast_results.tbl' ) )
)
exp.parseBlastTable <- read.table( text=
'"hit.name" "query.name" "bit.score"
"1" "Query_D3Z7J9" "P63039" 289
"2" "Query_D3Z7J9" "P63038" 289
"3" "Query_D3Z7J9" "D3Z2F2" 288
"4" "Query_D3Z7J9" "P10809" 286
"5" "Query_D3Z7J9" "E7EXB4" 285
"6" "Query_D3Z7J9" "E7ESH4" 285
"7" "Query_D3Z7J9" "B7Z712" 285
"8" "Query_D3Z7J9" "Q5ZL72" 283
"9" "Query_D3Z7J9" "Q803B0" 279
"10" "Query_D3Z7J9" "C9JL25" 265' )
# print( res.parseBlastTable )
# print( exp.parseBlastTable )
checkEquals( colnames( res.parseBlastTable ), colnames( exp.parseBlastTable ) )
checkEquals( rownames( res.parseBlastTable ), rownames( exp.parseBlastTable ) )
checkTrue( all( as.logical( 
  lapply( colnames( exp.parseBlastTable ), function( cn ) {
    checkEquals( res.parseBlastTable[ , cn ], exp.parseBlastTable[ , cn ] )
  } )
) ) )

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

# Test uniqueHomologs
print("Testing uniqueHomologs(...)")
uniqueHomologs( project.file.path( 'test', 'non_unique_hmlgs.fasta' ),
  'tmp.fasta', FALSE )
res.uniqueHomologs <- readAAStringSet( 'tmp.fasta' )
checkTrue( ! is.null( res.uniqueHomologs ) )
checkEquals( length( res.uniqueHomologs ), 167 )
exp.unique.hmlgs.names <- c( "Q8K4K1", "Q9S744", "Q9FYK2", "Q8VHC5", "Q9SRE6",
                            "O35648", "Q545L8", "Q9C8Y1", "O88751", "B1AKR1",
                            "O15182", "P08590", "P02600", "Q9NZU7", "Q8TD86",
                            "P05976", "Q9NPB3", "P57796", "Q9ZSA2", "P12829",
                            "Q9SRR7", "Q38868", "Q9NP86", "O23184", "P25070",
                            "Q9SU00", "Q3E9C0", "Q8RWL2", "A2AVY6", "Q9CQ19",
                            "Q9FMP5", "Q06850", "Q3THE2", "P18666", "Q8VZ50",
                            "Q9ZR02", "Q94AZ4", "Q9S9V0", "Q38870", "Q38872",
                            "Q39016", "Q9ZSA4", "O74435", "P13832", "Q6ZWQ9",
                            "Q38871", "Q38869", "Q42479", "Q9LE22", "Q5XJC3",
                            "Q9ZQE6", "P06704", "Q9ZSA3", "Q9M101", "E5RJF8",
                            "Q38873", "E5RK82", "Q9SSF8", "Q1PFH8", "Q9C9U8",
                            "Q6NLQ6", "Q99MJ8", "Q99MJ7", "Q09510", "Q42438",
                            "Q9FI19", "P93759", "P10916", "Q9XVI9", "P08733",
                            "P51667", "Q9FIH9", "P13833", "P48593", "P63100",
                            "Q63810", "P40423", "Q9M9V8", "Q95XF6", "O14008",
                            "Q8W4I7", "E9QNY3", "Q62082", "Q08331", "Q9SVG9",
                            "Q8CCS7", "Q52K82", "P04466", "P97457", "P22676",
                            "Q9SS31", "Q24214", "P53141", "P48451", "Q801M3",
                            "Q9UU93", "P28470", "Q5SVI8", "Q9QVP4", "Q5NCJ7",
                            "Q86V35", "Q91ZM8", "Q01449", "P19626", "P19625",
                            "Q5AK12", "P53014", "Q9BUA6", "Q55G87", "Q7XJR9",
                            "Q09196", "Q9FKW4", "P25296", "Q9BXU9", "P30188",
                            "E9Q8Y0", "B1AUQ7", "A1BN54", "Q542R1", "Q9JJG7",
                            "P05094", "F8WHE1", "P61601", "P61602", "Q5PQN0",
                            "Q6ZM98", "Q8CBJ9", "Q3B7N2", "P12814", "D3YUI7",
                            "Q7TPR4", "Q9Z1P2", "P09402", "Q9SRE7", "P35609",
                            "P42324", "D3ZCV0", "Q9JI91", "Q5FW75", "Q9UUG5",
                            "P42325", "O88990", "Q08043", "P12815", "P20111",
                            "P04354", "O75340", "P62748", "Q9QXQ0", "P57780",
                            "Q3ULT2", "O43707", "A8Y589", "E9PV73", "P84074",
                            "P84075", "P84076", "Q9FDX6", "Q54X77", "A6NER6",
                            "P18432", "P07171", "A2AS59", "P12658", "P05937",
                            "P06742", "Q8UUX9" )

checkEquals( names( res.uniqueHomologs ), exp.unique.hmlgs.names )
# clean up:
unlink( 'tmp.fasta' )

# Test msaEqual
print("Testing msaEqual(...)")
msa.a <- readAAStringSet( project.file.path( 'test', 'test_msa_equal_A.fasta' ) )
msa.b <- readAAStringSet( project.file.path( 'test', 'test_msa_equal_B.fasta' ) )
msa.c <- readAAStringSet( project.file.path( 'test', 'test_msa_equal_C.fasta' ) )
checkTrue( msaEqual( msa.a, msa.b ) )
checkTrue( ! msaEqual( msa.a, msa.c ) )
checkTrue( ! msaEqual( msa.b, msa.c ) )

# Test bestHits
print("Testing bestHits(...)")
seq.srch.rslt <- matrix( c( 'Query_A', 'Query_A', 'Query_A', 'Query_B',
    'Hit_1', 'Hit_2', 'Hit_3', 'Hit_4', 4, 3, 2, 4 ),
  ncol=3, nrow=4,
  dimnames=list( c(), c( 'query.name', 'hit.name', 'bit.score' ) )
)
res.bestHits <- bestHits( seq.srch.rslt, 'Query_A', 2 )
exp.bestHits <- matrix( c( 'Query_A', 'Query_A', 'Hit_1', 'Hit_2', 4, 3 ),
  nrow=2, ncol=3, dimnames=list( c(),
    c( 'query.name', 'hit.name', 'bit.score' ) )
)
checkEquals( res.bestHits, exp.bestHits ) 
res.bestHits <- bestHits( seq.srch.rslt, 'Query_A', 4 )
exp.bestHits <- matrix( c( 'Query_A', 'Query_A', 'Query_A', 'Hit_1', 'Hit_2', 'Hit_3',  4, 3, 2 ),
  nrow=3, ncol=3, dimnames=list( c(),
    c( 'query.name', 'hit.name', 'bit.score' ) )
)
checkEquals( res.bestHits, exp.bestHits ) 

# Test filterMultipleSequenceAlignment
print("Testing filterMultipleSequenceAlignment(...)")
msa <- readAAStringSet( project.file.path( 'test', 'test_msa.fasta' ) )
res.filterMultipleSequenceAlignment <- filterMultipleSequenceAlignment( msa )
exp.filterMultipleSequenceAlignment <- msa[ 1 ]
checkTrue( msaEqual( res.filterMultipleSequenceAlignment, exp.filterMultipleSequenceAlignment ) ) 
msa.empty <- readAAStringSet( project.file.path( 'test', 'test_msa_empty.fasta' ) )
res.filterMultipleSequenceAlignment <- filterMultipleSequenceAlignment( msa.empty )
exp.filterMultipleSequenceAlignment <- NULL
checkEquals( res.filterMultipleSequenceAlignment, exp.filterMultipleSequenceAlignment ) 

# Test chooseFilteredAlignment
print("Testing chooseFilteredAlignment(...)")
res.chooseAlignment <- chooseFilteredAlignment( msa )
exp.chooseAlignment <- FALSE
# expected to choose unfiltered MSA
checkTrue( ! exp.chooseAlignment )
res.chooseAlignment <- chooseFilteredAlignment( msa.empty )
exp.chooseAlignment <- FALSE
# expected to choose unfiltered MSA
checkTrue( ! exp.chooseAlignment )
msa.good <- readAAStringSet( project.file.path( 'test', 'test_msa_good.fasta' ) )
res.chooseAlignment <- chooseFilteredAlignment( msa.good )
exp.chooseAlignment <- TRUE
# expected to choose FILTERED MSA
checkTrue( exp.chooseAlignment )
