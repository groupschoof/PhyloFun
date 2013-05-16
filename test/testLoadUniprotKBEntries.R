library(RUnit)
library(tools)
library(RCurl)

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
src.project.file('src','loadUniprotKBEntries.R')
src.project.file('src','domainArchitectureSimilarity.R')

# Accessions of Proteins involved in the following tests:
test.accessions <- c('Q0KFR8','B2AFZ7','Q1LSI9','Protein_1')

# Test extractName
print("Testing extractName(...)")
res.extractName <- extractName(
  xmlInternalTreeParse( project.file.path( "test", "Q9ZZX1.xml" ) ),
  xpath.prefix='//xmlns:entry/'
)
exp.extractName <- "Q9ZZX1"
checkEquals( res.extractName, exp.extractName ) 
# With the result of getEntries:
res.extractName <- extractName( getEntries(
  readLines( project.file.path( "test", "Q9ZZX1.xml" ) ) )[[ 1 ]] )
exp.extractName <- "Q9ZZX1"
checkEquals( res.extractName, exp.extractName ) 
# NULL argument should return NULL accession:
checkTrue( is.null( extractName( NULL ) ) )

# Test downloadUniprotDocuments
print("Testing downloadUniprotDocuments(...)")
res.downloadUniprotDocuments <- downloadUniprotDocuments( test.accessions )
# print( res.downloadUniprotDocuments )
checkEquals( length( res.downloadUniprotDocuments ), 3 ) 

# Test retrieveAnnotationsBiomart
print("Testing retrieveAnnotationsBiomart(...)")
bm.annos <- retrieveAnnotationsBiomart(test.accessions)
# print( bm.annos )
checkEquals( class(bm.annos), 'matrix' )
checkEquals( nrow(bm.annos), 2 )
checkEquals( ncol(bm.annos), 4 )
# print(class(( bm.annos[ 'InterPro', 'Protein_1' ] ) ))
checkEquals( length( unlist( bm.annos[ 'InterPro', 'Protein_1' ] ) ), 0 )
checkEquals( length( unlist( bm.annos[ 'GO', 'Protein_1' ] ) ), 0 )
# print( bm.annos['Q0KFR8', 'InterPro'] )
checkTrue( 'IPR020591' %in% unlist( bm.annos[ 'InterPro', 'Q0KFR8' ] ) )
# print( bm.annos['Q0KFR8', 'GO'] )
checkTrue( 'GO:0005524' %in% unlist( bm.annos[ 'GO', 'Q0KFR8' ] ) )

# Test uniprotkb.url
print("Testing uniprotkb.url(...)")
checkEquals(uniprotkb.url('Q0KFR8'),
  'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/Q0KFR8/xml')
checkEquals(uniprotkb.url("sp|P34529|DCR1_CAEEL"),
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/sp%7cP34529%7cDCR1_CAEEL/xml")
checkEquals( uniprotkb.url('Q0KFR8', frmt='fasta'),
  'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/Q0KFR8/fasta'
)

# Test extractAnnotations
print("Testing extractAnnotations(...)")
d <- xmlInternalTreeParse( project.file.path( "test", "Q9ZZX1.xml" ) )
annos <- extractAnnotations(d, 'InterPro')
# print( d )
# print( annos )
checkTrue(length(annos) > 0)
checkTrue('IPR004860' %in% annos)
annos <- extractAnnotations(d, 'GO')
# print( annos )
checkTrue(length(annos) > 0)
checkTrue('GO:0005739' %in% annos)
annos <- extractAnnotations(d, 'Pfam')
checkTrue(length(annos) > 0)
checkTrue('PF00115' %in% annos)

# Test retrieveAnnotations
print("Testing retrieveAnnotations(...)")
res.retrieveAnnotations <- retrieveAnnotations( d, xpath.prefix='//' )
exp.retrieveAnnotations <- list(
  'GO'=c( "GO:0016021", "GO:0005739", "GO:6969696", "GO:0004129",
         "GO:0009055", "GO:0004519", "GO:0020037", "GO:0005506",
         "GO:0009060", "GO:0006316", "GO:0006397", "GO:0008380" ),
  'InterPro'=c( "IPR000883", "IPR023615", "IPR023616", "IPR004860" ),
  'Pfam'=c( "PF00115", "PF03161")
)
# print( res.retrieveAnnotations )
# print( exp.retrieveAnnotations )
checkEquals( res.retrieveAnnotations, exp.retrieveAnnotations ) 

# Test retrieveUniprotAnnotations
# Uniprot constantly changes the annotation of its proteins, hence we
# cannot test for correct annotations.
print("Testing retrieveUniprotAnnotations(...)")
res.retrieveUniprotAnnotations <- retrieveUniprotAnnotations( c( 'sp|Q5ZL72|CH60_CHICK', 'P63039' ) )
checkTrue( ! is.null( res.retrieveUniprotAnnotations ) )
checkEquals( class( res.retrieveUniprotAnnotations ), 'matrix' )
checkEquals( nrow( res.retrieveUniprotAnnotations ), 3 )
checkEquals( rownames( res.retrieveUniprotAnnotations ),
  c( 'GO','InterPro','Pfam') )
checkTrue( ncol( res.retrieveUniprotAnnotations ) > 0 )

# Test uniq.annotations
print("Testing uniq.annotations(...)")
fl <- file(project.file.path('test','test_annotations.tbl'),'r')
an.ma <- unserialize(fl)
close(fl)
u.an.ma <- uniq.annotations(an.ma, 'GO')
checkEquals(u.an.ma, c(NA,"GO:0003688","GO:0005524",
    "GO:0005737", "GO:0006270", "GO:0006275","GO:0017111"))
u.an.ma <- uniq.annotations(an.ma, 'GO', T)
checkEquals(u.an.ma, c("GO:0003688","GO:0005524",
    "GO:0005737", "GO:0006270", "GO:0006275","GO:0017111"))

# Test sharedFunction
print("Testing sharedFunction(...)")
f <- file( project.file.path( "test", "test_annotations.tbl" ), "r" )
shr.func.anno.mtrx <- unserialize( f )
close( f )
shrd.funk.res <- sharedFunction( shr.func.anno.mtrx )
# print( shrd.funk.res )
col.row.nms <- c("A0RLX8", "A0LE53", "Protein_1", "A0PKB2", "A0Q3U6", "A0AEI7",
  "A0K2M8", "A0KR35", "A0KEC3", "A0Q3U7", "A0L3I7")
bools <- c( TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
  TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE,
  FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
  TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
  FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE,
  TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE,
  FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE,
  TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE,
  TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
  TRUE, TRUE, TRUE, TRUE )
exp.shrd.funk.res <- matrix( bools, byrow=T, ncol=11, nrow=11,
    dimnames=list( col.row.nms, col.row.nms )
  )
checkEquals( exp.shrd.funk.res, shrd.funk.res )
# Test function when called to only compute partial shared functions:
shrd.funk.accs <- c( "A0AEI7", "A0KEC3" )
shrd.funk.res <- sharedFunction( shr.func.anno.mtrx, accessions=shrd.funk.accs )
# print( shrd.funk.res )
checkEquals( shrd.funk.res, exp.shrd.funk.res[ shrd.funk.accs, ] )

# Test sharedAnnotation
print("Testing sharedAnnotation(...)")
f <- file( project.file.path( "test", "test_annotations.tbl" ), "r" )
annos <- unserialize( f )
close( f )
shrd.annos.1 <- sharedAnnotation( annos, "GO:0017111" )
# print( shrd.annos.1 )
checkEquals( class(shrd.annos.1), "matrix" )
checkEquals( nrow(shrd.annos.1), 3 )
checkEquals( ncol(shrd.annos.1), 6 )
checkEquals( shrd.annos.1[[ "GO", "A0LE53" ]], "GO:0017111" )
checkEquals( shrd.annos.1[[ "GO", "A0PKB2" ]], annos[[ "GO", "A0PKB2" ]] )
shrd.annos.2 <- sharedAnnotation( annos, "GO:0005737" )
# print( shrd.annos.2 )
checkEquals( class(shrd.annos.2), "matrix" )
checkEquals( nrow(shrd.annos.2), 3 )
checkEquals( ncol(shrd.annos.2), 2 )
checkEquals( shrd.annos.2[[ "GO", "A0Q3U7" ]], annos[[ "GO", "A0Q3U7" ]] )

# Test extractExperimentallyVerifiedGoAnnos
print("Testing extractExperimentallyVerifiedGoAnnos(...)")
rslt <- extractExperimentallyVerifiedGoAnnos(
  xmlInternalTreeParse( project.file.path( "test", "A0AEI7.xml" ) )
)
checkTrue( is.null(rslt) )
rslt <- extractExperimentallyVerifiedGoAnnos(
  xmlInternalTreeParse( project.file.path( "test", "Q9ZZX1.xml" ) )
)
# print( rslt )
exp.rslt <- matrix( list(), ncol=1, nrow=1, dimnames=list( 'GO', 'Q9ZZX1' ) )
exp.rslt[[ 1, 1 ]] <- c( "GO:0004519", "GO:0006316" )
checkEquals( exp.rslt, rslt )
# including evidence codes TAS and IC
rslt <- extractExperimentallyVerifiedGoAnnos(
  xmlInternalTreeParse( project.file.path( "test", "Q9ZZX1.xml" ) ),
  evidence.codes=c( 'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC' )
)
exp.rslt <- matrix( list(), ncol=1, nrow=1, dimnames=list( 'GO', 'Q9ZZX1' ) )
exp.rslt[[ 1, 1 ]] <- c( "GO:0005739", "GO:6969696", "GO:0004519", "GO:0006316" )
checkEquals( exp.rslt, rslt )

# Test retrieveExperimentallyVerifiedGOAnnotations
print("Testing retrieveExperimentallyVerifiedGOAnnotations(...)")
exper.go.annos <- retrieveExperimentallyVerifiedGOAnnotations( c( "A0AEI7", "Q9ZZX1" ) )
# print( exper.go.annos )
checkEquals( exper.go.annos[[ 'GO', 'Q9ZZX1' ]], c( "GO:0004519", "GO:0006316" ) )
checkEquals( ncol(exper.go.annos), 1 )
# Test with 550 accessions, triggering recursive execution:
accs.550 <- as.character( read.table( project.file.path( "test", "550_exp_ver_accs.txt" ) )$V1 )
exper.go.annos.550 <- retrieveExperimentallyVerifiedGOAnnotations( accs.550 )
checkTrue( ! is.null( exper.go.annos.550 ) )
checkEquals( "matrix", class( exper.go.annos.550 ) )
checkEquals( 1, nrow( exper.go.annos.550 ) )
# Uniprot constantly changes their annotations, hence make a very conservative
# test:
checkTrue( ncol( exper.go.annos.550 ) > 400 )
checkEquals( ncol( exper.go.annos.550 ),
  length( unique( colnames( exper.go.annos.550 ) ) )
)

# Test extractRefSeqAccession
print("Testing extractRefSeqAccession(...)")
res.extractRefSeqAccession <- extractRefSeqAccession( 'gi|444722022|gb|ELW62727.1|' )
exp.extractRefSeqAccession <- 'ELW62727.1'
checkEquals( res.extractRefSeqAccession, exp.extractRefSeqAccession ) 
