require( PhyloFun )

# Accessions of Proteins involved in the following tests:
test.accessions <- c('Q0KFR8','B2AFZ7','Q1LSI9','Protein_1')

# Test extractName
print("Testing extractName(...)")
res.extractName <- extractName(
  xmlInternalTreeParse( project.file.path(  "Q9ZZX1.xml" ) ),
  xpath.prefix='//xmlns:entry/'
)
exp.extractName <- "Q9ZZX1"
checkEquals( res.extractName, exp.extractName ) 
# With the result of getEntries:
res.extractName <- extractName( getEntries(
  readLines( project.file.path(  "Q9ZZX1.xml" ) ) )[[ 1 ]] )
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
d <- xmlInternalTreeParse( project.file.path(  "Q9ZZX1.xml" ) )
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
fl <- file(project.file.path('test_annotations.tbl'),'r')
an.ma <- unserialize(fl)
close(fl)
u.an.ma <- uniq.annotations(an.ma, 'GO')
checkEquals(u.an.ma, c(NA,"GO:0003688","GO:0005524",
    "GO:0005737", "GO:0006270", "GO:0006275","GO:0017111"))
u.an.ma <- uniq.annotations(an.ma, 'GO', T)
checkEquals(u.an.ma, c("GO:0003688","GO:0005524",
    "GO:0005737", "GO:0006270", "GO:0006275","GO:0017111"))

# Test sharedAnnotation
print("Testing sharedAnnotation(...)")
f <- file( project.file.path(  "test_annotations.tbl" ), "r" )
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

# Test getEvidenceCode
ns <- c( xmlns="http://uniprot.org/uniprot" ) 
db.ref.tag <- '<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
<dbReference type="GO" id="GO:0009815">
<property type="term" value="F:1-aminocyclopropane-1-carboxylate oxidase activity"/>
<property type="evidence" value="IEA:UniProtKB-EC"/>
</dbReference>
</uniprot>'
print("Testing getEvidenceCode(...)")
res.getEvidenceCode <- getEvidenceCode(
  getNodeSet( xmlInternalTreeParse( db.ref.tag ),
    "//xmlns:dbReference", namespace=ns )[[1]] 
)
exp.getEvidenceCode <- 'IEA'
checkEquals( res.getEvidenceCode, exp.getEvidenceCode ) 
db.ref.tag <- '<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
<dbReference type="GO" id="GO:0009815">
<property type="term" value="F:1-aminocyclopropane-1-carboxylate oxidase activity"/>
</dbReference>
</uniprot>'
res.getEvidenceCode <- getEvidenceCode(
  getNodeSet( xmlInternalTreeParse( db.ref.tag ),
    "//xmlns:dbReference", namespace=ns )[[1]] 
)
exp.getEvidenceCode <- NA
checkEquals( res.getEvidenceCode, exp.getEvidenceCode ) 
db.ref.tag <- '<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">
<dbReference type="GO" id="GO:0009815">
<property type="term" value="F:1-aminocyclopropane-1-carboxylate oxidase activity"/>
<property type="evidence" value="IDA:UniProtKB-EC"/>
</dbReference>
</uniprot>'
res.getEvidenceCode <- getEvidenceCode(
  getNodeSet( xmlInternalTreeParse( db.ref.tag ),
    "//xmlns:dbReference", namespace=ns )[[1]] 
)
exp.getEvidenceCode <- 'IDA'
checkEquals( res.getEvidenceCode, exp.getEvidenceCode ) 

# Test extractExperimentallyVerifiedGoAnnos
print("Testing extractExperimentallyVerifiedGoAnnos(...)")
rslt <- extractExperimentallyVerifiedGoAnnos(
  xmlInternalTreeParse( project.file.path(  "A0AEI7.xml" ) )
)
checkTrue( is.null(rslt) )
rslt <- extractExperimentallyVerifiedGoAnnos(
  xmlInternalTreeParse( project.file.path(  "Q9ZZX1.xml" ) ),
  evidence.codes=c( 'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP' )
)
# print( rslt )
exp.rslt <- read.table( stringsAsFactors=FALSE, text="GO:0004519 IMP Q9ZZX1
GO:0006316 IMP Q9ZZX1" )
checkEquals( exp.rslt, rslt )
# including evidence codes TAS and IC
rslt <- extractExperimentallyVerifiedGoAnnos(
  xmlInternalTreeParse( project.file.path(  "Q9ZZX1.xml" ) ),
  evidence.codes=c( 'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC' )
)
exp.rslt <- read.table( stringsAsFactors=FALSE, text="GO:0005739 TAS Q9ZZX1
GO:6969696  IC Q9ZZX1
GO:0004519 IMP Q9ZZX1
GO:0006316 IMP Q9ZZX1" )
checkEquals( exp.rslt, rslt )

# Test retrieveExperimentallyVerifiedGOAnnotations
print("Testing retrieveExperimentallyVerifiedGOAnnotations(...)")
exper.go.annos <- retrieveExperimentallyVerifiedGOAnnotations( c( "A0AEI7", "Q9ZZX1" ) )
# print( exper.go.annos )
checkEquals( length(
    intersect( c( "GO:0004519", "GO:0006316" ), exper.go.annos[, 1] )
  ), 2 )
checkEquals( ncol(exper.go.annos), 3 )
# Test with 550 accessions, triggering recursive execution:
accs.550 <- as.character( read.table( project.file.path(  "550_exp_ver_accs.txt" ) )$V1 )
exper.go.annos.550 <- retrieveExperimentallyVerifiedGOAnnotations( accs.550 )
checkTrue( ! is.null( exper.go.annos.550 ) )
checkEquals( "data.frame", class( exper.go.annos.550 ) )
checkEquals( 3, ncol( exper.go.annos.550 ) )

# Test extractRefSeqAccession
print("Testing extractRefSeqAccession(...)")
res.extractRefSeqAccession <- extractRefSeqAccession( 'gi|444722022|gb|ELW62727.1|' )
exp.extractRefSeqAccession <- 'ELW62727.1'
checkEquals( res.extractRefSeqAccession, exp.extractRefSeqAccession ) 
