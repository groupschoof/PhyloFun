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

# Test extract.annotations
print("Testing extract.annotations(...)")
d <- xmlInternalTreeParse( project.file.path( "test", "Q9ZZX1.xml" ) )
annos <- extract.annotations(d, 'InterPro')
# print( d )
# print( annos )
checkTrue(length(annos) > 0)
checkTrue('IPR004860' %in% annos)
annos <- extract.annotations(d, 'GO')
# print( annos )
checkTrue(length(annos) > 0)
checkTrue('GO:0005739' %in% annos)
annos <- extract.annotations(d, 'Pfam')
checkTrue(length(annos) > 0)
checkTrue('PF00115' %in% annos)

# Test retrieve.annotations
anno.list <- retrieve.annotations(uniprotkb.url('Q0KFR8'))
checkTrue(length(anno.list) > 0)
anno.list <- retrieve.annotations(uniprotkb.url('Protein_1'))
checkTrue(length(anno.list) > 0)

# Test retrieve.annotations.parallel.t
print("Testing retrieve.annotations.parallel.t(...)")
test.accessions <- c('Q0KFR8','B2AFZ7','Q1LSI9','Protein_1')
res <- retrieve.annotations.parallel.t(test.accessions)
checkTrue(is.matrix(res))
checkTrue(nrow(res)==4)
checkTrue(ncol(res)==3)
checkEquals(colnames(res),c('GO','InterPro','Pfam'))
res <- retrieve.annotations.parallel.t(test.accessions, annotations=c('GO'))
checkTrue(is.matrix(res))
checkTrue(nrow(res)==4)
checkTrue(ncol(res)==1)
checkEquals(colnames(res),c('GO'))
res <- retrieve.annotations.parallel.t(test.accessions, annotations=c('GO','Pfam'))
checkTrue(is.matrix(res))
checkTrue(nrow(res)==4)
checkTrue(ncol(res)==2)
checkEquals(colnames(res),c('GO','Pfam'))

# Test retrieve.annotations.parallel
print("Testing retrieve.annotations.parallel(...)")
res <- retrieve.annotations.parallel(test.accessions)
checkTrue(is.matrix(res))
checkTrue(nrow(res)==3)
checkTrue(ncol(res)==4)
checkEquals(rownames(res),c('GO','InterPro','Pfam'))
res <- retrieve.annotations.parallel(test.accessions, annotations=c('GO'))
checkTrue(is.matrix(res))
checkTrue(nrow(res)==1)
checkTrue(ncol(res)==4)
checkEquals(rownames(res),c('GO'))
res <- retrieve.annotations.parallel(test.accessions, annotations=c('GO','Pfam'))
checkTrue(is.matrix(res))
checkTrue(nrow(res)==2)
checkTrue(ncol(res)==4)
checkEquals(rownames(res),c('GO','Pfam'))

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

# Test retrieveSequence
print("Testing retrieveSequence(...)")
aa.seq <- retrieveSequence( getURL('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/Q3EBC8/xml') )
checkEquals(aa.seq,
  'MTMDADAMETETTDQVSASPLHFARSYQVEALEKAIKQNTIVFLETGSGKTLIAIMLLRSYAYLFRKPSPCFCVFLVPQVVLVTQQAEALKMHTDLKVGMYWGDMGVDFWDSSTWKQEVDKYEVLVMTPAILLDALRHSFLSLSMIKVLIVDECHHAGGKHPYACIMREFYHKELNSGTSNVPRIFGMTASLVKTKGENLDSYWKKIHELETLMNSKVYTCENESVLAGFVPFSTPSFKYYQHIKIPSPKRASLVEKLERLTIKHRLSLGTLDLNSSTVDSVEKRLLRISSTLTYCLDDLGILLAQKAAQSLSASQNDSFLWGELNMFSVALVKKFCSDASQEFLAEIPQGLNWSVANINGNAEAGLLTLKTVCLIETLLGYSSLENIRCIIFVDRVITAIVLESLLAEILPNCNNWKTKYVAGNNSGLQNQTRKKQNEIVEDFRRGLVNIIVATSILEEGLDVQSCNLVIRFDPASNICSFIQSRGRARMQNSDYLMMVESGDLLTQSRLMKYLSGGKRMREESLDHSLVPCPPLPDDSDEPLFRVESTGATVTLSSSVSLIYHYCSRLPSDEYFKPAPRFDVNKDQGSCTLYLPKSCPVKEVKAEANNKVLKQAVCLKACIQLHKVGALSDHLVPDMVVAETVSQKLEKIQYNTEQPCYFPPELVSQFSAQPETTYHFYLIRMKPNSPRNFHLNDVLLGTRVVLEDDIGNTSFRLEDHRGTIAVTLSYVGAFHLTQEEVLFCRRFQITLFRVLLDHSVENLMEALNGLHLRDGVALDYLLVPSTHSHETSLIDWEVIRSVNLTSHEVLEKHENCSTNGASRILHTKDGLFCTCVVQNALVYTPHNGYVYCTKGVLNNLNGNSLLTKRNSGDQTYIEYYEERHGIQLNFVDEPLLNGRHIFTLHSYLHMAKKKKEKEHDREFVELPPELCHVILSPISVDMIYSYTFIPSVMQRIESLLIAYNLKKSIPKVNIPTIKVLEAITTKKCEDQFHLESLETLGDSFLKYAVCQQLFQHCHTHHEGLLSTKKDGMISNVMLCQFGCQQKLQGFIRDECFEPKGWMVPGQSSAAYSLVNDTLPESRNIYVASRRNLKRKSVADVVESLIGAYLSEGGELAALMFMNWVGIKVDFTTTKIQRDSPIQAEKLVNVGYMESLLNYSFEDKSLLVEALTHGSYMMPEIPRCYQRLEFLGDSVLDYLITKHLYDKYPCLSPGLLTDMRSASVNNECYALVAVKANLHKHILYASHHLHKHISRTVSEFEQSSLQSTFGWESDISFPKVLGDVIESLAGAIFVDSGYNKEVVFASIKPLLGCMITPETVKLHPVRELTELCQKWQFELSKAKDFDSFTVEVKAKEMSFAHTAKASDKKMAKKLAYKEVLNLLKNSLDY')
# checkEquals( class(retrieveSequence('')), 'try-error' )
# checkEquals( class(retrieveSequence(NA)), 'try-error' )
# checkEquals( class(retrieveSequence(NULL)), 'try-error' )
checkTrue( is.null( retrieveSequence( "Error 12. No entries found", return.error=F ) ) )
aa.doc <- xmlInternalTreeParse( project.file.path( "test", "Q9ZZX1.xml" ) )
checkEquals( retrieveSequence( aa.doc ),
  'MVQRWLYSTNAKDIAVLYFMLAIFSGMAGTAMSLIIRLELAAPGSQYLHGNSQLFNVLVVGHAVLMIFFLVMPALIGGFGNYLLPLMIGATDTAFPRINNIAFWVLPMGLVCLVTSTLVESGAGTGWTVYPPLSSIQAHSGPSVDLAIFALHLTSISSLLGAINFIVTTLNMRTNGMTMHKLPLFVWSIFITAFLLLLSLPVLSAGITMLLLDRNFNTSFFEVSGGGDPILYEHLFWFFGHPEVYILIIPGFGIISHVVSTYSKKPVFGEISMVYAMASIGLLGFLVWSHHMYIVGLDADTRAYFTSATMIIAIPTGIKIFSWLMNPFSKDKNKNKNKKLIRNYQKMNNNNMMKTYLNNNNMIMMNMYKGNLYDIYPRSNRNYIQPNNINKELVVYGYNLESCVGMPTYTNIVKHMVGIPNNILYIMTGILLTDGWIDYTSKKDLDKKTIMEINCRFRLKQSMIHSEYLMYVFMLLSHYCMSYPKMKIAKVKGKSYNQLEFYTRSLPCFTILRYMFYNGRVKIVPNNLYDLLNYESLAHMIMCDGSFVKGGGLYLNLQSFTTKELIFIMNILKIKFNLNCTLHKSRNKYTIYMRVESVKRLFPMIYKYILPSMRYKFDIMLWQKKYNMIN'
)

# Test retrieveSequences
print("Testing retrieveSequences(...)")
test.uris <- c( 
  lapply( c('Q0KFR8','B2AFZ7','Q1LSI9'), uniprotkb.url ),
  "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/NonExistingAccession/xml"
)
uniprot.docs <- getURL( test.uris )
seqs <- sapply( uniprot.docs, function( doc ) {
      try( retrieveSequence( doc ), silent=T )
    })
err.uris <- names( seqs[ grepl("^Error", seqs[], perl=T) ] )
checkEquals( length( err.uris ), 1 )

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

