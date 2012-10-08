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

# Test uniprotkb.url
print("Testing uniprotkb.url(...)")
checkEquals(uniprotkb.url('Q0KFR8'),
  'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/Q0KFR8/xml')
checkEquals(uniprotkb.url("sp|P34529|DCR1_CAEEL"),
"http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/sp%7cP34529%7cDCR1_CAEEL/xml")

# Test extract.annotations
print("Testing extract.annotations(...)")
d <- xmlInternalTreeParse(uniprotkb.url('Q0KFR8'))
annos <- extract.annotations(d, 'InterPro')
checkTrue(length(annos) > 0)
checkTrue('IPR003593' %in% annos)
annos <- extract.annotations(d, 'GO')
checkTrue(length(annos) > 0)
checkTrue('GO:0006270' %in% annos)
annos <- extract.annotations(d, 'Pfam')
checkTrue(length(annos) > 0)
checkTrue('PF11638' %in% annos)

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
test.accessions <- c('Q0KFR8','B2AFZ7','Q1LSI9','Protein_1')
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

# Test retrieveSequence
print("Testing retrieveSequence(...)")
aa.seq <- retrieveSequence(xmlInternalTreeParse(getURL('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/Q3EBC8/xml')))
checkEquals(aa.seq,
  'MTMDADAMETETTDQVSASPLHFARSYQVEALEKAIKQNTIVFLETGSGKTLIAIMLLRSYAYLFRKPSPCFCVFLVPQVVLVTQQAEALKMHTDLKVGMYWGDMGVDFWDSSTWKQEVDKYEVLVMTPAILLDALRHSFLSLSMIKVLIVDECHHAGGKHPYACIMREFYHKELNSGTSNVPRIFGMTASLVKTKGENLDSYWKKIHELETLMNSKVYTCENESVLAGFVPFSTPSFKYYQHIKIPSPKRASLVEKLERLTIKHRLSLGTLDLNSSTVDSVEKRLLRISSTLTYCLDDLGILLAQKAAQSLSASQNDSFLWGELNMFSVALVKKFCSDASQEFLAEIPQGLNWSVANINGNAEAGLLTLKTVCLIETLLGYSSLENIRCIIFVDRVITAIVLESLLAEILPNCNNWKTKYVAGNNSGLQNQTRKKQNEIVEDFRRGLVNIIVATSILEEGLDVQSCNLVIRFDPASNICSFIQSRGRARMQNSDYLMMVESGDLLTQSRLMKYLSGGKRMREESLDHSLVPCPPLPDDSDEPLFRVESTGATVTLSSSVSLIYHYCSRLPSDEYFKPAPRFDVNKDQGSCTLYLPKSCPVKEVKAEANNKVLKQAVCLKACIQLHKVGALSDHLVPDMVVAETVSQKLEKIQYNTEQPCYFPPELVSQFSAQPETTYHFYLIRMKPNSPRNFHLNDVLLGTRVVLEDDIGNTSFRLEDHRGTIAVTLSYVGAFHLTQEEVLFCRRFQITLFRVLLDHSVENLMEALNGLHLRDGVALDYLLVPSTHSHETSLIDWEVIRSVNLTSHEVLEKHENCSTNGASRILHTKDGLFCTCVVQNALVYTPHNGYVYCTKGVLNNLNGNSLLTKRNSGDQTYIEYYEERHGIQLNFVDEPLLNGRHIFTLHSYLHMAKKKKEKEHDREFVELPPELCHVILSPISVDMIYSYTFIPSVMQRIESLLIAYNLKKSIPKVNIPTIKVLEAITTKKCEDQFHLESLETLGDSFLKYAVCQQLFQHCHTHHEGLLSTKKDGMISNVMLCQFGCQQKLQGFIRDECFEPKGWMVPGQSSAAYSLVNDTLPESRNIYVASRRNLKRKSVADVVESLIGAYLSEGGELAALMFMNWVGIKVDFTTTKIQRDSPIQAEKLVNVGYMESLLNYSFEDKSLLVEALTHGSYMMPEIPRCYQRLEFLGDSVLDYLITKHLYDKYPCLSPGLLTDMRSASVNNECYALVAVKANLHKHILYASHHLHKHISRTVSEFEQSSLQSTFGWESDISFPKVLGDVIESLAGAIFVDSGYNKEVVFASIKPLLGCMITPETVKLHPVRELTELCQKWQFELSKAKDFDSFTVEVKAKEMSFAHTAKASDKKMAKKLAYKEVLNLLKNSLDY')
checkEquals( class(retrieveSequence('')), 'try-error' )
checkEquals( class(retrieveSequence(NA)), 'try-error' )
checkEquals( class(retrieveSequence(NULL)), 'try-error' )

# Test retrieveAnnotationsBiomart
print("Testing retrieveAnnotationsBiomart(...)")
bm.annos <- retrieveAnnotationsBiomart(test.accessions)
checkEquals( class(bm.annos), 'matrix' )
checkEquals( nrow(bm.annos), 4 )
checkEquals( ncol(bm.annos), 2 )
checkEquals( length( unlist( bm.annos['Protein_1', 'InterPro'] ) ), 0 )
checkEquals( length( unlist( bm.annos['Protein_1', 'GO'] ) ), 0 )
# print( bm.annos['Q0KFR8', 'InterPro'] )
checkTrue( 'IPR020591' %in% unlist( bm.annos['Q0KFR8', 'InterPro'] ) )
# print( bm.annos['Q0KFR8', 'GO'] )
checkTrue( 'GO:0005524' %in% unlist( bm.annos['Q0KFR8', 'GO'] ) )
