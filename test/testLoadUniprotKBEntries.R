library(RUnit)
library(tools)
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
print(res)
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
