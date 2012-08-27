library(XML)
library(parallel)

uniprotkb.url <- function(accession) {
  paste(
    'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/',
    accession,
    '/xml',
    sep=''
    )
}

extract.annotations <- function(doc, type) {
  ns <- c(xmlns="http://uniprot.org/uniprot");
  ndst <- getNodeSet(doc,
    paste("//xmlns:dbReference[@type='",type,"']",sep=""),
    namespaces=ns
    );
  vapply(ndst, xmlGetAttr, vector(mode='character', length=1), 'id');
}

retrieve.annotations <- function(url, annotations=c('GO','InterPro','Pfam')) {
  doc <- try( xmlInternalTreeParse(url), silent=T );
  sapply(annotations,
    function(type) {
      setNames(list(
          if(! identical(class(doc),"try-error")) {
            extract.annotations(doc,type)
          } else {
            NA
          }), type)
    },
    USE.NAMES=F)
}

retrieve.annotations.parallel <- function(accessions, ...) {
  do.call('rbind',
    mclapply(accessions,
    function(acc) {
      acc.annos <- retrieve.annotations(uniprotkb.url(acc),...)
      matrix(acc.annos, nrow=1, dimnames=list(c(acc),names(acc.annos)))
    },
    mc.preschedule=F, mc.cores=50)
  )
}

retrieve.annotations.parallel.t <- function(accessions, ...) {
  do.call('cbind',
    mclapply(accessions,
    function(acc) {
      acc.annos <- retrieve.annotations(uniprotkb.url(acc),...)
      matrix(acc.annos, nrow=length(acc.annos),
        dimnames=list(names(acc.annos),c(acc)))
    },
    mc.preschedule=F, mc.cores=50)
  )
}

uniq.annotations <- function(annotation.matrix, type) {
  sort(
    unique(do.call('c', (annotation.matrix[, type]))),
    na.last=F
    )
}
