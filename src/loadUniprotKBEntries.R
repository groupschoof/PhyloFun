library(XML)
library(utils)
library(RCurl)
library(parallel)
library(biomaRt)
library(stringr)

EVIDENCE.CODES <- c( 'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP' )

uniprotkb.url <- function( accession, frmt='xml' ) {
  # Returns valid URL to access Uniprot's RESTful Web-Service to download
  # data about the Protein as referenced by the argument 'accession'.
  # Note, that the accession is URL encoded before being pasted into the
  # Uniprot URL template.
  #
  # Args:
  #  accession : The Protein's Uniprot accession.
  #  frmt      : The format of the downloaded Uniprot Entry. Default is 'xml'.
  #
  # Returns: The Uniprot URL for the argument accession.
  #   
  paste(
    'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/',
    URLencode( accession ),
    '/', frmt, sep=''
  )
}

extractAnnotations <- function( doc, type, xpath.prefix='//' ) {
  # Parses XLM document and extracts the 'id' attributes of all
  # 'dbReference' tags where the 'type' attribute is set to the provided
  # argument.
  #
  # Args:
  #  doc : Downlaoded and parsed XMl document as returned by 'xmlInternalTreeParse'
  #  type : The type of 'dbReference' tags to extract. I.e. 'InterPro'.
  #
  # Returns: A character vector of all 'id' attributes of the
  # 'dbReference' tags having their 'type' attribute set to this
  # functions argument 'type'. 
  #   
  block <- function() {
    ns <- c( xmlns="http://uniprot.org/uniprot" )
    ndst <- getNodeSet( doc,
      paste( xpath.prefix, "xmlns:dbReference[@type='", type, "']", sep="" ),
      namespaces=ns
    )
    vapply( ndst, xmlGetAttr, vector( mode='character', length=1 ), 'id' )
  }
  tryCatch( block(), error=function( err ) {
    warning( err, " caused by parsed XML Uniprot document ", doc )
  })
}

retrieveAnnotations <- function( uniprot.entry,
    annotation.types=c( 'GO','InterPro','Pfam' ), xpath.prefix=''
  ) {
  # Parses the already downloaded XMl representation of a protein identified by
  # their uniprot.accessions. The specified annotation.types encoded in the
  # downloaded documents are then extracted and returned.
  #
  # Args:
  #   uniprot.entry : A single 'entry' XML node as a member of the result of
  #                   calling downloadUniprotDocuments(…) 
  #   annotation.types : Type of annotations to extract from the downloaded
  #                      document.
  #   xpath.prefix : The prefix to append the XPath query, passed to method
  #                  extractAnnotations(…).
  #
  # Returns: list with keys as in argument annotation.types and values vectors
  # containing the annotated domain IDs of the corresponding types.
  # Exmpl:
  # $GO
  # [1] "GO:0005737" "GO:0005634" "GO:0005524" "GO:0008026" "GO:0003725"
  # 
  # $InterPro
  # [1] "IPR005034" "IPR011545" "IPR001159" "IPR014001" "IPR001650" "IPR003100"
  # 
  # $Pfam
  # [1] "PF00270" "PF03368" "PF00271" "PF02170" "PF00636"
  #   
  if( ! is.null( uniprot.entry ) ) {
    setNames(
      lapply( annotation.types,
        function( type ) {
          extractAnnotations( uniprot.entry, type, xpath.prefix=xpath.prefix )
        }
      ),
      annotation.types
    )
  }
}

retrieveUniprotAnnotations <- function( uniprot.accessions,
  annotation.types=c( 'GO','InterPro','Pfam' ) ) {
  # Queries the Uniprot web service to download the XML encoded available
  # information for the provided uniprot.accessions. These are parsed to
  # extract functional annotations of the given annotation.types.
  #
  # Args:
  #  uniprot.accessions : The character vector of valid Uniprot accessions.
  #                       I.e. c( 'sp|Q5ZL72|CH60_CHICK', 'P63039' ) 
  #  annotation.types   : The types of function annotation to extract. Default
  #                       is c( 'GO','InterPro','Pfam' )
  #
  # Returns: A matrix of mode character. Row names are annotation.types and
  # column names are the uniprot.accessions. Each cell contains as a character
  # vector the parsed annotations of the repective type for the respective
  # Uniprot accession. If no document could be found for a given accession the
  # respective column is omitted. 
  #   
  entries <- downloadUniprotDocuments( uniprot.accessions )
  do.call( 'cbind',
    lapply( entries, function( ue ) {
      retrieveAnnotations( ue, annotation.types=annotation.types )
    })
  )
}

uniq.annotations <- function( annotation.matrix, type, exclude.NAs=FALSE ) {
  # Looks up all pairwise distinct annotations of argument type. NAs can be
  # excluded.
  #
  # Args:
  #  annotation.matrix : The matrix of protein annotations as returned by
  #                      function 'retrieveAnnotationsBiomart'.
  #  type              : The row name of 'annotation.matrix', i.e. 'GO' or
  #                      'InterPro'.
  #  exclude.NAs       : Switch indicating wether to include NA values in the
  #                      result.
  #
  # Returns: an alphabetically sorted vector of all pairwise distinct
  # annotations.
  #   
  ua <- sort(
    unique(do.call('c', (annotation.matrix[type,]))),
    na.last=F
  )
  if ( exclude.NAs ) {
    ua[ ! is.na(ua[]) ]
  } else {
    ua 
  }
}

retrieveAnnotationsBiomart <- function( accs,
    uni.mart=useDataset("uniprot",mart=useMart("unimart"))
  ) {
  # Uses library biomaRt to retrieve InterPro and Gene Ontology annotations for
  # the Proteins referenced by their accessions in argument 'accs'. Results are
  # returned as a matrix where the rows are the Protein accessions and the
  # columns 'InterPro' and 'GO'. NOTE: This function cannot retrieve Pfam
  # annotations.
  #
  # Args:
  #  accs       : The accessions of the proteins to retrieve annotations for.
  #  uni.mart   : The biomaRt mart object to use.
  #
  # Returns: A matrix with all InterPro and GO annotations for the query
  # proteins. Columns are the protein accessions and rows 'InterPro' and 'GO'.
  #   
  annos <- getBM( c("accession", "interpro_id", "go_id"),
    filters=c("accession"), values=accs, mart=uni.mart )
  do.call( 'cbind',
    lapply( accs, function(a) {
      matrix(
        list( 
          'InterPro'=unique( annos[ annos["accession"] == a, "interpro_id" ] ),
          'GO'=unique( annos[ annos["accession"] == a, "go_id" ] )
        ),
        nrow=2, ncol=1,
        dimnames=list( c('InterPro', 'GO'), a )
      )
    })
  )
}

extractName <- function( uniprot.entry, xpath.prefix='./', noverbose=T,
  return.error=F ) {
  # Finds and returns the content of the first accession tag inside the
  # argument 'uniprot.entry'.
  #
  # Args:
  #  uniprot.entry : Either the character content of a Uniprot XML document or
  #                  a _single_ result of function 'getEntries'.
  #  xpath.prefix  : Depends on type of argument 'uniprot.entry'. If first
  #                  option was chosen, use '//xmlns:entry/', if second option
  #                  was chosen, use the default './'.
  #  noverbose     : Set to FALSE, if error messages should be printed to
  #                  stderr.
  #  return.error  : If TRUE and an error occured the message will be returned,
  #                  in case of an error and this switch set to FAlSE NULL will
  #                  be returned.
  #
  # Returns: A character, content of the first accession tag in argument
  # 'uniprot.entry', or NULL.
  #   
  block <- function() {
    if ( ! is.null(uniprot.entry) ) {
      ns <- c(xmlns="http://uniprot.org/uniprot");
      if ( is.character(uniprot.entry) ) {
        uniprot.entry <- xmlInternalTreeParse( uniprot.entry )
      }
      xmlValue(
        getNodeSet(
          uniprot.entry,
          paste( xpath.prefix, "xmlns:accession", sep='' ),
          ns
        )[[ 1 ]]
      )
    }
  }
  accession <- try( block(), silent=noverbose )
  if ( class( accession ) == "try-error" && ! return.error ) {
    accession <- NULL
  }
  accession
}

downloadUniprotDocuments <- function( uniprot.accessions, frmt='xml',
  uniprot.webfetch.max.ids=200 ) {
  # Downloads the documents from the RESTful Uniprot web service.
  #
  # Args:
  #  uniprot.accessions       : A vector or list of uniprot accessions
  #  frmt                     : The format in which the uniprot documents shall
  #                             be encoded
  #  uniprot.webfetch.max.ids : The current maximum number of IDs in a batch
  #                             fetch allowed by the Uniprot webfetch service.
  #
  # Returns: Returns a named character vector of the downloaded XML documents.
  # Names are the 'uniprot.accessions'.
  if ( length( uniprot.accessions ) > uniprot.webfetch.max.ids ) {
    # Recursive fetch of uniprot.webfetch.max.ids sized batches:
    c( 
      downloadUniprotDocuments( 
        uniprot.accessions[ 1:uniprot.webfetch.max.ids ] ),
      downloadUniprotDocuments(
        uniprot.accessions[
          ( uniprot.webfetch.max.ids + 1 ):length(uniprot.accessions)
        ]
      )
    )
  } else {
    # Fetch max uniprot.webfetch.max.ids in a single batch:
    fetch.url <- uniprotkb.url(
      paste( uniprot.accessions, collapse=",", sep=""),
      frmt=frmt
    )
    uniprot.entries <- getEntries( getURL( fetch.url ) )
    if ( ! is.null( uniprot.entries ) )
      names( uniprot.entries ) <- unlist(
        lapply( uniprot.entries, extractName )
      )
    uniprot.entries
  }
}

retrieveSequence <- function( doc, noverbose=T, return.error=T ) {
  # Parses the XML Uniprot document 'doc' and returns the content of the
  # contained amino acid sequence.
  #
  # Args:
  #  doc : The document downloaded from Uniprot using the URL as generated by
  #        function 'uniprotkb.url'.
  #  noverbose : If TRUE error messages won't be printed to stderr.
  #  return.error : If TRUE and an error occured the message will be returned,
  #                 in case of an error and this switch set to FAlSE NULL will
  #                 be returned.
  #
  # Returns: The content of the contained '<sequence>...</sequence>' tag.
  #   
  block <- function() {
    if ( ! is.null(doc) ) {
      ns <- c(xmlns="http://uniprot.org/uniprot");
      if ( is.character(doc) ) {
        doc <- xmlInternalTreeParse( doc )
      }
      gsub("\\s", "",
        xmlValue(
          getNodeSet( doc, "//xmlns:entry/xmlns:sequence", ns)[[1]]
        )
      )
    }
  }
  seq <- try( block(), silent=noverbose )
  if ( class( seq ) == "try-error" && ! return.error ) {
    seq <- NULL
  }
  seq
}

retrieveSequences <- function( downloaded.uniprot.docs, accessions=names( downloaded.uniprot.docs ),
  max.retries=10, noverbose=F ) {
  # Parses each downloaded argument Uniprot document and extracts the content
  # of the '<sequence>' tag. If an error occurs doing so, re-tries downloading
  # and parsing the document, after sleeping a random amount of time out of
  # 1:90 seconds.
  #
  # Args:
  #  downloaded.uniprot.docs : The documents downloaded from Uniprot's RESTful
  #                            web service. 
  #  accessions : The names of the returned list. Default names( downloaded.uniprot.docs )
  #  max.retries : The maximum number of times another download is attempted.
  #  noverbose : switch indicating wether to print out errors
  #
  # Returns: A named vector of extracted sequences or error messages.
  #   
  seqs <- sapply( downloaded.uniprot.docs, function( doc ) {
      try( retrieveSequence( doc ), silent=noverbose )
    })
  err.indxs <- grepl("^Error", seqs[], perl=T) 
  err.uris <- names( seqs[ err.indxs ] )
  if ( length(err.uris) > 0 && max.retries > 0 ) {
    Sys.sleep( sample(1:90, 1) )
    print( paste( "Retry number", as.character( 11 - max.retries) ) )
    seqs <<- c(
      seqs[ ! err.indxs ],
      retrieveSequences( downloadUniprotDocuments( err.uris ),
        (max.retries - 1) )
    )
  }
  names( seqs ) <- accessions
  # Return
  seqs
}

sharedFunction <- function( annotation.matrix,
  accessions=colnames(annotation.matrix), annotation.type='GO' ) {
  # Constructs a Boolean matrix with column names as in the argument
  # 'annotation.matrix' columns. Rows are the argument accessions. Each matrix
  # entry is TRUE if and only if the two corresponding proteins share at least
  # a single annotation. 
  #
  # Args:
  #  annotation.matrix : Matrix holding different function annotation for its
  #                      column proteins
  #  accessions : The accessions to compute the Boolean entries for, default
  #               all proteins in annotation.matrix
  #  annotation.type : The type of functional annotation to use, default 'GO'
  #
  # Returns: Boolean matrix with rows as in argument 'accessions' and columns
  # the same as in argument 'annotation.matrix'.
  #   
  do.call( 'cbind',
    setNames(
      mclapply( colnames(annotation.matrix), function( a ) {
          sapply( accessions, function( b ) {
              length(
                intersect(
                  annotation.matrix[[ annotation.type, a ]],
                  annotation.matrix[[ annotation.type, b ]]
                )
              ) > 0
          } )
      }, mc.cores=detectCores(), mc.preschedule=T ),
      colnames(annotation.matrix)
    )
  )
}

sharedAnnotation <- function( annotation.matrix,
  annotation, annotation.type='GO' ) {
  # Filters argument 'annotation.matrix' where cells in row 'annotation.type'
  # include the argument 'annotation'. 
  #
  # Args:
  #  annotation.matrix : The matrix of functional annotations as returned by
  #                      i.e. 'retrieveAnnotationsBiomart'.
  #  annotation        : The annotation to filter for, i.e. "GO:0005634"
  #  annotation.type   : The row of the annotation.matrix to check entries,
  #                      default 'GO'.
  #
  # Returns: Returns matching sub-matrix.
  #   
  annotation.matrix[ ,
    mapply(
      function(x){ any( x == annotation & ! is.na(x) ) },
      annotation.matrix[ annotation.type, ]
    ),
    drop=F
  ]
}

intersectAnnotations <- function( annotation.matrix, acc.a, acc.b, annotation.type="GO" ) {
  # Intersects the annotations of protein 'acc.a' with protein 'acc.b' using
  # the argument 'annotation.type'.
  #
  # Args:
  #  annotation.matrix : The annotations as returned i.e. by
  #                      'retrieveAnnotationsBiomart'
  #  acc.a : Accession of first protein 
  #  acc.b : Accession of second protein 
  #
  # Returns: The set intersection of the two proteins' annotation sets.
  #   
  intersect( annotation.matrix[[ annotation.type, acc.a ]],
    annotation.matrix[[ annotation.type, acc.b ]]
  )
}

shareAnnotation <- function( annotation, annotation.matrix, acc.a, acc.b,
  annotation.type="GO" ) {
  # Looks up the annotations of both arguments 'acc.a' and 'acc.b' and checks
  # if both have argument 'annotation' in their 'annotation.type' set. 
  #
  # Args:
  #  annotation : The GO term, InterPro ID, or Pfam ID to look up.
  #  annotation.matrix : The annotations as returned i.e. by
  #                      'retrieveAnnotationsBiomart'
  #  acc.a : Accession of first protein 
  #  acc.b : Accession of second protein 
  #
  # Returns: TRUE if and only if, both proteins share the argument
  # 'annotation', FALSE otherwise.
  #   
  annotation %in% intersectAnnotations( annotation.matrix, acc.a, acc.b,
    annotation.type
  )
}

extractExperimentallyVerifiedGoAnnos <- function( doc, xpath.prefix='//',
  evidence.codes=EVIDENCE.CODES ) {
  # Uses XPath to extract those GO annotations that are experimentally
  # verified. Note, that warnings generated by calls to the XML library are
  # suppressed to not confuse the user, when no experimentally verified GO
  # annotations could be found.
  #
  # Args:
  #  doc            : A XML tag of type entry as returned i.e. by function
  #                   'getEntries'
  #  xpath.prefix   : Define the scope in which to apply the XPath search. '//'
  #                   means whole document, './' means children of the current
  #                   node.
  #  evidence.codes : The evidence codes accpeted as 'experimentally verified',
  #                   default all *truly* experimentally verified annotations.
  #                   Add 'TAS' and 'IC' to also accept 'Traceable Author
  #                   Statement' and 'Inferred by Curator', respectively. See
  #                   http://www.geneontology.org/GO.evidence.shtml for more.
  #
  # Returns: A 1*1 matrix containing a character vector of the extracted
  # experimentally verified GO annotations. Row name is 'GO' and column name is
  # the proteins FIRST accession as appearing in the XML document. Returns
  # NULL, if no experimentally verified Go annotations can be found.
  #   
  block <- function() {
    ns <- c( xmlns="http://uniprot.org/uniprot" )
    xpath.ev.cds <- paste( lapply( evidence.codes, function( ec ) paste( "contains(@value, '", ec, "')", sep="" ) ), collapse=" or " )
    xpath.query <- paste( xpath.prefix, 
      "xmlns:dbReference[@type='GO']//xmlns:property[@type='evidence' and ", 
      "( ", xpath.ev.cds, " ) ]/..",
      sep=''
    )
    ndst <- suppressWarnings( getNodeSet( doc, xpath.query, namespaces=ns ) )
    acc <- suppressWarnings(
      xpathApply( doc, paste( xpath.prefix, "xmlns:accession", sep="" ),
        xmlValue, namespaces=ns )[[1]]
    )
    if ( ! is.null( ndst ) && length( ndst ) > 0 ) {
      mtrx <- matrix( list(), ncol=1, nrow=1, dimnames=list( 'GO', acc ) )
      mtrx[[ 1, 1 ]] <- vapply( ndst, xmlGetAttr, vector( mode='character', length=1 ), 'id' )
      mtrx
    } else
      NULL
  }
  tryCatch( block(), error=function( err ) {
    warning( err, " caused by document ", doc )
  })
}

getEntries <- function( uniprot.xml, uniprot.error.msg.regex='^ERROR' ) {
  # Uniprot's dbfetch can be asked to return several entry tags in the same XML
  # document. This function uses XPath queries to extract all complete uniprot
  # tags. 
  #
  # Args:
  #  uniprot.xml             : The result of a web fetch to Uniprot i.e. using
  #                            getURL.
  #  uniprot.error.msg.regex : A regular expression to avoid parsing an error
  #                            returned from Uniprot.
  #
  # Returns: A list of extracted uniprot-entry-tags as returned by function
  # 'getNodeSet'. 
  #   
  if ( ! is.null( uniprot.xml ) && '' != uniprot.xml &&
    ! grepl( uniprot.error.msg.regex, uniprot.xml )
  ) {
    ns <- c( xmlns="http://uniprot.org/uniprot" )
    getNodeSet( 
      xmlInternalTreeParse( uniprot.xml ), 
      "//xmlns:entry", namespaces=ns
    )
  } else {
    NULL
  }
}

retrieveExperimentallyVerifiedGOAnnotations <- function( uniprot.accessions,
  uniprot.webfetch.max.ids=200, evidence.codes=EVIDENCE.CODES ) {
  # Downloads and parses XML documents from Uniprot for each accession in
  # argument. Extracts all experimentally verified GO annotations.
  #
  # Args:
  #  uniprot.accessions       : A character vector of Uniprot accessions.
  #  uniprot.webfetch.max.ids : The current maximum number of IDs in a batch
  #                             fetch allowed by the Uniprot webfetch service.
  #
  # Returns: A matrix with row 'GO' and one column for each Uniprot accession.
  # Each cell is either NULL or a character vector holding all experimentally
  # verified GO annotations. NULL annotations are excluded, so the returned
  # matrix can be of zero columns and a single row.
  #   
  uniprot.entries <- downloadUniprotDocuments( uniprot.accessions )  
  if ( ! is.null(uniprot.entries) && length( uniprot.entries ) > 0 ) {
    annos <- do.call( 'cbind',
      lapply( uniprot.entries , function( d ) {
        extractExperimentallyVerifiedGoAnnos( d, xpath.prefix='./',
          evidence.codes=evidence.codes )
      })
    )
    # Exclude NULL columns:
    annos[ , as.character( annos[ 'GO', ] ) != 'NULL' , drop=F ]
  }
}

extractRefSeqAccession <- function( ref.seq.prot.name,
  acc.regex='^[^\\|]+\\|[^\\|]+\\|[^\\|]+\\|([^\\|]+)\\|' ) {
  str_match( ref.seq.prot.name, acc.regex )[[ 1, 2 ]]
}
