library(XML)
library(utils)
library(RCurl)
library(parallel)
library(biomaRt)
library(stringr)

EVIDENCE.CODES <- c( 'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC' )

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

downloadUniprotDocuments <- function( uniprot.accessions, 
  uniprot.webfetch.max.ids=200 ) {
  # Downloads the documents from the RESTful Uniprot web service.
  #
  # Args:
  #  uniprot.accessions       : A vector or list of uniprot accessions
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
      paste( uniprot.accessions, collapse=",", sep=""), 'xml'
    )
    uniprot.entries <- getEntries( getURL( fetch.url ) )
    if ( ! is.null( uniprot.entries ) )
      names( uniprot.entries ) <- unlist(
        lapply( uniprot.entries, extractName )
      )
    uniprot.entries
  }
}

downloadSequences <- function( uniprot.accessions, fasta.file,
  uniprot.webfetch.max.ids=200, sanitize.uniprot.accessions=TRUE ) {
  # Uses the Uniprot web services to download each FASTA document for the
  # argument uniprot.accessions. 
  # 
  # Args:
  #  uniprot.accessions : The valid Uniprot accessions to download the
  #               amino acid sequences for, i.e. c(
  #               'sp|O08601|MTP_MOUSE', 'P55158' ).
  #  fasta.file : The path to the FASTA file the downloaded sequences will be
  #               stored in.
  #  sanitize.uniprot.accessions : If set to TRUE, the default, the accessions
  #               of the downloaded sequences will be sanitized using the
  #               'sanitizeUniprotAccession(…)' function.
  #
  # Returns: TRUE if and only if no error has occurred.
  #   
  if ( length( uniprot.accessions ) > uniprot.webfetch.max.ids ) {
    # Recursive fetch of uniprot.webfetch.max.ids sized batches:
    all(
      c(
        downloadSequences( 
          uniprot.accessions[ 1:uniprot.webfetch.max.ids ], fasta.file,
          uniprot.webfetch.max.ids
        ),
        downloadSequences(
          uniprot.accessions[
            ( uniprot.webfetch.max.ids + 1 ):length(uniprot.accessions)
          ], fasta.file, uniprot.webfetch.max.ids
        )
      )
    )
  } else {
    # Fetch max uniprot.webfetch.max.ids in a single batch:
    fetch.url <- uniprotkb.url(
      paste( uniprot.accessions, collapse=",", sep=""), 'fasta'
    )
    fastas <- getURL( fetch.url ) 
    if ( ! is.null( fastas ) && ! is.na( fastas ) && length( fastas ) > 0 ) {
      # Sanitize Uniprot Accessions, if requested:
      if ( sanitize.uniprot.accessions ) {
        fastas <- paste( as.character(
          lapply( strsplit( fastas, '\n', fixed=TRUE )[[1]],
            function(x) {
              if( grepl( '^>', x )[[1]] ) {
                paste( '>', sanitizeUniprotAccession( x ), sep='' )
              } else {
                x
              }
            }
          )
        ), collapse="\n" )
      }
      # Append downloaded sequences to file:
      write( sub( '\\n$', '', fastas ), file=fasta.file, append=T )
      TRUE
    }
  }
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
  #                   If set to NULL, NA, 'all', or 'ALL', no filtering for
  #                   specific evidence codes will be performed and all
  #                   available GO annotations will thus be returned.
  #
  # Returns: A three column data.frame where each row represents a single
  # protein GO term annotation. The first column of holds the GO term
  # accessions, the second the evidence codes, and the final third column holds
  # the protein accessions.
  #
  block <- function() {
    ns <- c( xmlns="http://uniprot.org/uniprot" )
    xpath.ev.cds <- if ( is.null( evidence.codes ) || is.na( evidence.codes )
      || evidence.codes == 'ALL' || evidence.codes == 'all' ) {
      ''
    } else {
      paste(  ' and (',
        paste( lapply( evidence.codes, function( ec ) {
            paste( "contains(@value, '", ec, "')", sep="") 
          } ),
          collapse=" or "
        ),
        ') '
      )
    }
    xpath.query <- paste( xpath.prefix, 
      "xmlns:dbReference[@type='GO']//xmlns:property[@type='evidence' ", 
      xpath.ev.cds, " ]/..",
      sep=''
    )
    ndst <- suppressWarnings( getNodeSet( doc, xpath.query, namespaces=ns ) )
    acc <- suppressWarnings(
      xpathApply( doc, paste( xpath.prefix, "xmlns:accession", sep="" ),
        xmlValue, namespaces=ns )[[1]]
    )
    if ( ! is.null( ndst ) && length( ndst ) > 0 ) {
      as.data.frame( cbind(
        as.character( lapply( ndst, xmlGetAttr, 'id' ) ),
        as.character( lapply( ndst, getEvidenceCode ) ),
        'V3'=acc
      ), stringsAsFactors=FALSE )
    } else
      NULL
  }
  tryCatch( block(), error=function( err ) {
    warning( err, " caused by document ", doc )
  })
}

getEvidenceCode <- function( db.ref.tag, xpath.prefix='./' ) {
  # Uses XPATH to find and extract the GO annotation 'db.ref.tag' evidence
  # code.
  #
  # Args:
  #  db.ref.tag   : An instance of class 'XMLAbstractNode' representing the
  #                 Uniprot 'dbReference' tag.
  #  xpath.prefix : The prefix to put at the beginning of the XPATH query,
  #                 default is './'.
  #
  # Returns: Returns the extracted and parsed evidence code, shortened to the
  # three letter abbreviation. For example 'IEA'. If no matching evidence code
  # can be found or argument db.ref.tag is NULL the return value of this
  # function is NULL.
  #   
  if ( ! is.null( db.ref.tag ) ) {
    ns <- c( xmlns="http://uniprot.org/uniprot" ) 
    ev.ns <- getNodeSet( db.ref.tag,
      paste( xpath.prefix, "xmlns:property[@type='evidence']", sep="" ),
      namespace=ns )
    if ( length( ev.ns ) > 0 ) {
      ev.cd <- xmlGetAttr( ev.ns[[ 1 ]], 'value' )
      sub( ':.*$', '', ev.cd )
    } else {
      NA
    }
  } else {
    NA
  }
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
    do.call( 'rbind',
      lapply( uniprot.entries , function( d ) {
        extractExperimentallyVerifiedGoAnnos( d, xpath.prefix='./',
          evidence.codes=evidence.codes )
      })
    )
  }
}

extractRefSeqAccession <- function( ref.seq.prot.name,
  acc.regex='^[^\\|]+\\|[^\\|]+\\|[^\\|]+\\|([^\\|]+)\\|' ) {
  str_match( ref.seq.prot.name, acc.regex )[[ 1, 2 ]]
}
