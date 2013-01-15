parseJackhmmerTable <- function( jr,
  skip.lines=3, parse.line=list( '1'='hit.name', '3'='query.name' ) ) {
  # Parses the table output of HMMER-3's Jackhmmer. Reads out hit accessions
  # and query accessions.
  #
  # Args:
  #  jr : List of lines of the tabular Jackhmmer output. As a result of
  #       function 'scan(...)'
  #  skip.lines : Skip so many lines of the argument list 'jr', default: 3.
  #  parse.line : Read out position 'i' and return as list entry, named as
  #               value of 'i'. I.e. position '1' is read out as 'hit.name'. 
  #
  # Returns: Matrix with colnames as the values in argument 'parse.line' and
  # one row per line in jr, the first 'skip.lines' are ignored.
  #   
  jr.body <- jr[ ( skip.lines + 1 ) : length(jr) ]
  do.call( 'rbind',
    lapply( jr.body, function( ln ) {
      split.line <- strsplit( ln, "\\s+", perl=T )[[1]]
      setNames(
        lapply( names( parse.line ), function( pos ) {
          split.line[ as.integer( pos ) ]
        }),
      parse.line )
    })
  )
}

extractUniprotAccessionFromUniprotName <- function( 
  uniprot.name ) {
  # Extracts the string in between "|" and returns it.
  #
  # Args:
  #  uniprot.name : The complete Uniprot name of a protein as encountered in
  #                 the swissprot or trEMBL fasta databases.
  #                 I.e. "sp|B5YXA4|DNAA_ECO5E" 
  #
  # Returns: "B5YXA4" for input "sp|B5YXA4|DNAA_ECO5E".
  #   
  sub( "\\S+\\|(\\S+)\\|\\S+", "\\1", uniprot.name )
}

parseInterProScanTable <- function( ipr.scn.lines,
  prot.acc.regex="^(\\S+)\\s+.*", ipr.regex=".*(IPR\\d{6}).*" ) {
  # Each line in the InterProScan result that contains an InterPro Entry ID is
  # coerced into a an annotation matrix.
  #
  # Args:
  #  ipr.scn.lines : Lines of the InterProScan result file to parse.
  #  prot.acc.regex : Regular expression to match the Query Protein accession
  #                   in an InterProScan result line.
  #  ipr.regex : Regular expression to match the InterPro Entry ID in an
  #              InterProScan result line.
  #
  # Returns: A Matrix with the Proteins' InterPro annotations, where the
  # columns are the protein accessions and a single row "InterPro" points to
  # matrix cells with lists of InterPro Entry IDs for each protein.
  #   
  ipr.res <- list()
  for ( ln in ipr.scn.lines ) {
    if( grepl( ipr.regex, ln ) ) {
      prot.acc <- sub( prot.acc.regex, "\\1", ln )
      ipr.id   <- sub( ipr.regex, "\\1", ln )
      if( ! ipr.id %in% ipr.res[[ prot.acc ]]$InterPro ) {
        ipr.res[[ prot.acc ]]$InterPro <- append( 
          ipr.res[[ prot.acc ]]$InterPro,
          ipr.id
        )
      }
    }
  }
  do.call( 'cbind', ipr.res )
}

commandLineArguments <- function( trailing.args, default.args ) {
  # Reads out command line arguments like '-foo bar -franky jr' passed as
  # argument trailing.args. Merges named lists trailing.args and default.args
  # giving precedence to arguments in trailing.args over default.args.
  #
  # Args:
  #  trailing.args : Character, Rscript command line arguments as returned by
  #                  commandArgs(trailingOnly = TRUE).
  #  default.args  : Named list of required default arguments, to fill in for
  #                  missing arguments in trailing.args.
  #
  # Returns: A named list of arguments for the R code to be run.
  #   
  sub.funk <- function( a ) sub('^-', '', a)
  arg.names <- sapply( trailing.args[ which( grepl( '^-', trailing.args[] ) ) ],
    sub.funk, USE.NAMES=F
  ) 
  arg.values <- sapply( trailing.args[ which( grepl( '^[^-]', trailing.args[] ) ) ],
    sub.funk, USE.NAMES=F
  ) 
  names( arg.values ) <- arg.names
  all.arg.names <- unique( c( arg.names, names( default.args ) ) )
  setNames( 
    lapply( all.arg.names, function( arg.nm ) {
      if ( arg.nm %in% arg.names ) {
        arg.values[[ arg.nm ]]
      } else if ( arg.nm %in% names( default.args ) ) {
        default.args[[ arg.nm ]]
      }
    }),
    all.arg.names
  )
}

sanitizeUniprotAccession <- function( protein.name ) {
  # Sanitized protein accessions to be used with PhyloFun's pipeline programs,
  # i.e. 'GBlocks'. Trims whitespaces and extracts word between pipes, if pipes
  # are found.
  #
  # Args:
  #  protein.name : The name of the protein, i.e. sp|MyAccession|MOUSE_SHEEP
  #
  # Returns: Returns the sanitied protein accession as character vector of
  # length one, i.e. 'MyAccession'.
  #   
  if ( is.null( protein.name ) )
    return( protein.name )
  no.blanks <- '^\\s*(\\S+)\\s*'
  ua <- str_match( protein.name, no.blanks )[[ 1, 2 ]]
  between.pipes <- '\\S+\\|(\\S+)\\|\\S+'
  if ( grepl( between.pipes, ua, perl=T ) ) {
    str_match( ua, between.pipes )[[ 1, 2 ]]
  } else {
    ua
  }
}
