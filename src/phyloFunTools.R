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

uniqueHomologs <- function( path.2.homlgs.fasta,
  path.2.unique.hmlgs.fasta=path.2.homlgs.fasta, print.warning=T ) {
  # Looks up duplicated accessions in file 'path.2.homlgs.fasta' and removes
  # those duplicated entries, saving the result into
  # 'path.2.unique.hmlgs.fasta'.
  #
  # Args:
  #  path.2.homlgs.fasta       : The file path to the FASTA input.
  #  path.2.unique.hmlgs.fasta : The file path to the FASTA output. Default is
  #                              to overwrite input.
  #  print.warning             : If TRUE duplicated accessions are printed out
  #                              in a warning message.
  #
  # Returns: TRUE if no error occurred.
  #   
  hmlgs <- read.AAStringSet( path.2.homlgs.fasta )
  dplcts <- duplicated( names( hmlgs ) )
  if ( any( dplcts ) ) {
    uniq.hmlgs <- hmlgs[ - which( dplcts ) ]
    write.XStringSet( uniq.hmlgs, path.2.unique.hmlgs.fasta )
    if ( print.warning ) {
      warning( "Removed duplicated Accessions:",
        names( hmlgs )[ which( dplcts ) ]
      )
    }
  }
  T
}

filterMultipleSequenceAlignment <- function( msa.xstring.set, min.chars=30 ) {
  # Discards all Amino Acid sequences shorter than 30 Amino Acids, that is NOT
  # counting gaps.
  #
  # Args:
  #  msa.xstring.set : An object of class XStringSet holding all AA seqs of the
  #                    MSA ( see package Biostrings for details on
  #                    readAAStringSet() )
  #  min.chars       : The minimum number of Amino Acids to be present in any
  #                    sequence of the MSA to be retained after filtering.
  #
  # Returns: The filtered MSA from which all AA sequences with less then
  # min.chars AAs have been discarded. Or null, if the whole MSA is shorter
  # than min.chars.
  #   

  # In a MSA all sequences have the same length:
  if ( nchar( msa.xstring.set[ 1 ] ) < min.chars ) {
    NULL # return
  } else {
    msa.aa.lens <- lapply( msa.xstring.set, function( aa.seq ) {
      nchar( gsub( '-|\\s', '', toString( aa.seq ) ) )
    })
    msa.xstring.set[ which( msa.aa.lens[] > min.chars ) ]
  }
}

chooseAlignment <- function( msa.unfiltered,
  msa.filtered=filterMultipleSequenceAlignment( msa.unfiltered ),
  min.fraction=0.5, min.no.seqs=2 ) {
  # Chooses the Multiple Sequence Alignment to be used in further analyses.
  # Discards the filtered MSA, if it is NULL or has not at least min.no.seqs or
  # has just retained less than min.fraction sequences of the original
  # msa.unfiltered. 
  #
  # Args:
  #  msa.unfiltered : The unfiltered MSA as XStringSet
  #  msa.filtered   : The filtered MSA as XStringSet
  #  min.fraction   : The minimum fraction of the number of sequences that
  #                   msa.filtered needs to be selectable.  
  #  min.no.seqs    : The minimum number of sequences msa.filtered has to have
  #                   to be selectable.
  #
  # Returns: The filtered MSA, if it meets the above criteria, or the
  # unfiltered one.
  #   
  if ( is.null( msa.filtered ) ||
    length( msa.filtered ) < min.no.seqs || 
    length( msa.filtered ) / length( msa.unfiltered ) < min.fraction
  ) {
    msa.unfiltered
  } else {
    msa.filtered
  }
}

filterPhylogeneticTree <- function( phyl.tree, leaves.to.retain,
  leaves.to.check=1:length( phyl.tree$tip.label ), leave.annos,
  annotation.type='GO' ) {
  curr.leaf.parnt <- phyl.tree$edge[
    phyl.tree$edge[ , 2 ] == leaves.to.check[[ 1 ]], 1, drop=T
  ][[ 1 ]]
  zero.len.edges <- phyl.tree$edge[
    which( phyl.tree$edge.length[] == 0.0 ), , drop=F
  ]
  curr.leaves <- zero.len.edges[
    which( zero.len.edges[ , 1 ] == curr.leaf.parnt ), 2
  ]
  curr.annos <- leave.annos[
    annotation.type,
    intersect( colnames( leave.annos ), phyl.tree$tip.label[ curr.leaves ] ),
    drop=F
  ]
  curr.anno.space <- annotationSpace( curr.annos )
  curr.leave.labels.to.retain <- union(
    as.character(
      lapply( curr.anno.space, function( anno.spc ) {
        leave.ind <- which( as.logical(
          lapply(
            curr.annos[ annotation.type, ],
            function( iter.anno ) identical( iter.anno, anno.spc )
          )
        ) )[[ 1 ]]
       colnames( curr.annos[ annotation.type, leave.ind ] )
      })
    ),
    intersect( leaves.to.retain, phyl.tree$tip.label[ curr.leaves ] )
  )
}
