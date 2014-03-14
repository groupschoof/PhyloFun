proteinPairsSharingAnnotation <- function( annotation, protein.pairs.tbl,
  annotation.tbl, pairs.first.col=1, pairs.secnd.col=2, annot.prot.col=3,
  annot.col=1 ) {
  # Returns unique pairs where at least one member has the annotation
  # 'annotation' an additional boolean column indicates wether both proteins
  # share the annotation or only a single member of the pair is annotated with
  # it.
  #
  # Args:
  #  annotation        : The annotation term to lookup protein pairs for, i.e.
  #                      "GO:0001234"
  #  protein.pairs.tbl : The table of protein pairs sharing a significant
  #                      sequence similarity
  #  annotation.tbl    : The table of protein annotations. Format should be
  #                      like the one used by GOstats, a three column matrix,
  #                      in which the first column holds GO terms, the second
  #                      evidence codes, and the third column protein
  #                      accessions.
  #  pairs.first.col   : The column of 'protein.pairs.tbl' in which to lookup
  #                      the first member of protein pairs.
  #  pairs.secnd.col   : The column of 'protein.pairs.tbl' in which to lookup
  #                      the second member of protein pairs.
  #  annot.prot.col    : The column of 'annotation.tbl' in which to lookup
  #                      protein accessions.
  #  annot.col         : The column of 'annotation.tbl' in which to lookup GO
  #                      term accessions.
  #
  # Returns: A subset of protein.pairs.tbl in which each pair has at least one
  # member (protein) annotated with argument 'annotation'. Can be an empty
  # subset (matrix with same number of columns as argument 'protein.pairs.tbl'
  # and 0 rows). The returned table holds an additional boolean column
  # indicating wether both protein pair members are annotated with the argument
  # 'annotation' or not.
  #   
  annot.ps <- annotation.tbl[ which( annotation.tbl[ , annot.col ] ==
                                    annotation ), annot.prot.col ]
  prs <- protein.pairs.tbl[ which(
    protein.pairs.tbl[ , pairs.first.col ] %in% annot.ps |
    protein.pairs.tbl[ , pairs.secnd.col ] %in% annot.ps
  ), , drop=FALSE ]

  lhs <- prs[ , 1 ] %in% annot.ps
  rhs <- prs[ , 2 ] %in% annot.ps

  cbind( prs, 'annotation.shared'=lhs & rhs )
}

parsePhmmerTable <- function( jr,
  skip.lines=3, parse.line=list( '1'='hit.name', '3'='query.name',
    '6'='bit.score' )
  ) {
  # Parses the table output of HMMER-3's PHMMER. Reads out hit accessions
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

parseBlastTable <- function( br,
  parse.line=list( '1'='query.name', '2'='hit.name',
    '12'='bit.score' )
  ) {
  # Parses a tabular Blast output ( option -m 8 ) reading Query and Hit
  # accessions, as well as the bit scores into a data frame.
  #
  # Args:
  #  br : The result of calling read.table on the path to a tabular Blast
  #               result file.
  #  parse.line : Read out position 'i' and return as list entry, named as
  #               value of 'i'. I.e. position '1' is read out as 'hit.name'. 
  #
  # Returns: Data frame with colnames as the values in argument 'parse.line'
  # and one row per line in br.
  #   
  bt <- br[ , as.integer( names( parse.line ) ) ]
  colnames( bt ) <- as.character( parse.line )
  bt
}

bestHits <- function( seq.search.reslt.mtrx, query.acc, n.best.hits=1000,
  query.acc.column='query.name', sort.column.name='bit.score' ) {
  # Filters the sequence similarity search result table for the n.best.hits
  # using column 'sort.column.name' to sort the results.
  #
  # Args:
  #  seq.search.reslt.mtrx : The matrix holding the sequence similarity search
  #                          results. 
  #  query.acc             : The accession of the query protein to obtain hits
  #                          for.
  #  n.best.hits           : The number of best hits to select, default is 1000
  #  sort.column.name      : The name or index of the column to use for sorting
  #                          the query's hits by.
  #
  # Returns: A matrix with the same columns as seq.search.reslt.mtrx and the
  # subset of maximum n.best.hits rows which hold the sequence similarity
  # search hits for query.acc. The returned Hits are the all or the n.best.hits
  # best ones according to column sort.column.name.
  #   
  query.rslts <- seq.search.reslt.mtrx[
    which( seq.search.reslt.mtrx[ , query.acc.column ] == query.acc ), , drop=F
  ]
  if ( nrow( query.rslts ) > n.best.hits ) {
    query.rslts[
      order( as.numeric( query.rslts[ , sort.column.name ] ), decreasing=T ),
      , drop=F
    ][ 1:n.best.hits, ]
  } else {
    query.rslts
  }
}

sanitizeUniprotAccessions <- function( seq.search.rslt.tbl,
  col.2.san='hit.name', sanitize.funk=sanitizeUniprotAccession ) {
  # Sanitizes each protein accession in column col.2.san of argument
  # seq.search.rslt.tbl.
  #
  # Args:
  #  seq.search.rslt.tbl : The table of sequence similarity search results as
  #                        i.e. obtained by functions parsePhmmerTable(…) or
  #                        parseBlastTable(…)
  #  col.2.san           : The column name or index of seq.search.rslt.tbl's
  #                        column holding the protein accessions to be
  #                        sanitized.
  #  sanitize.funk       : The sanitizing function applied to each single entry
  #                        in col.2.san.
  #
  # Returns: A copy of seq.search.rslt.tbl in which each entry of col.2.san has
  # been replaced with the result of calling sanitize.funk with it.
  #   
  rslt <- seq.search.rslt.tbl
  rslt[ , col.2.san ] <- as.character(
    lapply( as.character( rslt[ , col.2.san ] ), sanitize.funk )
  )
  rslt
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
  hmlgs <- readAAStringSet( path.2.homlgs.fasta )
  dplcts <- duplicated( names( hmlgs ) )
  if ( any( dplcts ) ) {
    uniq.hmlgs <- hmlgs[ - which( dplcts ) ]
    writeXStringSet( uniq.hmlgs, path.2.unique.hmlgs.fasta )
    if ( print.warning ) {
      warning( "Removed duplicated Accessions:",
        paste( names( hmlgs[ dplcts ] ), collapse=", " )
      )
    }
  }
  TRUE
}

filterMultipleSequenceAlignment <- function( msa.xstring.set, min.chars=5 ) {
  # Discards all Amino Acid sequences shorter than 'min.chars' Amino Acids,
  # that is NOT counting gaps.
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

  # In a MSA all sequences ( class "AAString" ) have the same length:
  if ( length( msa.xstring.set[[ 1 ]] ) < min.chars ) {
    NULL # return
  } else {
    msa.aa.lens <- lapply( names( msa.xstring.set ), function( accsn ) {
      nchar( gsub( '-|\\s', '', toString( msa.xstring.set[[ accsn ]] ) ) )
    })
    msa.xstring.set[ which( msa.aa.lens[] > min.chars ) ]
  }
}

chooseFilteredAlignment <- function( msa.unfiltered,
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
  # Returns: TRUE if the filtered Alignment meets the requirements and thus
  # should be used, FALSE otherwise.
  #   
  if ( is.null( msa.filtered ) ||
    length( msa.filtered ) < min.no.seqs || 
    length( msa.filtered ) / length( msa.unfiltered ) < min.fraction
  ) {
    FALSE
  } else {
    TRUE
  }
}

msaEqual <- function( msa.one, msa.two ) {
  # Compares two multiple sequence alignments 'msa.one' and 'msa.two' for
  # identity.
  #
  # Args:
  #  msa.one : The first MSA 
  #  msa.two : The second MSA
  #
  # Returns: TRUE if both MSAs have the same number of sequences and if each
  # position has the same sequences.
  #   
  if ( length( msa.one ) != length( msa.two ) ) {
    FALSE
  } else {
    all( as.logical( lapply( 1:length( msa.one ), function( i ) {
      toString( msa.one[ i ] ) == toString( msa.two[ i ] )
    }) ) )
  }
}

homologsStats <- function( homologs.table, query.protein.accession ) {
  # Generates some statistics about the distribution of bit scores in the
  # current query protein's sequence similarity search results. This function
  # is used to generate PhyloFun's output.
  #
  # Args:
  #  homologs.table : The table of sequence homologs as returned by either
  #                   calling parsePhmmerTable(…) or parseBlastTable(…)
  #  query.protein.accession : The query protein's accession to be used as row
  #                   name of the resulting statistics.
  #
  # Returns: A matrix with the query.protein.accession as row name and the
  # number of sequence homologs in homologs.table and the distribution of bit
  # scores as returned by function summary.
  #   
  m <- t( as.matrix( summary( as.numeric( homologs.table[ , 'bit.score' ] ) ) ) )
  colnames( m ) <- paste( colnames( m ), 'bit.score', sep='-' )
  rownames( m ) <- query.protein.accession
  cbind( matrix( nrow( homologs.table ),
    dimnames=list( query.protein.accession, 'Homologs' ) ),
    m 
  )
}

msaStats <- function( orig.msa, filtered.msa, query.protein.accession ) {
  # Reports statistics on the multiple sequence alignment generated from the
  # query protein and its homologs. Differences in number of retained sequences
  # and positions between the original MSA and the filtered are reported.
  #
  # Args:
  #  orig.msa     : The original MSA as returned by Biostrings' function
  #                 readAAStringSet(…) 
  #  filtered.msa : The filtered MSA as returned by Biostrings' function
  #                 readAAStringSet(…) 
  #  query.protein.accession : The query protein's accession to be used as row
  #                 name of the returned report matrix.
  #
  # Returns: A single row matrix with query.protein.accession as the row name
  # and one column each for the number of sequences and positions in the
  # respective MSAs.
  #   
  npos <- function( msa ) {
    if ( ! is.null( msa ) ) {
      nchar( gsub( '\\s', '', msa[[ 1 ]] ) )
    } else {
      0
    }
  } 
  oln  <- length( orig.msa )
  fln  <- length( filtered.msa )
  opos <- npos( orig.msa )
  fpos <- npos( filtered.msa )
  matrix( c( oln, fln, opos, fpos ), nrow=1,
    dimnames=list( query.protein.accession,
      c( 'Orig.MSA.N.Seqs', 'Filtered.MSA.N.Seqs',
        'Orig.MSA.N.Pos', 'Filtered.MSA.N.Pos' )
    )
  )
}

mergeQueryPredictionsAndHomologAnnotations <- function( query.accession,
  query.predictions, homolog.go.type.annos,
  go.types=c( 'biological_process', 'cellular_component', 'molecular_function')
  ) {
  # Merges the Gene Ontology (GO) term predictions made by PhyloFun for the
  # query protein 'query.accession' with the GO term annotations available for
  # the homologs of the query. 
  #
  # Args:
  #  query.accession       : The accession of the query protein.
  #  query.predictions     : The PhyloFun GO term predictions for the query as
  #                          returned by function goTermPredictionTable(…)
  #  homolog.go.type.annos : The GO term annotations available for the query's
  #                          homologous proteins, i.e. as returned by calling
  #                          goTypeAnnotationMatrices(
  #                            retrieveUniprotAnnotations( homologs.accessions )
  #                          )
  #  go.types              : The names used in returned list, the three GO
  #                          types: BP, CC, and MF.
  #
  # Returns: A named list with GO term annotation matrices. One for each type
  # of GO term. Each matrix has a single row 'GO' and one column for each
  # protein that has GO term anntotations for the respective GO type.
  #   
  if ( is.null( query.accession ) ||
     is.null( query.predictions ) || is.null( homolog.go.type.annos ) ) {
    warning( paste( 'At least one argument for function", 
      "mergeQueryPredictionsAndHomologAnnotations(…) is NULL.", 
      "Returning NULL.' )
    )
    NULL
  } else {
    setNames(
      lapply( go.types, function( go.type ) {
        if ( go.type %in% query.predictions[ , 'term_type' ] ) {
          preds <- sort( as.character(
            query.predictions[
              which( query.predictions[ , 'term_type' ] == go.type ), , drop=F
            ][ , 'acc' ]
          ) )
          m <- matrix( list(), ncol=1, nrow=1,
            dimnames=list( 'GO', query.accession ) )
          m[[ 1, 1 ]] <- preds
          cbind( m, homolog.go.type.annos[[ go.type ]] )
        } else {
          homolog.go.type.annos[[ go.type ]]
        }
      }),
      go.types
    )
  }
}

replaceSelenocysteinInFasta <- function( source.fasta.file,
  filtered.fasta.file=source.fasta.file ) {
  # Reads in a FASTA file of amino acid sequences and filters each sequence
  # with the function replaceSelenocystein( … ). The filtered AA-Sequence
  # set is then stored in 'filtered.fasta.file'.
  #
  # Args:
  #  source.fasta.file : Valid path to the FASTA file.
  #  filtered.fasta.file : Valid path to the FASTA file the filtered
  #                        AA-Sequences shall be saved in.
  #
  # Returns: TRUE if and only if no error occurred.
  #   
  aa.seqs <- readAAStringSet( source.fasta.file )
  aa.seqs.fltrd <- setNames(
    AAStringSet( as.character(
      lapply( names( aa.seqs ), function( accsn ) {
        replaceSelenocystein( aa.seqs[[ accsn ]] )
      })
    ) ),
    names( aa.seqs )
  )
  writeXStringSet( aa.seqs.fltrd, filtered.fasta.file )
  TRUE
}

project.file.path <- function( ..., dir.sep="/" ) {
  # Returns the full file path to a file in subdirectories given in arguments
  # '…'
  #
  # Args:
  #  ...     : The path of subdirectories and finally the file to generate the
  #            complete file path for.
  #  dir.sep : The character to divide dirs with, '/' by default.
  #
  # Returns: The full path
  #   
  paste( path.package( "PhyloFun" ), ..., sep=dir.sep )
}

joinGOTermMutationProbabilityTables <- function( binary.tbl.paths ) {
  # Loads all binary RData in files listed in argument 'binary.tbl.paths' and
  # appends the loaded 'pmds.no.null' lists into a single list
  # GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE.
  #
  # Args:
  #  binary.tbl.paths : A character vector of valid paths to lists of GO term
  #                     mutation probability tables stored as binary RData.
  #                     Each such list should be named 'pmds.no.null'. This
  #                     argument can be created from all files resulting from
  #                     parallel calibration of the mentioned GO term mutation
  #                     probabilities - initialization of this argument can be
  #                     achieved by the following R expression using the
  #                     directory 'res.dir' where all resulting binary lists of
  #                     mutation tables have been stored:
  #                     binary.tbl.paths <-
  #                     system( "find res.dir -name '*.RData'", intern=TRUE )
  #
  # Returns: The named list 'GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE'
  # holding GO term mutation tables for each GO term available.
  # After having called this function the calibration of PhyloFun has been
  # successfully done and the resulting GO term mutation probability
  # distributions should be saved in the PhyloFun package itself, using:
  # save(
  #   GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE,
  #   file="path/2/PhyloFun/data/go_term_mutation_prob_distribs.RData"
  # )
  #   
  GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE <- list()
  lapply( binary.tbl.paths, function( pth ) {
    load( pth )
    GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE <<- append(
      GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE,
      pmds.no.null
    )
  } )
  GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE
}

retrieveGOAnnotations <- function( prot.accs, evidence.codes=EVIDENCE.CODES,
  go.con=connectToGeneOntology(), close.db.con=TRUE ) {
  # Retrieves and merges all GO annotations available for argument protein
  # accessions 'prot.accs'. Data sources are the UniprotKB web services and the
  # Gene Ontology (GO) database as available using argument db connection
  # 'go.con'.
  #
  # Args:
  #  prot.accs      : A character vector of valid UniprotKB protein accessions.
  #                   Usage of sanitizeUniprotAccessions(prot.accs) is strongly
  #                   recommended.
  #  evidence.codes : A character vector of evidence codes to accept as valid
  #                   (trustworthy) to be included in the resulting
  #                   annotations.
  #  go.con         : A valid and active MySQL database connection to an
  #                   instance of the GO database - default is retrieved by
  #                   calling
  #                   connectToGeneOntology(…).
  #  close.db.con   : If set to TRUE the database connection 'go.con' will
  #                   automatically be closed before this function returns.
  #
  # Returns: A data.frame with three columns where each row holds a single
  # Protein's GO annotation. Column one holds the GO term accessions, column
  # two the Evidence Codes, and column three the protein accessions.
  #    
  go.db.annos <- goTermsForProteinAccessionAndEvidenceCodes( prot.accs,
    evidence.codes, go.con ) 
  if ( close.db.con ) {
    dbDisconnect( go.con )
  }
  unipr.go.annos <- if ( is.null( evidence.codes ) || is.na( evidence.codes )
    || evidence.codes == 'ALL' || evidence.codes == 'all' ) {
    retrieveUniprotAnnotations( prot.accs )
  } else {
    retrieveExperimentallyVerifiedGOAnnotations( prot.accs,
      evidence.codes=evidence.codes )
  }
  if ( ! is.null( go.db.annos ) && ! is.na( go.db.annos ) &&
    nrow( go.db.annos) > 0 ) {
    rbind(
      go.db.annos[ , c( 'acc', 'code', 'xref_key' ) ],
    )
  }
}
