precision <- function( predicted.gos, true.gos,
  go.con=connectToGeneOntology() ) {
  # 'precision' is a statistical quality measure of predicted Gene Ontology
  # (GO) terms:
  # precision( predicted.gos ) = | true_positives | / ( | true_positives | + | false_negatives | )
  # Predicted GO terms that are parents of reference terms are NOT counted as
  # false negatives.
  #
  # Args:
  #  predicted.gos : The set of predicted GO terms
  #  true.gos      : The set of experimentally verified GO terms.
  #  go.con        : A database connection to Gene Ontology as returned by
  #                  function connectToGeneOntology().
  #
  # Returns: The precision of the prediction as a numeric between Zero and One.
  # Precision for empty true.gos is always One.
  #   
  pgs <- unique( predicted.gos )
  if ( length( pgs ) == 0 )
    1
  else {
    tgs <- unique( true.gos )
    tp <- intersect( tgs, pgs )
    fps <- falsePositives( pgs, tgs, go.con=go.con )
    denom <- length( tp ) + length( fps )
    if ( denom == 0 )
      1
    else
      length( tp ) / ( denom )
  }
}

falsePositives <- function( predicted.gos, true.gos, 
  go.con=connectToGeneOntology() ) {  
  # Identifies those GO terms in 'predicted.gos' that are not equal to any term
  # in 'true.gos' nor are any of their parent terms.
  #
  # Args:
  #  predicted.gos : The GO terms that have been predicted
  #  true.gos      : The reference GO terms
  #  go.con        : A database connection to the Gene Ontology as returned by
  #                  function connectToGeneOntology().
  #
  # Returns: Those predicted GO terms that are false positives.
  #   
  true.gos.tbl <- goTermsForAccessionWithLevel( true.gos, con=go.con )
  all.parents <- do.call( 'rbind', lapply( true.gos.tbl$id, parentGoTerms,
    con=go.con ) )
  fn.candidates <- setdiff( predicted.gos, true.gos )
  setdiff( fn.candidates, all.parents$acc )
}

recall <- function( predicted.gos, true.gos ) {
  # 'recall' is a statistical quality measure of predicted Gene Ontology (GO)
  # terms:
  # recall( predicted.gos ) = | true_positives | / | true.gos |
  #
  # Args:
  #  predicted.gos : The set of predicted GO terms
  #  true.gos      : The set of experimentally verified GO terms.
  #
  # Returns: The recall of the prediction as a numeric between Zero and One. 
  #
  pgs <- unique( predicted.gos )
  tgs <- unique( true.gos )
  tp <- intersect( tgs, pgs )
  if ( length( tgs ) == 0 ) 1 else length( tp ) / length( tgs )
}

falsePositiveRate <- function( predicted.gos, true.gos ) {
  # The false positive rate is defined as:
  # fpr( predicted.gos ) = | predicted.gos \ true.gos | / | predicted.gos |
  #
  # Args:
  #   : 
  #
  # Returns:
  #   
  pgs <- unique( predicted.gos )
  tgs <- unique( true.gos )
  length( setdiff( pgs, tgs ) ) / length( pgs )
}

fScore <- function( predicted.gos, true.gos, beta.param=1,
  go.con=connectToGeneOntology() ) {
  # Computes the weighted harmonic mean of precision and recall as a quality
  # measure of predicted.gos. The higher beta.param the more emphasis is put on
  # recall than on precision.
  #
  # Args:
  #  predicted.gos : The set of predicted GO terms
  #  true.gos      : The set of experimentally verified GO terms.
  #  beta.param    : The factor of how much more emphasis should be put on
  #                  recall than on precision.
  #  go.con        : database connection to the Gene Ontology as returned by
  #                  function connectToGeneOntology()
  #
  # Returns: The weighted harmonic mean of precision and recall as a numeric
  # between Zero and One.
  #   
  prcsn <- precision( predicted.gos, true.gos, go.con=go.con )
  rcll <- recall( predicted.gos, true.gos )
  bp <- beta.param^2
  if ( 0 == (prcsn + rcll) )
    0
  else
    ( 1 + bp ) * ( prcsn * rcll ) / ( bp * prcsn + rcll )
}

parseBlast2GOresults <- function( b2g.res ) {
  # Parses the output of a Blast2GO annot file and returns an annotation matrix
  # as it is used throughout PhyloFun.
  #
  # Args:
  #  b2g.res : Result of calling readLines on a Blast2GO annot file.
  #
  # Returns: An annotation matrix in which the columns are the GO annotated
  # protein accessions and the single row 'GO' holds in each of its cells the
  # character vector o GO terms each protein has been annotated with.
  #   
  b2g.sanitized <- as.character(
    lapply(  b2g.res, function( l ) {
      str_match( l, '^[^\\t]+\\tGO:\\d{7}' )[[ 1, 1 ]]
    })
  )
  b2g.tbl <- read.table( text=b2g.sanitized[ ! is.na( b2g.sanitized[] ) ] )
  do.call( 'cbind', lapply( unique( b2g.tbl$V1 ), function( acc ) {
      mtrx <- matrix( list(), ncol=1, nrow=1, dimnames=list( 'GO', acc ) )
      mtrx[[ 1, 1 ]]  <- as.character( b2g.tbl[ b2g.tbl$V1 == acc, 2 ] )
      mtrx
    })
  )
}

parseInterProScan2GOresults <- function( ipr.scn.res ) {
  # Parses a typical InterProScan result file for GO annotations. Uses function
  # parseInterProScanTable with a different ipr.regex argument. 
  #
  # Args:
  #  ipr.scn.lines : Lines of the InterProScan result file to parse as result of 'readLines'.
  #
  # Returns:  An annotation matrix in which the columns are the GO annotated
  # protein accessions and the single row 'GO' holds in each of its cells the
  # character vector o GO terms each protein has been annotated with.
  #   
  ipr.go.annos <- parseInterProScanTable( ipr.scn.res, ipr.regex='.*(GO:\\d{7}).*' )
  rownames( ipr.go.annos ) <- 'GO'
  ipr.go.annos
}

rates <- function( protein.accessions, predicted.annotations,
    rate.funk, annotation.type='GO',
    reference.annotations=retrieveExperimentallyVerifiedGOAnnotations(
      protein.accessions
    )
  ) {
  # Computes for each predicted annotation the rates of interest, specified by
  # function argument 'rate.funk'.
  #
  # Args:
  #  protein.accessions    : The query proteins' accessions.
  #  predicted.annotations : The annotations predicted for the query proteins
  #                          'protein.accessions'.
  #  rate.funk             : The function to call to compute the rate of
  #                          interest. This function MUST accept two arguments:
  #                          'pa' a vector of predicted annotations and 'oa' a
  #                          vector of original reference annotations.
  #  annotation.type       : The row name of argument annotation.type to
  #                          evaluate the FPRs for.
  #  reference.annotations : The TRUE annotations to be used as reference.
  #
  # Returns: A named list of recall rates for each predicted annotation. Names
  # are the query protein accessions. Every time no reference annotations are
  # found the returned rate is 'NA'.
  #   
  setNames(
    lapply( protein.accessions, function( a ) {
      if ( a %in% colnames( reference.annotations ) ) {
        # predicted annos for 'a'
        pa <- if ( a %in% colnames( predicted.annotations ) )
         predicted.annotations[[ annotation.type, a ]]
        else 
         c()
        # experimentally verified annos for 'a'
        oa <- reference.annotations[[ annotation.type, a ]]
        # fScore for predictions on 'a'
        if ( ! is.null( oa ) && length( oa ) > 0 )
          rate.funk( pa, oa )
        else
          NA
      } else {
        NA
      }
    }),
    protein.accessions
  )
}

fScores <- function( protein.accessions, predicted.annotations,
  annotation.type='GO', beta.param=1,
  reference.annotations=retrieveExperimentallyVerifiedGOAnnotations(
    protein.accessions ), go.con=connectToGeneOntology(), close.go.con=TRUE ) {
  # Computes the fScores for all predicted annotations.
  #
  # Args:
  #  protein.accessions    : The reference proteins to compute the fScores for
  #  predicted.annotations : The predictions, if any, for the reference
  #                          proteins. Note, that some might be missing.
  #                          Argument should be an annotation matrix with the
  #                          protein accessions as columns and cells the
  #                          predicted annotations.
  #  annotation.type       : The type of annotation predicted, default 'GO'.
  #  beta.param            : Passed to function fScore, see its docu for more
  #                          details.
  #  reference.annotations : The annotations of the reference proteins, default
  #                          experimentally verified GO annotations as
  #                          retrievable from UniProt.
  #  go.con                : Database connection to the Gene Ontology as
  #                          returned by function connectToGeneOntology()
  #  close.go.con          : If set to TRUE the database connection 'go.con'
  #                          will be automatically closed after executing this
  #                          function and before returning its results.
  #
  # Returns: A named list with the predictions' fScores of each reference
  # protein. The fScore NA is assigned whenever no reference annotations are
  # found.
  #   
  rate.funk <- function ( pa, oa ) fScore( pa, oa, beta.param=beta.param,
    go.con=go.con )
  rts <- rates( protein.accessions, predicted.annotations, rate.funk,
    annotation.type, reference.annotations )
  if ( close.go.con ) dbDisconnect( go.con )
  # return
  rts
}

falsePositiveRates <- function( protein.accessions, predicted.annotations,
    annotation.type='GO',
    reference.annotations=retrieveExperimentallyVerifiedGOAnnotations(
      protein.accessions
    )
  ) {
  # Computes for each predicted annotation the false positive rate (FPR).
  #
  # Args:
  #  protein.accessions    : The query proteins' accessions.
  #  predicted.annotations : The annotations predicted for the query proteins
  #                          'protein.accessions'.
  #  annotation.type       : The row name of argument annotation.type to
  #                          evaluate the FPRs for.
  #  reference.annotations : The TRUE annotations to be used as reference.
  #
  # Returns: A named list of false positive rates for each predicted
  # annotation. Names are the query protein accessions.
  #   
  rates( protein.accessions, predicted.annotations, falsePositiveRate,
    annotation.type, reference.annotations )
}

recallRates <- function( protein.accessions, predicted.annotations,
    annotation.type='GO',
    reference.annotations=retrieveExperimentallyVerifiedGOAnnotations(
      protein.accessions
    )
  ) {
  # Computes for each predicted annotation the recall rate.
  #
  # Args:
  #  protein.accessions    : The query proteins' accessions.
  #  predicted.annotations : The annotations predicted for the query proteins
  #                          'protein.accessions'.
  #  annotation.type       : The row name of argument annotation.type to
  #                          evaluate the FPRs for.
  #  reference.annotations : The TRUE annotations to be used as reference.
  #
  # Returns: A named list of recall rates for each predicted
  # annotation. Names are the query protein accessions.
  #   
  rates( protein.accessions, predicted.annotations, recall,
    annotation.type, reference.annotations )
}

