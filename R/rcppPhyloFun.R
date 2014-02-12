findMatchingRow <- function( sTable, val, colInd ) {
  # Using the argument matrix 'sTable's column identified by its index 'colInd'
  # this function finds the index of the columns value >= than argument 'val'.
  # The row of this index is returned. If none can be found the last row is
  # returned.
  #
  # Args:
  #  sTable : Numeric matrix to find matching row in.
  #  val    : Value to match with cells in 'sTable's column 'colInd'.
  #  colInd : The index of 'sTable's column whose cells to match val with.
  #           Remember that colInd has to be counted from column 0 not 1 as
  #           normal in R. This is because colInd is passed into C++ which
  #           indexes from 0, not from 1, as does R.
  #
  # Returns: Numeric vector the matched row of matrix sTable. If none matches
  # returns the last row.
  #   
  ret <- .Call( "findMatchingRow", as.matrix( sTable ), val, colInd )
  return( ret )
}

mutationProbability <- function( compositeAnnotation, branchLength,
  annotsMutationProbTables, distanceColumnIndx ) {
  # Calls Rcpp function of same name to infer the maximum mutation probability
  # for any annotation present in 'compositeAnnotation'. Note that the column
  # in which to lookup the distance measure 'branchLength' in the respective
  # mutation probability tables, held in 'annotsMutationProbTables', is
  # specified by 'distanceColumnIndx'. This has to be in C++ format, that is
  # starting with 0 not with 1 as in R. For example column 5 should be passed
  # as 4 to be working in C++.
  #
  # Args:
  #  compositeAnnotation      : Character vector of annotations, for example c(
  #                             "GO_1", "GO_2", "GO_3" )
  #  branchLength             : The length of the current branch as numeric.
  #  annotsMutationProbTables : The list of numeric matrices holding the
  #                             mutation probability tables for the annotations
  #                             present in compositeAnnotation. See constant
  #                             GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE
  #                             in this projects data.
  #  distanceColumnIndx       : The index of the column of above mutation
  #                             probability tables in which to lookup the
  #                             distance measure matching 'branchLength'.
  #                             Remember that it will be passed into C++ and
  #                             thus has to be counted from column 0, not 1 as
  #                             in R.
  #
  # Returns: The numeric vector of length one holding the looked up mutation
  # probability.
  #   

  # Pass only the required mutation probability tables as matrices:
  annotsMutProbTbls <- filterAnnotationMutationProbabilityTableList( compositeAnnotation,
    annotsMutationProbTables )
  ret <- .Call( "mutationProbability", compositeAnnotation, branchLength,
               annotsMutProbTbls, distanceColumnIndx, PACKAGE="PhyloFun" )
  return( ret )
}

conditionalProbabilityTable <- function( branchLength, annos,
  stringifiedAnnotations, annotsMutationProbTables, mutTblLengthColIndx ) {
  # Generates the numeric probability matrix which holds the mutation
  # probabilities for each composite annotation to each other in argument list
  # 'annos'. Calls the respective Rcpp function of same name.
  #
  # Args:
  #  branchLength             : The length of the current phylogenetic branch
  #                             in expected character change per site.
  #  annos                    : The list of composite annotations present in
  #                             the current phylogenetic tree.
  #  stringifiedAnnotations   : The list 'annos' converted to strings using
  #                             function annotationToString(…)
  #  annotsMutationProbTables : The list of mutation probability tables for
  #                             each atomic annotation term present in the
  #                             composite annotations in 'annos'.
  #  mutTblLengthColIndx      : The index of the column in the mutation
  #                             probability tables in which to lookup the
  #                             distance measure matching 'branchLength'.
  #
  # Returns: A numeric matrix holding the respective mutation probabilities for
  # each composite annotation (CA) to each CA in 'annos'. Column and row names
  # are the respective stringified annotations 'stringifiedAnnotations'.
  #   

  # Pass only the required mutation probability tables as matrices:
  annotsMutProbTbls <- filterAnnotationMutationProbabilityTableList( annos,
    annotsMutationProbTables )
  ret <- .Call( "conditionalProbabilityTable", branchLength, annos,
               stringifiedAnnotations, annotsMutProbTbls, mutTblLengthColIndx,
               PACKAGE="PhyloFun" )
  return( ret )
}

conditionalProbabilityTables <- function( uniqueEdgeLengths, annos,
  stringifiedAnnotations, annotsMutationProbTables, mutTblLengthColIndx ) {
  # Generates the conditional probability tables for each unique branch length
  # in 'uniqueEdgeLengths'. Uses Rcpp function conditionalProbabilityTables.
  # Returned probability matrices hold mutation probabilities from each
  # composite annotation to each composite annotation found in 'annos'.
  #
  # Args:
  #  uniqueEdgeLengths        : The numeric vector of pairwise distinct branch
  #                             lengths as found in the current phylogenetic
  #                             tree.
  #  annos                    : The list of unique composite annotations.
  #  stringifiedAnnotations   : 'annos' stringified as done by function
  #                             annotationToString(…).
  #  annotsMutationProbTables : The list of mutation probability tables for
  #                             each atomic annotation term present in the
  #                             composite annotations in 'annos'.
  #  mutTblLengthColIndx      : The index of the column in the mutation
  #                             probability tables in which to lookup the
  #                             distance measure matching 'branchLength'.
  #
  # Returns: List of conditional probability tables one for each unique branch
  # length in 'uniqueEdgeLengths'. 'uniqueEdgeLengths' converted to strings are
  # the names of the returned list and values are numeric matrices.
  # 

  # Pass only the required mutation probability tables as matrices:
  annotsMutProbTbls <- filterAnnotationMutationProbabilityTableList( annos,
    annotsMutationProbTables )
  ret <- .Call( "conditionalProbabilityTables", uniqueEdgeLengths, annos,
               stringifiedAnnotations, annotsMutProbTbls,
               mutTblLengthColIndx, PACKAGE="PhyloFun" )
  return( ret )
}

filterAnnotationMutationProbabilityTableList <- function( annotations,
  annotsMutationProbTables ) {
  # Filters the large list 'annotsMutationProbTables' for all unique annotation
  # terms found in 'annotations'. The so shortened list is processed to contain
  # numeric matrices and not data.frames. In short: This function shortens the
  # argument list of tables and ensures the contained types are numeric
  # matrices.
  #
  # Args:
  #  annotations :              A list or character vector of arbitrary depth
  #                             holding atomic or composite annotation terms.
  #  annotsMutationProbTables : A list of measured mutation probabilities.
  #                             Names are the atomic annotation terms ( for
  #                             example "GO:012344" ) and the values are the
  #                             tables as data.frames or numeric matrices. See
  #                             generation of these tables in PhyloFun's
  #                             documentation ( README.textile ).
  #
  # Returns: A list of measured mutation probability tables for the atomic
  # annotation terms found in 'annotations'. Names are these atomic annotations
  # and values are numeric matrices as generated by funnction
  # pMutationMinMaxParentValues(…). These precomputed tables are stored in
  # PhyloFun's data list "GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE".
  #   
  uniq.annos <- intersect(
    unique( unlist( annotations ) ), 
    names( annotsMutationProbTables )
  )
  lapply( annotsMutationProbTables[ uniq.annos ], as.matrix )
}

uniqueProteinPairs <- function( proteinPairsTbl,
  pairFirstMemberColIndex=0, pairSecondMemberColIndex=1 ) {
  # Identifies unique pairs in the set of un-ordered protein pairs in
  # 'proteinPairsTbl'. As these pairs are symmetric, pairs (A,B) and (B,A) are
  # identical and only one will be retained in the resulting set. Also self
  # matched, that is identity pairs like (A,A) will be excluded.
  #
  # Args:
  #  proteinPairsTbl          : A two column matrix in which each row
  #                             represents a pair of parotein pair accessions.
  #  pairFirstMemberColIndex  : The index of the column in which to lookup the
  #                             first member of protein pairs.
  #  pairSecondMemberColIndex : The index of the column in which to lookup the
  #                             second member of protein pairs.
  #
  # Returns: A subset of argument proteinPairsTbl holding unique, pairwise
  # distinct, protein pairs.
  #   
  prot.tbl <- if ( class( proteinPairsTbl ) == "data.frame" ) {
    as.matrix( proteinPairsTbl )
  } else {
    proteinPairsTbl
  }

  prs.list <- .Call( "uniqueProteinPairs", prot.tbl, pairFirstMemberColIndex,
    pairSecondMemberColIndex )

  # Convert to a 'character'-data.frame:
  as.data.frame(
    matrix( unlist( prs.list ), ncol=2, nrow=length( prs.list ), byrow=TRUE ),
    stringsAsFactors=FALSE
  )
}
