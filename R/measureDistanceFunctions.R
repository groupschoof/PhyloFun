measureDistances <- function( annotation, annotation.matrix,
  blast.result.tbl, aa.seqs, domain.weights.table,
  annotation.type="GO", round.2.digits=2,
  blast.tbl.query.col=1, blast.tbl.hit.col=2,
  lapply.funk=lapply ) {
  #
  # DEPRECATED
  #
  # For each protein pair, where at least one partner has the argument
  # 'annotation' the sequence and domain architecture distances are
  # computed. The distances are returned along with a boolean indicating
  # wether the pair shares the annotation or just one member has it. This
  # data can be used to compute the function mutation pobability
  # distribution for any given distance. 
  #
  # Args:
  #  annotation           : The GO or InterPro annotation to find protein
  #                         pairs for.
  #  blast.result.tbl     : Protein pairs of significant sequence
  #                         similarity.
  #  aa.seqs              : The amino acid sequences of the proteins
  #                         referenced in the other arguments.
  #  domain.weights.table : Weights of InterPro domains used to compute the
  #                         domain architecture distances.
  #  annotation.type      : Wether to use InterPro or GO (default).
  #  round.2.digits       : Round the computed distances to these digits
  #                         (default 2).
  #  blast.tbl.query.col  : The column index in which to find the Blast
  #                         query accessions.
  #  blast.tbl.hit.col    : The column index in which to find the Blast hit
  #                         accessions.
  #  lapply.funk          : Function to use when iterating over the protein
  #                         pairs. Set to 'mclapply' if parallel computation
  #                         is wanted.
  #
  # Returns: a matrix where each row is labelled with the considered protein
  # pair and columns are sequence and domain architecture distance,
  # respectively. Another column holds a boolean indicating wether the pairs
  # shares the argument annotation. The last column hold the pairwise euclidean
  # distance of the distance pair coordinates to (0,0). The last column is used
  # to sort the measured distances when computing the mutation probability
  # depending on _both_ distances: p( "Mutation" | Sequence.Distance,
  # Domain.Architecture.Distance )
  #   
  annot.annot.mtrx <- sharedAnnotation( annotation.matrix, annotation,
    annotation.type=annotation.type
  )
  annot.pairs.tbl <- pairsForAccessionsAssumingSymmetry( blast.result.tbl,
    colnames( annot.annot.mtrx)
  )
  # Measure distances:
  if( ! is.null( annot.pairs.tbl ) && nrow( annot.pairs.tbl ) > 0 ) {
    do.call('rbind', 
      lapply.funk( 1:nrow(annot.pairs.tbl), function(i) {
        acc.a <- as.character( annot.pairs.tbl[[ i, blast.tbl.query.col ]] ) 
        acc.b <- as.character( annot.pairs.tbl[[ i, blast.tbl.hit.col ]] ) 
        seq.dist <- round(
          pairwiseSequenceDistance( as.character(aa.seqs[ acc.a ]),
            as.character(aa.seqs[ acc.b ])
          ), 
          round.2.digits
        )
        das.dist <- round(
          partialDomainArchitectureDistances( annotation.matrix,
            domain.weights.table, list( c(acc.a, acc.b) )
          )[[1]],
          round.2.digits
        )
        shrd.annot <- shareAnnotation( annotation, annotation.matrix,
          acc.a, acc.b
        )
        eucl.dist <- round( 
          dist( matrix( c(0, seq.dist, 0, das.dist), nrow=2 ) )[[1]],
          round.2.digits
        )
        matrix(
          c( seq.dist, das.dist, shrd.annot, eucl.dist ), nrow=1,
          dimnames=list(
            paste( acc.a, acc.b, sep="_" ),
            c( "Sequence.Distance", "Domain.Architecture.Distance",
              paste("Share.", annotation, sep=""),
              "Euclidean.Distance.To.Origin"
            )
          )
        )
      })
    )
  } else {
    write( paste( "No protein pairs of significant sequence similarity found for", annotation ), stderr() )
    NULL
  }
}

pMutation <- function( pairs.shared, pairs.diff, min.p.mut=0 ) {
  # Computes the fraction pairs.diff / ( pairs.shared + pairs.diff ).
  #
  # Args:
  #  pairs.shared : The number of protein pairs sharing a certain annotation. 
  #  pairs.diff : The number of protein pairs NOT sharing a certain annotation. 
  #  min.p.mut : The minimum of the computed p.mutation to return.
  #
  # Returns: The computed fraction or zero, if denominator is zero. If computed
  # fraction is smaller than argument 'min.p.mut' min.p.mut is returned.
  #   
  denominator <- pairs.shared + pairs.diff
  p.mut <- if ( denominator == 0 ) 0 else pairs.diff / denominator
  if ( p.mut < min.p.mut ) min.p.mut else p.mut
}

mutationProbabilityDistribution <- function( distances.tbl, annotation.tbl,
  annotation, pairs.first.col=1, pairs.secnd.col=2, distance.col=3,
  annot.col=1, prot.col=3, p.mut.colname='p.mutation|Sequence.Distance',
  distance.colname = 'Sequence.Distance' ) {
  # Measure function mutation probability distribution for argument
  # 'annotation' depending on the distances given for protein pairs where at
  # least one member has the argument 'annotation'. These pairs are sorted by
  # their ascending distances, and the annotation mutation probability p is
  # iteratively calculated given how many pairs sharing and not sharing the
  # current function annotation have been found that have a distance <= to the
  # current iteration's. This computation only returns increasing values for p,
  # if more pairs sharing the annotation are found in later iterations and thus
  # p would decrease for the current distance value, the last and higher value
  # is used as p for the current pair again.
  # 
  # Args:
  #  distances.tbl     : A three column table in which the first two columns hold
  #                      the respective accessions of proteins forming a pair,
  #                      and a third column with the precomputed distance
  #                      value.
  #  annotation.tbl    : A table holding the annotations for each protein given in
  #                      argument 'distances.tbl'. Required format is one row
  #                      per annotation.
  #  annotation        : The annotation to compute its mutation probability
  #                      distribution for, e.g. 'GO:0001234'.
  #  pairs.first.col   : The column of 'distances.tbl' in which to lookup
  #                      the first member of protein pairs.
  #  pairs.secnd.col   : The column of 'distances.tbl' in which to lookup
  #                      the second member of protein pairs.
  #  distance.col      : The column of 'distances.tbl' in which to lookup the
  #                      precomputed distance.
  #  prot.col          : The column of 'annotation.tbl' in which to lookup the
  #                      protein accessions.
  #  annot.col         : The column of 'annotation.tbl' in which to lookup the
  #                      annotation.
  #  p.mut.colname     : The name of the column to hold the mutation
  #                      probabilities. Default is
  #                      'p.mutation|Sequence.Distance'.
  #  distance.colname  : The name of the column to hold the distance measure
  #                      upon which the mutation probability is computed.
  #                      Default is 'Sequence.Distance'.
  #  
  #
  # Returns: A sorted subset of argument matrix 'distances.tbl' with an
  # additional column holding the measured annotation mutation probabilities
  # for each row. Sorting is done on each pair's distance.
  #   
  prot.pairs.wt.anno <- proteinPairsSharingAnnotation( annotation,
    distances.tbl, annotation.tbl, pairs.first.col, pairs.secnd.col, prot.col,
    annot.col
  )
  if ( nrow( prot.pairs.wt.anno ) > 0 ) {
    srtd <- prot.pairs.wt.anno[
      sort.list( as.numeric( prot.pairs.wt.anno[ , distance.col ] ) ), ,
      drop=FALSE
    ]
    pairs.sharing <- 0; pairs.diff <- 0; p.mut.last <- 0;
    dists <- unique( as.numeric( srtd[ , distance.col ] ) )
    matrix( byrow=TRUE, ncol=2, 
      dimnames=list( c(), c( p.mut.colname, distance.colname ) ),
      unlist( lapply( dists, function(dist) {
        # Count pairs sharing and not sharing annotation for current value:
        candidates <- if( is.na(dist) )
          srtd[ which( is.na( srtd[ , distance.col ] ) ), , drop=F ]
        else
          srtd[ which( srtd[ , distance.col ] == dist ), , drop=F ]

        no.cand.sharing.anno <- nrow(
          candidates[ which( candidates[ , 'annotation.shared' ] == TRUE ), ,
            drop=F ]
        )
        pairs.sharing <<- pairs.sharing + no.cand.sharing.anno
        pairs.diff <<- pairs.diff + ( nrow(candidates) - no.cand.sharing.anno )

        p.mut <- pMutation( pairs.sharing, pairs.diff, p.mut.last )
        # Mutation probability must not decline with increasing distances:
        p.mut.last <<- p.mut

        # Current row:
        c( p.mut, dist )
      } ) ) 
    )
  } else {
    NULL
  }
}

measureEuclideanDists <- function( dists.mtrx, dists.to.coordinates=c( 0,0 ),
  dists.mtrx.columns=c( 1,2 ), lapply.funk=lapply ) {
  # Measures pairwise euclidean distances to argument 'dists.to.coordinates'
  # for each row of argument 'dists.mtrx', coordinates are extracted using
  # indices defined in argument 'dists.mtrx.columns'. 
  #
  # Args:
  #  dists.mtrx : Matrix, where each row is a set of euclidean coordinates.
  #               I.e. as returned by function 'measureDistances'.
  #  dists.to.coordinates : Measure distances to this point in vector space.
  #  dists.mtrx.columns : The column indices of the argument 'dists.mtrx' to
  #                       measure euclidean distances with.
  #  lapply.funk : Set to 'mclapply' if computation is wanted in parallel.
  #
  # Returns: An extended copy of the argument matrix 'dists.mtrx' with an
  # additional column holding the pairwise euclidean distances to argument
  # 'dists.to.coordinates'. 
  #   
  paste.funk <- function( ... ) { paste( ..., sep="." ) }
  euc.dist.col <- do.call( 'paste.funk',
    as.list( c( "Euclidean.Distance.To", dists.to.coordinates ) )
  )
  
  do.call( 'rbind',
    lapply.funk( 1:nrow(dists.mtrx), function(i) {
      matrix( c( dists.mtrx[ i, ],
          dist(
            matrix( 
              c(
                dists.mtrx[ i, dists.mtrx.columns ],
                dists.to.coordinates
              ), nrow=2, byrow=T
            )
          )[[1]]
        ),
        nrow=1,
        dimnames=list(
          rownames(dists.mtrx)[[i]],
          c( colnames(dists.mtrx), euc.dist.col )
        ) 
      )
    })
  )
}

gridPMutation <- function( dists.mtrx, grid.by=0.1,
  round.2.digits=( nchar( as.character(grid.by) ) - 2 ),
  p.column.index=1, lapply.funk=lapply ) {
  # Fits the measured probabilites in argument column 'p.column.index' of
  # argument matrix 'dists.mtrx' to the grid sequence from argument 'grid.by'
  # to one, by 'grid.by' steps. Fitting is done by rounding to nearest grid
  # value. If no suitable value for a given grid.point is encountered in
  # interval ( grid.point - grid.by / 2, grid.point + grid.by / 2 ] NA is
  # assigned to that grid.point. 
  #
  # Args:
  #  dists.mtrx : The matrix of mutation probabilities and distances measured
  #               for each pair of proteins, as returned by function
  #               'mutationProbabilityDistribution'.
  #  grid.by    : The steps to create the grid sequence with as in seq(
  #               grid.by, 1, by=grid.by ).
  #  round.2.digits : The grid sequence is rounded to these digits to ensure
  #                   rounding errors in comparisons. Has to have as many
  #                   decimal digits as argument 'grid.by'.
  #  p.column.index : The index of argument 'dists.mtrx' column in which to
  #                   find the measured probabilites.
  #  lapply.funk : Set to 'mclapply' if parallel execution is wanted.
  #
  # Returns: A matrix with as many rows as grid points are created, where each
  # row hold the first matching row of argument 'dists.mtrx' for its
  # corresponding grid.point. The grid.points are stored in the first column of
  # the matrix, while the other columns are unchanged as they appear in
  # 'dists.mtrx'.
  #   
  grd <- round( seq( grid.by, 1, by=grid.by ),
    round.2.digits
  )
  interval <- grid.by
  gridded <- do.call( 'rbind', 
    lapply.funk( grd, function( x ) {
      hits <- dists.mtrx[ 
        -1 * ( dists.mtrx[ , p.column.index ] - x ) <= interval & 
        -1 * ( dists.mtrx[ , p.column.index ] - x ) >= 0, , drop=F 
      ]
      if ( nrow(hits) > 0 ) {
        hits[ nrow(hits), , drop=F ]
      } else {
        NA
      }
    })
  )
  gridded <- cbind( grd, gridded )
  rownames( gridded ) <- c()
  colnames( gridded ) <- c(
    paste( colnames(dists.mtrx)[p.column.index], "grid", sep="." ),
    colnames( dists.mtrx )
  )
  gridded
}

pMutationMinMaxParentValues <- function( p.mut.dists.mtrx, p.column,
  parent.columns=c("Sequence.Distance", "Domain.Architecture.Distance",
    "Euclidean.Distance.To.Origin"), round.2.digits=2 ) {
  # Rounds the measured probabilities to 'round.2.digits' and generates a
  # unique sorted list of them. For each of these p-values the matching rows in
  # argument 'p.mut.dists.mtrx' are selected and the minimum and maximum values
  # of columns 'parent.columns' extracted. These values are composed into a new
  # matrix.
  #
  # Args:
  #  p.mut.dists.mtrx : The matrix of mutation probabilites and measured
  #                     distances as returned by function
  #                     'mutationProbabilityDistribution'.
  #  p.column : The name of the column in which the probabilities are stored.
  #  parent.columns : Vector of column names the mutation probability depend
  #                   on, the min and max values are computed for these
  #                   columns.
  #  round.2.digits : The number of decimal digits to round the p-values to.
  #
  # Returns: A matrix in which for each unique rounded p-value the found
  # matching min and max values of the 'parent.columns' are stored.
  #   
  rnd <- p.mut.dists.mtrx
  rnd[ , p.column ] <- round( rnd[ , p.column ], round.2.digits )
  ps <- unique( rnd[ , p.column ] )
  do.call( 'rbind',
    lapply( ps, function( p ) {
      candidates <- rnd[ rnd[ , p.column ] == p, , drop=F ] 
      parent.mins <- lapply( parent.columns, function( p.parent ) {
        min( candidates[ , p.parent ], na.rm=T )
      })
      parent.maxs <- lapply( parent.columns, function( p.parent ) {
        max( candidates[ , p.parent ], na.rm=T )
      })
      matrix( c( p, unlist(parent.mins), unlist( parent.maxs ) ), nrow=1,
        dimnames=list( c(), c( p.column, paste( "min", parent.columns, sep="." ),
            paste( "max", parent.columns, sep="." ) )
        )
      )
    })
  )
}
