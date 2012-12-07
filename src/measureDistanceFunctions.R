measureDistances <- function( annotation, annotation.matrix,
  blast.result.tbl, aa.seqs, domain.weights.table,
  annotation.type="GO", round.2.digits=2,
  blast.tbl.query.col=1, blast.tbl.hit.col=2,
  lapply.funk=lapply ) {
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
  # pair (accessions in alphabetical order) and columns are sequence and
  # domain architecture distance, respectively. Last column holds a boolean
  # indicating wether the pairs shares the argument annotation.
  #   
  annot.annot.mtrx <- sharedAnnotation( annotation.matrix, annotation,
    annotation.type=annotation.type
  )
  annot.pairs.tbl <- pairsForAccessionsAssumingSymmetry( blast.result.tbl,
    colnames( annot.annot.mtrx)
  )
  # Measure distances and sort by pairs sharing and not sharing the annotation,
  # respectively:
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
      matrix(
        c( seq.dist, das.dist, shrd.annot ), nrow=1,
        dimnames=list(
          sortedPairName( acc.a, acc.b ),
          c( "Sequence.Distance", "Domain.Architecture.Distance",
            paste("Share.", annotation, sep="")
          )
        )
      )
    })
  )
}

mutationProbabilityDistribution <- function( distances.tbl, 
  p.column.name, dep.column.name, 
  p.grid=seq(0.1,1,by=0.1)
  ) {
  p.mut.colname <- paste( "p.mutation|", p.column.name, sep="" )

  srtd <- distances.tbl[ sort.list( distances.tbl[ , p.column.name ] ), ]
  pairs.sharing <- 0; pairs.diff <- 0; row.ind <- 1; p.mut.last <- 0;
  p.mut.mtrx <- do.call('rbind',
    lapply( 1:nrow(srtd), function(i) {
      if ( srtd[[ i, 3 ]] ) {
        pairs.sharing <<- pairs.sharing + 1
      } else {
        pairs.diff <<- pairs.diff + 1
      }
      p.mut <- if( pairs.diff > pairs.sharing ) 1 else pairs.diff / pairs.sharing
      
      # Mutation probability can not decline with increasing distances:
      if ( p.mut < p.mut.last )
        p.mut <- p.mut.last
      p.mut.last <<- p.mut

      setNames(
        list( p.mut, srtd[[ i, p.column.name ]], srtd[[ i, dep.column.name ]], srtd[[ i, 3 ]] ),
        list( p.mut.colname, p.column.name, dep.column.name, colnames( srtd )[3] )
      )
    })
  )
  rownames( p.mut.mtrx ) <- rownames( srtd )
  
  # return
  p.mut.mtrx
}
