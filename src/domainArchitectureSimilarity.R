library(Biostrings)
library(phangorn)

constructVectorSpaceModel <- function( annotation.matrix, type='InterPro' ) {
  # Constructs a character vector of alphabetically sorted annotations IDs,
  # excluding NA annotations. Uses function 'uniq.annotations' to process the
  # 'annotation.matrix'.
  #
  # Args:
  #  annotation.matrix : Rows are annotation types and columns the annotated
  #                      proteins, as returned by function
  #                      'retrieve.annotations.parallel'
  #  type : Column of 'annotation.matrix' to be processed.
  #
  # Returns: A character vector of alphabetically sorted pairwise distinct
  # annotation IDs coerced from the 'annotation.matrix' corresponding row as
  # defined by 'type'.
  #   
  uas <- uniq.annotations( annotation.matrix[ type, , drop=F ] )
  uas[ ! is.na( uas[] ) ]
}

generateDomainArchitectureSpaceVectors <- function( vector.space.model,
  annotation.matrix, domain.weights.table, annotation.type='InterPro', 
  vectors.4.accessions=colnames(annotation.matrix),
  lapply.funk=lapply ) {
  # All proteins as in the column names of 'annotation.matrix' are processed.
  # For each of these proteins a list is generated where the list positions
  # hold the domain weights as in table 'domain.weights.table' if the protein
  # is annotated with the corresponding domain or 0.0 otherwise.
  #
  # Args:
  #  vector.space.model : A alphabetically sorted character vector of all
  #                       unique annotations. See function 'uniq.annotations'
  #                       for details.
  #  annotation.matrix  : Columns are the annotated proteins, cells hold each
  #                       protein's annotations.
  #  domain.weights.table : The matrix of domain weights. Rows are expected to
  #                         be the annotated domains and cells expected to hold
  #                         their corresponding domain weights. See project
  #                         'generate_domain_weights' for details.
  #  annotation.type    : The row of the annotation.matrix to access. If NULL
  #                       row 1 is used.
  #  vectors.4.accessions : Compute the DAS vectors for these accessions. This
  #                         enables partial computation of DAS vectors. By default DAS vectors for all
  #                         accessions are generated.
  #  lapply.funk : Set to mclapply, if computation is to be done in parallel.
  #                Obviously default mode is serial ('lapply'). 
  #
  # Returns: List of domain architecture space vectors one for each annotated
  # protein.
  #   
  amr <- if ( is.null(annotation.type) ) 1 else annotation.type

  do.call('cbind',
    setNames(
      lapply.funk( vectors.4.accessions, function( protein.id ) {
        sapply( vector.space.model, function( domain.id ) {
            prot.annos <- annotation.matrix[[ amr, protein.id ]]
            if ( domain.id %in% prot.annos ) {
              dw <- domain.weights.table[ domain.id, 1 ]
              # Domain IDs not present in the domain.weights.table will be
              # ignored:
              if ( is.na(dw) ) 0.0 else dw
            } else {
              0.0
            }
        })
      }),
      vectors.4.accessions )
  )
}

pairwiseDomainArchitectureDistance <- function( dom.space.vector.one,
  dom.space.vector.two ) {
  # Computes 1.0 minus the similarity score of two vectors in the domain
  # architecture space. Where this similarity score is computed as the cosine
  # of the angle between those two vectors. The distance between two such
  # vectors is then returned as 1 - similarity.
  #
  # Args:
  #   dom.space.vector.one, dom.space.vector.two : Vectors in the domain
  #                                                architecture space.
  #
  # Returns: 1.0 minus the cosine of the angle between the two argument
  # vectors. A vector with elements 'NaN' if one or both arguments are vectors
  # of length zero. 
  #   
  a <- as.numeric( dom.space.vector.one )
  b <- as.numeric( dom.space.vector.two )
  if ( identical(a,b) ) {
    # Identical vectors have no distance to themselves:
    0.0
  } else {
    # Vectors are NOT identical:
    1.0 - sum( a * b ) / ( sqrt(sum(a^2)) * sqrt(sum(b^2)) )
  }
}

domainArchitectureDistances <- function( domain.architecture.space.vectors ) {
  # Method computes pairwise distances of the domain architecture space vectors
  # given in argument matrix 'domain.architecture.space.vectors'. Column names
  # are expected to be the Proteins and row names the Domains, hence the ith
  # column is the vector of protein i. 
  #
  # Args:
  #  domain.architecture.space.vectors : The matrix of vectors in domain
  #                                      architecture space.
  #
  # Returns: Returns an object of class dist, holding the pairwise distances of
  # the argument vectors.
  #   
  
  # All proteins to compute pairwise distances for
  ps <- colnames( domain.architecture.space.vectors )
  # Iterativly compute cells of this distance matrix:
  m <- matrix( nrow=length(ps), ncol=length(ps),
    dimnames=list(ps, ps) )
  lapply( ps, function( protein.acc ) {
    # Proteins to compute distances to:
    i <- which( ps == protein.acc )
    ps2c <- ps[ i : length( ps ) ]
    lapply( ps2c, function( p2c ) {
      m[[ p2c, protein.acc ]] <<- pairwiseDomainArchitectureDistance(
        domain.architecture.space.vectors[ , protein.acc ],
        domain.architecture.space.vectors[ , p2c ]
        )
    })
  })
  # Return an object of type 'dist':
  as.dist(m)
}

partialDomainArchitectureDistances <- function( annotation.matrix,
  domain.weights.table, prot.pairs, annotation.type="InterPro",
  lapply.funk=lapply ) {
  # Computes pairwise distances in domain architecture space (DAS) for all
  # protein pairs in argument 'prot.pairs'.
  #
  # Args:
  #  annotation.matrix    : Matrix of domain annotations
  #  domain.weights.table : domain weight (column) for each annotated domain
  #                         (row)
  #  prot.pairs           : Pairs of proteins to compute distances in DAS for.
  #  lapply.funk          : 'lapply' function to use. Set to mclapply, if
  #                         parallel computation is wanted.
  #  annotation.type      : The type of domain annotation to use, default
  #                         'InterPro'
  #
  # Returns: Returns a named list of domain architecture distances for each
  # pair of proteins in argument 'prot.pairs'. List names are generated by
  # function 'pairwiseDistanceKey'. 
  #   
  setNames(
    lapply.funk( prot.pairs, function( prot.pair ) {
      accs <- c( prot.pair[[1]], prot.pair[[2]] )
      vsm <- constructVectorSpaceModel(
        annotation.matrix[ , accs, drop=F ] )
      das.vects <- generateDomainArchitectureSpaceVectors( vsm,
        annotation.matrix, domain.weights.table,
        annotation.type=annotation.type, vectors.4.accessions=accs )
      pairwiseDomainArchitectureDistance(
        das.vects[ , prot.pair[[1]] ],
        das.vects[ , prot.pair[[2]] ]
      )
    }),
    lapply.funk( prot.pairs, function( prot.pair ) {
      pairwiseDistanceKey( prot.pair[[1]], prot.pair[[2]] )
    })
  )
}

pairwiseSequenceDistance <- function( aa.seq.pattern, aa.seq.subject,
  sub.matrix='PAM250', gap.open.pnlty=-10, gap.extension.pnlty=-0.1,
  distance.model='Dayhoff') {
  # For the two argument amino acid sequences this function computes first a
  # global pairwise alignment based on the supplied substitution matrix with
  # the argument gap opening and extension penalties. Then the sequence
  # distance is computed based on the argument model. See libraries Biostrings
  # and phangorn for details on functions pairwiseAlignment and dist.ml.
  #
  # Args:
  #  aa.seq.pattern, aa.seq.subject : argument amino acid sequences
  #  sub.matrix                     : substitution matrix to use in the global
  #                                   alignment
  #  gap.open.pnlty                 : score penalty to apply for opening a gap
  #                                   in the global alignment
  #  gap.extension.pnlty            : score penalty to apply for extending a
  #                                   gap in the global alignment
  #  distance.model                 : model to base the sequence distance
  #                                   computation on
  #
  # Returns: The computed sequence distance for the two argument amino acid
  # sequences as instance of 'dist'. 
  #   
  funk.body <- function( aa.seq.pattern, aa.seq.subject, sub.matrix="PAM250",
    gap.open.pnlty=-10, gap.extension.pnlty=-0.1, distance.model="Dayhoff" ) {
    pa <- pairwiseAlignment(
      AAString( replaceSelenocystein( aa.seq.pattern ) ),
      AAString( replaceSelenocystein( aa.seq.subject ) ),
      substitutionMatrix=sub.matrix, gapOpening=gap.open.pnlty,
      gapExtension=gap.extension.pnlty
    )
    pd <- phyDat(
      list(
        "P1"=strsplit( toString( pattern(pa) ), NULL ),
        "P2"=strsplit( toString( subject(pa) ), NULL )
        ),
      type="AA"
      )
    # Return only the numeric distance value:
    as.matrix( dist.ml( pd, model=distance.model ) )[[ 1,2 ]]
  }
  # Run it safe:
  tryCatch(
    funk.body( aa.seq.pattern, aa.seq.subject, sub.matrix,
      gap.open.pnlty, gap.extension.pnlty, distance.model ),
    error=function( e ) {
      msg <- paste(
        "Could not compute pairwise sequence distance for arguments",
        aa.seq.pattern, aa.seq.subject, sub.matrix, gap.open.pnlty,
        gap.extension.pnlty, distance.model, e, sep="\n"
      )
      write( msg, stderr() )
      NA
    }
  )
}

replaceSelenocystein <- function( aa.seq ) {
  # Function replaces all appearances of Selenocystein codes U and u with X and
  # x respectively.
  #
  # Args:
  #  aa.seq : amina acid sequence to sanitize
  #
  # Returns: Sanitized amino acid sequence where U or u is replaced with X and
  # x respectively.
  #  
  aa.seq.sanitized <- gsub( "U", "X", aa.seq, fixed=T )
  gsub( "u", "x", aa.seq.sanitized, fixed=T )
}

sequenceDistances <- function( protein.list ) {
  # Function computes pairwise sequence distances between the argument
  # proteins. It uses function 'pairwiseSequenceDistance' to achieve this. 
  #
  # Args:
  #  protein.list : A _named_ list of Proteins.
  #
  # Returns: An object of class 'dist' holding the pairwise sequence distances
  # for each pair of argument proteins.
  #   
  
  ps <- names( protein.list )
  # Iterativly compute cells of this distance matrix:
  m <- matrix( nrow=length(ps), ncol=length(ps),
    dimnames=list(ps, ps) )
  lapply( ps, function( protein.acc ) {
    # Proteins to compute distances to:
    i <- which( ps == protein.acc )
    ps2c <- ps[ i : length( ps ) ]
    lapply( ps2c, function( p2c ) {
      m[[ p2c, protein.acc ]] <<- pairwiseSequenceDistance(
        as.character( protein.list[ protein.acc ] ),
        as.character( protein.list[ p2c ] )
      )
    })
  })
  # Return an object of type 'dist':
  as.dist(m)
}

partialSequenceDistances <- function( protein.list, prot.pairs,
  lapply.funk=lapply ) {
  # Function computes pairwise sequence distances between the argument
  # proteins. It uses function 'pairwiseSequenceDistance' to achieve this.
  # Only computes distances for those protein pairs specified in argument
  # 'port.pairs'.
  #
  # Args:
  #  protein.list : A _named_ list of Proteins.
  #  prot.pairs   : Set of pairs to compute distances for.
  #  lapply.funk  : 'lapply' function to use. Set to mclapply, if parallel
  #                 computation is wanted.
  #
  # Returns: An object of class 'list' holding the pairwise sequence distances
  # for each pair of proteins speciefied in argument 'prot.pairs'. The list
  # names are goin to be the results of function 'pairwiseDistanceKey' for each
  # pair.
  #   
  setNames( 
    lapply.funk( prot.pairs, function( prot.pair ) {
      pairwiseSequenceDistance(
        as.character( protein.list[[ prot.pair[[1]] ]] ),
        as.character( protein.list[[ prot.pair[[2]] ]] )
      )
    }),
    lapply.funk( prot.pairs, function( prot.pair ) {
      pairwiseDistanceKey( prot.pair[[1]], prot.pair[[2]],
        distance.type="seq_dist" )
    })
  )
}

pairsFromBlastResult <- function( tabular.blast.out,
  query.column=1, hit.column=2 ) {
  # Constructs a list of protein pairs from the Blast results in argument
  # 'tabular.blast.out'. The first argument of each pair is the query accession
  # and the second is the hit accession. Cells are converted to type character.
  #
  # Args:
  #  tabular.blast.out : The Blast output table read using 'read.table'. 
  #  query.column      : The column index of the Blast output table which holds
  #                      the query accessions.
  #  hit.column        : The column index of the Blast output table which holds
  #                      the hit accessions.
  #
  # Returns: A list of protein pairs from the Blast results in argument
  # 'tabular.blast.out'.
  #   
  tbo.red <- tabular.blast.out[ as.character( tabular.blast.out[ , query.column ] ) != as.character( tabular.blast.out[ , hit.column ] ) , , drop=F ][ , 1:2 ] 
  unique(
    lapply( 1:nrow(tbo.red), function(i) {
      c(
        as.character( tbo.red[[ i, 1 ]] ),
        as.character( tbo.red[[ i, 2 ]] )
      )
    })
  )
}

distanceMatrixIndices <- function( accessions ) {
  # For a list of accessions - or other list names - this function
  # generates the paired indices of a distance matrix. 
  #
  # Args:
  #  accessions : to generate the paired indices of an object of type
  #               'dist' 
  #
  # Returns: The list of accession pairs to compute distances for.
  #   
  unlist(
    mclapply( accs, function( protein.acc ) {
      # Proteins to compute distances to:
      i <- which( accs == protein.acc ) + 1
      if ( i <= length( accs ) ) {
        ps2c <- accs[ i : length( accs ) ]
        lapply( ps2c, function( p2c ) {
          c( protein.acc, p2c )
        })
      }
    }),
    recursive=F
  )
}

distanceIndices <- function( batch.no, batch.size, accessions ) {
  # For argument 'batch.no' and 'batch.size' computes all protein pairs
  # of the argument batch. 
  #
  # Args:
  #  batch.no : the nth subset of 'batch.size' protein pairs
  #  batch.size : size of each subset
  #  accessions : All protein accessions to generate pairs of.
  #
  # Returns: List of indicated protein pairs.
  #   
  
  # In which row and column of the distance matrix should we start the distance
  # computation?
  row.sizes <- seq( length(accessions) - 1, 1, -1 )
  row.ind <- 0; col.ind <- 0;
  if ( batch.no == 1 ) {
    # No distance matrix cells to skip for first batch
    row.ind <- 1
    col.ind <- 1
  } else {
    # Find row and column in which to start this batch's computation:
    skip <- ( batch.no - 1 ) * batch.size
    while ( skip >= 0 ) {
      row.ind <- row.ind + 1
      col.ind <- skip + 1
      if ( row.ind <= length( row.sizes ) ) {
        skip <- skip - row.sizes[[ row.ind ]]
      } else {
        # last batch
        break
      }
    }
  }
  # Generate batch.size pairs to compute distances for
  curr.row <- accessions[ ( row.ind + col.ind ) : length(accessions) ]
  b.ind <- 1
  dist.pairs <- list()
  while ( length(dist.pairs) < batch.size ) {
    acc.a <- accessions[[ row.ind ]]
    acc.b <- curr.row[[ b.ind ]]
    dist.pairs[[ length(dist.pairs) + 1 ]] <- c( acc.a, acc.b)
    b.ind <- b.ind + 1
    if( b.ind > length( curr.row ) ) {
      row.ind <- row.ind + 1
      if ( row.ind < length(accessions) ) {
        curr.row <- accessions[ row.ind : length(accessions) ]
        b.ind <- 2
      } else {
        # Last batch
        break
      }
    }
  }
  # return
  dist.pairs
}

pairsForAccessions <- function( all.pairs, accessions, lapply.funk=lapply ) {
  # Filters the argument matrix of protein pairs 'all.pairs' for all those
  # where one pair member is contained in argument set 'accessions'. Pairs are
  # filtered for uniqueness and identity pairs like ( a, a ) are excluded.
  #
  # Args:
  #  all.pairs : The 2 column matrix or dataframe in which the set of protein
  #              pairs is stored. Partner one is supposed to be in column one
  #              and partner two in column two.
  #  accessions : The set of unique protein accessions to extract all existing
  #               pairs for.
  #  lapply.funk : Set to 'mclapply' if parallel execution is wanted.
  #
  # Returns: The subset of argument 'all.pairs' where each pair contains at
  # least a single accession of argument set 'accessions'. Returned pairs are
  # unique and do not contain identity pairs.
  #   
  uniquePairs(
    all.pairs[ all.pairs[,1] %in% accessions | all.pairs[,2] %in% accessions, , drop=F ],
    lapply.funk=lapply.funk
  )
}

pairsForAccessionsAssumingSymmetry <- function( all.pairs, accessions,
  lapply.funk=lapply ) {
  # Filters the argument matrix of protein pairs 'all.pairs' for all those where
  # pair member in column one is contained in argument set 'accessions'. This
  # method assumes that the table 'all.pairs' is i.e. the result of a
  # symmetrical 'all vs all' sequence similarity search, whose result would
  # have each pair twice once with member A in the first column and the second
  # time in the second column. Pairs are filtered for uniqueness and identity
  # pairs like ( a, a ) are excluded.
  #
  # Args:
  #  all.pairs : The 2+ column matrix or dataframe in which the set of protein
  #              pairs is stored. Partner one is supposed to be in column one
  #              and partner two in column two.
  #  accessions : The set of unique protein accessions to extract all existing
  #               pairs for.
  #  lapply.funk : Set to 'mclapply' if parallel execution is wanted.
  #
  # Returns: The subset of argument 'all.pairs' where each pair contains at
  # least a single accession of argument set 'accessions'. Returned pairs are
  # unqique and do not contain identity pairs.
  #   
  uniquePairs(
    all.pairs[ all.pairs[,1] %in% accessions, , drop=F ],
    lapply.funk=lapply.funk
  )
}

uniquePairs <- function( pairs.tbl, lapply.funk=lapply ) {
  # Each row - column 1 and 2 - of argument matrix 'pairs.tbl' is sorted
  # alphabetically. The set of so sorted pairs is then filtered for _unique_
  # pairs. So a symmetrical pair ( a, b ) and ( b, a ) will only appear once in
  # the result as pair ( a, b ). Also identity pairs like ( a, a ) will be
  # discarded.
  #
  # Args:
  #  pairs.tbl : The table of symmetrical pairs with member one in column 1 and
  #              member two in column 2.
  #  lapply.funk : Set to mclapply, if parallel execution is wanted.
  #
  # Returns: A two column matrix of pairs, whose members are alphabetically
  # sorted.
  #   
  if( ! is.null( pairs.tbl ) &&
    nrow( pairs.tbl ) > 0 && ncol( pairs.tbl ) >= 2 ) {
    do.call( 'rbind', 
      unique(
        lapply.funk( 1:nrow(pairs.tbl), function(i) {
          acc.a <- as.character( pairs.tbl[[ i, 1 ]] )
          acc.b <- as.character( pairs.tbl[[ i, 2 ]] )
          if( ! identical( acc.a, acc.b ) )
            sort( c( acc.a, acc.b ) )
          else
            NULL
        })
      )
    )
  }
}

