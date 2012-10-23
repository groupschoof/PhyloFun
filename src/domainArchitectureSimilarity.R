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
  annotation.matrix, domain.weights.table, annotation.type='InterPro' ) {
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
  #
  # Returns: List of domain architecture space vectors one for each annotated
  # protein.
  #   
  amr <- if ( is.null(annotation.type) ) 1 else annotation.type

  do.call('cbind',
    setNames(
      mclapply( colnames(annotation.matrix), function( protein.id ) {
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
        },
      mc.cores=detectCores(), mc.preschedule=T),
    colnames( annotation.matrix ))
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

partialDomainArchitectureDistances <- function( domain.architecture.space.vectors,
  accessions ) {
  # Method computes pairwise distances of the domain architecture space vectors
  # given in argument matrix 'domain.architecture.space.vectors'. Column names
  # are expected to be the Proteins and row names the Domains, hence the ith
  # column is the vector of protein i. 
  #
  # Args:
  #  domain.architecture.space.vectors : The matrix of vectors in domain
  #                                      architecture space.
  #
  # Returns: Returns an object of class 'matrix', holding the pairwise distances of
  # the argument vectors.
  #   
  
  # All proteins to compute pairwise distances for
  ps <- colnames( domain.architecture.space.vectors )
  # Iterativly compute cells of this distance matrix:
  m <- matrix( nrow=length(ps), ncol=length(ps),
    dimnames=list(ps, ps) )
  lapply( accessions, function( protein.acc ) {
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
  # Return an object of type 'matrix':
  m
}

pairwiseSequenceDistance <- function( aa.seq.pattern, aa.seq.subject, sub.matrix="PAM250",
  gap.open.pnlty=-10, gap.extension.pnlty=-0.1, distance.model="Dayhoff") {
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
  pa <- pairwiseAlignment(
    AAString( replaceSelenocystein( aa.seq.pattern ) ),
    AAString( replaceSelenocystein( aa.seq.subject ) ),
    substitutionMatrix=sub.matrix, gapOpening=gap.open.pnlty, gapExtension=gap.extension.pnlty )
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

partialSequenceDistances <- function( protein.list, accessions ) {
  # Function computes pairwise sequence distances between the argument
  # proteins. It uses function 'pairwiseSequenceDistance' to achieve this. 
  # Only computes distances between proteins specified by accessions and all
  # proteins in list.
  #
  # Args:
  #  protein.list : A _named_ list of Proteins.
  #
  # Returns: An object of class 'matrix' holding the pairwise sequence distances
  # for each pair of proteins where one member is in the argument 'accessions'.
  #   
  
  ps <- names( protein.list )
  # Iterativly compute cells of this distance matrix:
  m <- matrix( 0, nrow=length(accessions), ncol=length(ps),
    dimnames=list(accessions, ps) )
  lapply( accessions, function( protein.acc ) {
    # Proteins to compute distances to:
    i <- which( ps == protein.acc )
    ps2c <- ps[ ps != protein.acc ]
    lapply( ps2c, function( p2c ) {
      m[[ protein.acc, p2c ]] <<- pairwiseSequenceDistance(
        as.character( protein.list[ protein.acc ] ),
        as.character( protein.list[ p2c ] )
      )
    })
  })
  # Return an object of type 'matrix':
  m
}
