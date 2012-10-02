constructVectorSpaceModel <- function( annotation.matrix, type='InterPro' ) {
  # Constructs a character vector of alphabetically sorted annotations IDs,
  # excluding NA annotations. Uses function 'uniq.annotations' to process the
  # 'annotation.matrix'.
  #
  # Args:
  #  annotation.matrix : Rows are annotation types and columns the annotated
  #                      proteins, as returned by function
  #                      'retrieve.annotations.parallel.'
  #  type : Row name of 'annotation.matrix' to be processed.
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

  sapply( colnames(annotation.matrix), function( protein.id ) {
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
  })
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
