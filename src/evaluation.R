precision <- function( predicted.gos, true.gos ) {
  # 'precision' is a statistical quality measure of predicted Gene Ontology
  # (GO) terms:
  # precision( predicted.gos ) = | true_positives | / | predicted.gos |
  #
  # Args:
  #  predicted.gos : The set of predicted GO terms
  #  true.gos      : The set of experimentally verified GO terms.
  #
  # Returns: The precision of the prediction as a numeric between Zero and One.
  # Precision for empty true.gos is always One.
  #   
  pgs <- unique( predicted.gos )
  tgs <- unique( true.gos )
  tp <- intersect( tgs, pgs )
  if ( length( pgs ) == 0 ) 1 else length( tp ) / length( pgs )
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

fScore <- function( predicted.gos, true.gos, beta.param=1 ) {
  # Computes the weighted harmonic mean of precision and recall as a quality
  # measure of predicted.gos. The higher beta.param the more emphasis is put on
  # recall than on precision.
  #
  # Args:
  #  predicted.gos : The set of predicted GO terms
  #  true.gos      : The set of experimentally verified GO terms.
  #  beta.param    : The factor of how much more emphasis should be put on
  #                  recall than on precision.
  #
  # Returns: The weighted harmonic mean of precision and recall as a numeric
  # between Zero and One.
  #   
  prcsn <- precision( predicted.gos, true.gos )
  rcll <- recall( predicted.gos, true.gos )
  bp <- beta.param^2
  if ( 0 == (prcsn + rcll) )
    0
  else
    ( 1 + bp ) * ( prcsn * rcll ) / ( bp * prcsn + rcll )
}
