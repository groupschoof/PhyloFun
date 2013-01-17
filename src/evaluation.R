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
