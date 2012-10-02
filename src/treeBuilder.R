buildVectorSpaceModel <- function(uniq.annotations) {
  # Just removes the NA entries from the set of unique annotations.
  #
  # Args:
  #  uniq.annotations : Result from calling the function uniq.annotations.
  #
  # Returns: Argument set without NA entries.
  uniq.annotations[ ! is.na(uniq.annotations[]) ]
}

domainSpaceVectors <- function(vector.space.model, annotation.matrix,
  annotation.type="InterPro") {
  # Constructs a list of vectors in the vector space of conserved protein
  # domains, where each annotated protein is assigned a vector with the
  # appropriate domain weight at position i, if and only if the protein is
  # annotated with domain i.
  #
  # Args:
  #  vector.space.model : A character vector of alphabetically sorted unique
  #                       conserved protein domains (i.e. InterPro or Pfam )
  #  annotation.matrix  : The output of function retrieve.annotations.parallel
  #                       which is a matrix with the Protein Accessions as 
  #                       columns and the annotated conserved domains as rows.
  #  annotation.type    : Flag indicating what type of conserved protein
  #                       domains to use 'InterPro' or 'Pfam' ..
  #
  # Returns: A list with the protein accessions as names and their domain
  # vectors as values.
  

}
