library(gRain)
library(Matrix)

# Returns unique node, that has no edge leading to it. Thus, it must be a node
# in the first column whereas mustn't appear in the second column of the
# phylTree's edge matrix.
get.root.node <- function(phylTree) {
  as.character(
    unique(phylTree$edge[!(phylTree$edge[, 1] %in% phylTree$edge[, 2]), ][, 1])[[1]]
    )
}

get.node.label <- function(phylTree, node.index) {
  res <- try(phylTree$tip.label[[node.index]], silent=T)
  if(identical("try-error", class(res))) as.character(node.index) else res
}

edge.to.formula <- function(phyloTree, edge.index) {
  prnt <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index, 1]])
  child <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index, 2]])
  eval(bquote(~ .(child) | .(prnt) ))
}

findMatchingColumn <- function( mutation.probability.table, value, column.index ) {
  # Looks up the _last_ row of argument 'mutation.probability.table' who's
  # entry in column 'column.index' is smaller or equal to argument 'value'. If
  # none can be found returns the last row of mutation.probability.table. 
  #
  # Args:
  #  mutation.probability.table : The matrix holding mutation probabilities as
  #                               returned i.e. by function
  #                               'mutationProbabilityDistribution'.
  #  value                      : The current distance value, i.e. the branch
  #                               length to find the correct mutation
  #                               probability for.
  #  column.index               : The index of the column holding the distance
  #                               measure to compare argument value with.
  #
  # Returns: The argument mutation.probability.table's matching row as matrix.
  #   
  hits <- mutation.probability.table[
    which( mutation.probability.table[ , column.index ] >= value ), , drop=F
  ]
  if ( nrow( hits ) > 0 ) {
    hits[ 1, , drop=F ]
  } else {
    mutation.probability.table[ nrow(mutation.probability.table), , drop=F ]
  }
}

validAnnotations <- function( annotation.matrix, unknown.annot='unknown',
  annotation.type='GO', annotations.with.mutation.probability.tables=names(
  GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE ) ) {
  # Filters all annotations of argument type 'annotation.type' to be pairwise
  # unique, to exclude NAs, and to have mutation probability tables set in
  # argument 'mutation.probability.tables.list'. To this list the UNKNOWN
  # annotation is appended.
  #
  # Args:
  #  annotation.matrix : The matrix of protein function annotations as returned
  #                      i.e. by 'retrieveAnnotationsBiomart'.
  #  unknown.annot : The UNKNOWN annotation, can be set to any convenient name.
  #  annotation.type : The type of annotation to extract, must be one of
  #                    argument 'annotation.matrix' row names.
  #  annotations.with.mutation.probability.tables : The annotations for which
  #                    mutation probability tables have been generated. Default
  #                    is the names of constant list
  #                    GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE.
  #
  # Returns: A character vector of legal and valid annotations including the
  # UNKNOWN.
  #   
  c(
    intersect(
      uniq.annotations( annotation.matrix, annotation.type, exclude.NAs=T ),
      annotations.with.mutation.probability.tables
    ),
    unknown.annot
  )
}

conditional.probs.tbl <- function( edge.length, annos,
  annots.mut.prob.table.list, mut.tbl.length.col.indx=5, p.mut.col.indx=1,
  unknown.annot='unknown' ) {
  # Looks up the mutation probability tables for each annotation and constructs
  # the transition probabilities based on the current branch's length.
  # Introduced a new annotation UNKNOWN which can mutate to every other
  # annotation in argument 'uniq.annotations' and does have a zero probability
  # of self retainment.
  #
  # Args:
  #  edge.length : sequence distance as the current branch's length.
  #  annos : The unique anntotations as they appear in the current phylogenetic
  #          tree. Should NOT contain NULL or NA values. Use function
  #          'validAnnotations' to generate this argument.
  #  annots.mut.prob.table.list : The list of empirical mutation probabilities
  #                               for each annotation.
  #  mut.tbl.length.col.indx : The index of the mutation probability table's
  #                            column in which to find the corresponding
  #                            sequence distance.
  #  p.mut.col.indx : The mutation probability column's index.
  #  unknown.annot : The UNKNOWN annotation, set to any convenient name.
  #
  # Returns: 
  #   A length(annos)*length(annos) probability matrix with
  #   the m[i, j] := P(j | i), i = parent function, j = child function
  #   
  do.call( 'rbind', 
    lapply( annos, function( a ) {
      p.mut.tbl <- annots.mut.prob.table.list[[ a ]]
      if ( identical( a, unknown.annot ) ) {
        trans.probs <- matrix( 1 / ( length( annos ) - 1 ),
          nrow=1, ncol=length( annos ),
          dimnames=list( a, annos )
        )
        # unknown retainment probability:
        trans.probs[[ a, a ]] <- 0
        trans.probs
      } else if ( ! is.null( p.mut.tbl ) ) {
        p.mut <- findMatchingColumn(
          p.mut.tbl, edge.length, mut.tbl.length.col.indx
          )[[ 1, p.mut.col.indx ]]
        trans.probs <- matrix( ( p.mut / ( length( annos ) - 1 ) ),
          nrow=1, ncol=length( annos ),
          dimnames=list( a, annos )
        )
        # annotation retainment probability:
        trans.probs[[ a, a ]] <- 1 - p.mut
        trans.probs
      } else {
        write( paste( "WARNING! Could not find mutation probability table for", a ),
          stderr() )
        NULL
      }
    })
  )
}

bayesNodes <- function( phylo.tree, annotation.matrix, annotation.type='GO',
  mutation.probability.tables.list=GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE,
  unknown.annot='unknown'
  ) {
  # Constructs a bayesian network from the argument phylogenetic tree assigns
  # each node the annotation mutation probability distribution. 
  #
  # Args:
  #  phylo.tree                       : Argument phylogenetic tree as returned
  #                                     by funtion 'read.tree'.
  #  annotation.matrix                : The function annotations for the leaves
  #                                     in the phylogenetic tree.
  #  mutation.probability.tables.list : The list of conditional annotation
  #                                     mutation probabilities as generated by
  #                                     function 'mutationProbabilityDistribution'.
  #  unknown.annot                    : UNKNOWN annotation, can be set to any
  #                                     convenient name.
  #
  # Returns: a list of bayesian network nodes as generated by function
  # 'cptable' from library 'gRain'.
  #   
  annotations <- validAnnotations( annotation.matrix, unknown.annot,
    annotation.type, names( mutation.probability.tables.list )
  )
  no.evidence <- as.numeric( matrix( 1.0, nrow=1, ncol=length( annotations ) ) )
  nds <- list(
    cptable( eval( bquote( ~ .( get.root.node( phylo.tree ) ) ) ),
      values=no.evidence, levels=annotations )
  )
  for( i in 1:nrow( phylo.tree$edge ) ) {
    nds <- append( nds,
       list( cptable( edge.to.formula( phylo.tree, i ),
           values=conditional.probs.tbl( phylo.tree$edge.length[[i]], annotations,
             mutation.probability.tables.list, unknown.annot=unknown.annot ),
           levels=annotations ) )
    )
  }
  nds
}

# Returns the labels of phylo.tree's tips, that have no annotation of specified
# annotation.type, or those that have non NA annotations, if negate=TRUE.
getTipsWithNaAnnotation <- function(phylo.tree, annotation.matrix, annotation.type='GO', negate=F) {
  bool.indxs <- is.na(annotation.matrix[annotation.type, ])
  surroundEachWithQuotes(
    colnames(annotation.matrix)[which( if(negate) ! bool.indxs else bool.indxs)]
  )
}

# Each entry of char.vector is surrounded by "\"".
surroundEachWithQuotes <- function(char.vector) {
  sapply(char.vector,
    function(x){
      paste("\"", x, "\"", sep='')
    },
    USE.NAMES=F)
}

queryPhylBayesNetwork <- function(
  phylo.tree,
  annotation.matrix,
  annotation.type='GO',
  bayes.netw=grain(
    compileCPT(bayesNodes(phylo.tree,
        annotation.matrix,
        annotation.type))),
  type='distribution') {
  # Provide diagnostic evidence matrix (true 'function annotation'):
  evidence.matrix <- annotation.matrix[!is.na(annotation.matrix[, annotation.type]), annotation.type]
  colnames(evidence.matrix) <- surroundEachWithQuotes(colnames(evidence.matrix))
  list(prediction=predict(bayes.netw,
    getTipsWithNaAnnotation(phylo.tree, annotation.matrix, annotation.type),
    getTipsWithNaAnnotation(phylo.tree, annotation.matrix, annotation.type, negate=T),
    evidence.matrix,
    type),
    bayesian.network=bayes.netw)
}
