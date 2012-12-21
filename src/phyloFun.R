library(gRain)
library(Matrix)

# Returns unique node, that has no edge leading to it. Thus, it must be a node
# in the first column whereas mustn't appear in the second column of the
# phylTree's edge matrix.
get.root.node <- function(phylTree) {
  as.character(
    unique(phylTree$edge[!(phylTree$edge[,1] %in% phylTree$edge[,2]),][,1])[[1]]
    )
}

get.node.label <- function(phylTree,node.index) {
  res <- try(phylTree$tip.label[[node.index]],silent=T)
  if(identical("try-error",class(res))) as.character(node.index) else res
}

edge.to.formula <- function(phyloTree,edge.index) {
  prnt <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index,1]])
  child <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index,2]])
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

conditional.probs.tbl <- function( edge.length, uniq.annotations,
  annots.mut.prob.table.list, mut.tbl.length.col.indx=5, p.mut.col.indx=1 ) {
  # Looks up the mutation probability tables for each annotation and constructs
  # the transition probabilities based on the current branch's length.
  #
  # Args:
  #  edge.length : sequence distance as the current branch's length.
  #  uniq.annotations : The unique anntotations as they appear in the current
  #                     phylogenetic tree.
  #  annots.mut.prob.table.list : The list of empirical mutation probabilities
  #                               for each annotation.
  #  mut.tbl.length.col.indx : The index of the mutation probability table's
  #                            column in which to find the corresponding
  #                            sequence distance.
  #  p.mut.col.indx : The mutation probability column's index.
  #
  # Returns: 
  #   A length(uniq.annotations)*length(uniq.annotations) probability matrix with
  #   the m[i,j] := P(j | i), i = parent function, j = child function
  #   
  p.retain.lst <- as.numeric( 
    lapply( uniq.annotations, function( a ) {
      p.mut.tbl <- annots.mut.prob.table.list[[ a ]]
      if ( ! is.null( p.mut.tbl ) ) {
        1 - findMatchingColumn(
          p.mut.tbl, edge.length, mut.tbl.length.col.indx
        )[[ 1, p.mut.col.indx ]]
      } else {
        0.5
      }
    })
  )
  names( p.retain.lst ) <- uniq.annotations
  print( p.retain.lst )
  # ToDo: Do we really need a double iteration?
  do.call( 'rbind', 
    lapply( uniq.annotations, function( a ) {
      p.ret <- p.retain.lst[[ a ]]
      p.mut <- 1 - p.ret
      norm.trans.probs <- p.mut / sum( as.numeric (
        p.retain.lst[ which( names(p.retain.lst[]) != a ) ]
      ) )
      trans.probs <- p.retain.lst * norm.trans.probs
      trans.probs[[ a ]] <- p.ret
      matrix( trans.probs, nrow=1, dimnames=list( a, uniq.annotations ) )
    })
  )
}

bayes.nodes <- function(phylo.tree,annotation.matrix,annotation.type='GO') {
  annotations <- uniq.annotations(annotation.matrix,annotation.type)
  no.evidence <- as.numeric(matrix(1.0,nrow=1,ncol=length(annotations)))
  nds <- list(
    cptable(
      eval(bquote(~ .(get.root.node(phylo.tree)))),
      values=no.evidence,
      levels=annotations))
  for(i in 1:nrow(phylo.tree$edge)) {
    nds <- append(nds,
      list(cptable(edge.to.formula(phylo.tree,i),
        values=conditional.probs.tbl(
          phylo.tree$edge.length[[i]],
          annotations),
        levels=annotations))
      )
  }
  nds
}

# Returns the labels of phylo.tree's tips, that have no annotation of specified
# annotation.type, or those that have non NA annotations, if negate=TRUE.
getTipsWithNaAnnotation <- function(phylo.tree,annotation.matrix,annotation.type='GO',negate=F) {
  bool.indxs <- is.na(annotation.matrix[annotation.type,])
  surroundEachWithQuotes(
    colnames(annotation.matrix)[which( if(negate) ! bool.indxs else bool.indxs)]
  )
}

# Each entry of char.vector is surrounded by "\"".
surroundEachWithQuotes <- function(char.vector) {
  sapply(char.vector,
    function(x){
      paste("\"",x,"\"",sep='')
    },
    USE.NAMES=F)
}

queryPhylBayesNetwork <- function(
  phylo.tree,
  annotation.matrix,
  annotation.type='GO',
  bayes.netw=grain(
    compileCPT(bayes.nodes(phylo.tree,
        annotation.matrix,
        annotation.type))),
  type='distribution') {
  # Provide diagnostic evidence matrix (true 'function annotation'):
  evidence.matrix <- annotation.matrix[!is.na(annotation.matrix[,annotation.type]),annotation.type]
  colnames(evidence.matrix) <- surroundEachWithQuotes(colnames(evidence.matrix))
  list(prediction=predict(bayes.netw,
    getTipsWithNaAnnotation(phylo.tree,annotation.matrix,annotation.type),
    getTipsWithNaAnnotation(phylo.tree,annotation.matrix,annotation.type,negate=T),
    evidence.matrix,
    type),
    bayesian.network=bayes.netw)
}
