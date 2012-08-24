library(gRain)
library(Matrix)

bayes.nodes <- function(phylo.tree,annotation.matrix,type='GO') {
  annotations <- uniq.annotations(annotation.matrix,type)
  no.evidence <- as.numeric(matrix(1.0,nrow=1,ncol=length(annotations)))
  nds <- list(
    cptable(
      eval(bquote(~ .(as.integer(get.root.node(phylo.tree))))),
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


# Returns unique node, that has no edge leading to it. Thus, it must be a node
# in the first column whereas mustn't appear in the second column of the
# phylTree's edge matrix.
get.root.node <- function(phylTree) {
  as.integer(
    unique(phylTree$edge[!(phylTree$edge[,1] %in% phylTree$edge[,2]),][,1])[[1]]
    )
}

get.node.label <- function(phylTree,node.index) {
  res <- try(phylTree$tip.label[[node.index]],silent=T)
  if(identical("try-error",class(res))) node.index else res
}

edge.to.formula <- function(phyloTree,edge.index) {
  prnt <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index,1]])
  child <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index,2]])
  eval(bquote(~ .(child) | .(prnt) ))
}

conditional.probs.tbl <- function(
  edge.length,
  uniq.annotations,
  base.transition.probs=matrix(0.03,
    nrow=length(uniq.annotations),
    ncol=length(uniq.annotations),
    dimnames=list(uniq.annotations,uniq.annotations)
    )) {
  # Validate arguments:
  if(identical(uniq.annotations,NA) &&
    identical(base.transition.probs,NA))
    stop( "Illegal arguments: Either argument 'uniq.annotations' or 'base.transition.probs' must be set.")
  # Arguments are fine:
  expm(edge.length*base.transition.probs)
}
