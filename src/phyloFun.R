phylo.fun <- function(phylo.tree,annotation.matrix) {

}

# Returns unique node, that has no edge leading to it. Thus, it must be a node
# in the first column whereas mustn't appear in the second column of the
# phylTree's edge matrix.
getRootNode <- function(phylTree) {
  unique(phylTree$edge[!(phylTree$edge[,1] %in% phylTree$edge[,2]),][,1])
}

