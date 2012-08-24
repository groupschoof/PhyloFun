library(RUnit)
library(tools)
library(ape)
# In R sourcing other files is not trivial, unfortunately.
# WARNING:
# This method ONLY works for project files in depth one sub dirs!
project.file.path <- function(...) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dir <- dirname(file_path_as_absolute(script.name))
  project.dir <- sub(basename(script.dir),'',script.dir)
  normalizePath(file.path(project.dir,...))
}
src.project.file <- function(...) {
  source(project.file.path(...))
}
# We set-up required libraries in the test case, not in the R file, as path
# problems will be resolved, as soon as this R package is loaded as such.
src.project.file('src','loadUniprotKBEntries.R')
src.project.file('src','phyloFun.R')

# Test tree is midpoint rooted!
phylo.tree <- read.tree(project.file.path('test', 'test_tree.newick'))
fl <- file(project.file.path('test','test_annotations.tbl'),"r")
annotation.matrix <- unserialize(fl)
close(fl)

# Test conditional.probs.tbl
print("Testing conditional.probs.tbl(...)")
err <- try(conditional.probs.tbl(0.45,NA,NA),silent=T)
checkTrue(identical(class(err),"try-error"))
# err <- try(conditional.probs.tbl(0.45,NA),silent=T)
# checkTrue(identical(class(err),"try-error"))
ua <- uniq.annotations(annotation.matrix,'GO')
res <- conditional.probs.tbl(0.45,
  uniq.annotations(annotation.matrix,'GO'))
checkTrue(identical(dim(res),as.integer(c(length(ua),length(ua)))))
checkEquals(rownames(res),ua)
checkEquals(colnames(res),ua)

# Test get.node.label
print("Testing get.node.label(...)")
checkEquals(get.node.label(phylo.tree,21),21L)
checkTrue(identical(class(get.node.label(phylo.tree,11)),"character"))

# Test edge.to.formula
print("Testing edge.to.formula(...)")
# Test a tip's formula
indx <- which(phylo.tree$edge[,2] == 1)
frml <- edge.to.formula(phylo.tree,indx)
checkTrue(identical(class(frml),'formula'))
checkTrue(length(as.character(frml)) == 2)
checkTrue(grepl('[a-zA-Z]+',as.character(frml)[[2]],perl=T))
indx <- which(phylo.tree$edge[,1] == (length(phylo.tree$tip.label)+1))[1]
frml <- edge.to.formula(phylo.tree,indx)
checkTrue(! grepl('[a-zA-Z]+',as.character(frml)[[2]],perl=T))

# Test bayes.nodes
print("Testing bayes.nodes(...)")
bys.nds <- bayes.nodes(phylo.tree,annotation.matrix)
checkTrue(length(bys.nds) == nrow(phylo.tree$edge)+1)
root.bys.nd <- bys.nds[[1]]
checkTrue(length(root.bys.nd$values) == length(ua))
checkEquals(root.bys.nd$levels,ua)
# print(bys.nds[[1]])

# Test create bayesian Network
print("Testing create bayesian Network(...)")
grain(compileCPT(bys.nds))


