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

phylo.tree <- read.tree(project.file.path('test', 'test_tree.newick'))
print(phylo.tree$edge)
annotation.matrix <- read.table(project.file.path('test', 'test_annotations.tbl'))
print(annotation.matrix)

