library(RUnit)
library(tools)
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
src.project.file('src','loadUniprotKBEntries.R')
src.project.file('src','treeBuilder.R')

# Test buildVectorSpaceModel
print("Testing buildVectorSpaceModel(...)")
f <- file( project.file.path("test", "test_annotations.tbl"), "r" )
annos <- unserialize(f)
close(f)
ua <- uniq.annotations(annos, type='InterPro')
vsm <- buildVectorSpaceModel(ua)
checkEquals( vsm,
  c( "IPR001957", "IPR003593", "IPR010921", "IPR013159", "IPR013317",
    "IPR018312", "IPR020591", "IPR024633" ))

# Read test domain weights 
domain.weights <- read.table(project.file.path("test", "domainWeights.tbl"))
