# Test isConnectionAlive
library(RUnit)
library(tools)
library(RMySQL)

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
src.project.file( 'src', 'geneOntologySQL.R' )

print("Testing isConnectionAlive(...)")
go.con <- connectToGeneOntology() 
res.isConnectionAlive <- isConnectionAlive( go.con )
checkTrue( res.isConnectionAlive )
dbDisconnect( go.con )
checkTrue( ! isConnectionAlive( go.con ) )
