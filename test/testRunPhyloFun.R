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

system( 
  paste( 'Rscript', project.file.path( 'src', 'runPhyloFun.R' ),
    '-q', project.file.path( 'test', 'protein_1.fasta' ),
    '-j', project.file.path( 'test', 'jackhmmer_out_10_prots.tbl' )
  )
)
