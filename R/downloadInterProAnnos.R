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

# Usage:
print( "Usage: Rscript downloadInterProAnnos.R path/2/accession_per_line.txt path/2/output.tbl")

# Input
trailing.args <- commandArgs(trailingOnly = TRUE)

# Read Accessions:
accs <- as.character( read.table( trailing.args[[1]] )$V1 )

# Iteratively download InterPro Annotations for 100 Accessions:
u <- useDataset( "uniprot", mart=useMart("unimart") )
no.batches <- ceiling( (length(accs) / 100) )

all.annos <- do.call( 'cbind', 
  lapply(1:no.batches, function(i) { 
    star.ind <- if ( i == 1 ) i else (i - 1) * 100
    stop.ind <- if ( (i * 100 - 1) < length(accs) ) (i * 100 - 1) else length(accs)
    accs.batch <- accs[star.ind:stop.ind]
    retrieveAnnotationsBiomart( accs.batch, uni.mart=u ) 
  })
)

# Save result:
f <- file( trailing.args[[2]], "w" )
serialize( all.annos, f )
close( f )

print( "DONE" )
