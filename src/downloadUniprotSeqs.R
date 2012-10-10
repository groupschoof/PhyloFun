library(tools)
library(RCurl)
library(XML)

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
src.project.file( 'src','loadUniprotKBEntries.R' )

input.args <- commandArgs( trailingOnly = TRUE )

print("Usage: Rscript downloadUniprotSeqs.R path/2/accession_per_line.txt path/2/output.fasta")

# Read input
accs <- scan( input.args[[1]], what=character(), sep="\n" )

# Generate Uniprot-URLs
uris <- lapply( accs, uniprotkb.url )

# Download docs
docs <- getURL( uris )

# Attempt to download 'Server Too Busy' docs again:
busy.uris <- names( findServerBusyResults( docs ) )
if ( length(busy.uris) > 0 ) {
  # Wait a couple of seconds:
  Sys.sleep( sample(1:90, 1) )
  # Now try again:
  docs <<- c(
    docs[ names(docs) != busy.uris ],
    getURL( busy.uris )
  )
}

# Extract sequences
seqs <- sapply( 
  docs, function(d) retrieveSequence( xmlInternalTreeParse(d) ),
  USE.NAMES=F
)

# Generate output
o <- sapply( 1:length(accs), function(i) {
    paste(
      paste( ">", accs[[i]], sep="" ),
      seqs[[i]], sep="\n"
    )
  },
  USE.NAMES=F 
)

# Write output
fasta <- file( input.args[[2]], "w" )
writeLines( o, fasta )
