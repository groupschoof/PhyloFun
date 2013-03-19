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

fastTreeCall <- if( try( system( 'FastTreeMP' ), silent=T ) == 0 ) 'FastTreeMP' else 'FastTree'

system( 
  paste( 'Rscript', project.file.path( 'src', 'runPhyloFun.R' ),
    '-q', project.file.path( 'test', 'protein_1.fasta' ),
    '-j', project.file.path( 'test', 'protein_1_jackhmmer_out.tbl' ),
    '-f', fastTreeCall
  )
)

checkTrue( file.exists( project.file.path( 'Protein_1', 'go_term_predictions.tbl' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'ml_tree.newick' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta-gb' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta-gb-phylo_filtered' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'phyloFun_R_serialized.txt' ) ) )

# clean up:
unlink( "Protein_1", recursive=T )
