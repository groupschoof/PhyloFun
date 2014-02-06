require( PhyloFun )

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
    '-q', project.file.path(  'protein_1.fasta' ),
    '-p', project.file.path(  'protein_1_jackhmmer_out.tbl' ),
    '-f', fastTreeCall, '-h true', '-m true'
  )
)

checkTrue( file.exists( project.file.path( 'Protein_1', 'go_term_predictions.tbl' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'ml_tree.newick' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta-gb' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs_stats.txt' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa_stats.txt' ) ) )

# clean up:
unlink( "Protein_1", recursive=T )

# Test PhyloFun with Blast results and not filtering for EVIDENCE CODES:
system( 
  paste( 'Rscript', project.file.path( 'src', 'runPhyloFun.R' ),
    '-q', project.file.path(  'protein_1.fasta' ),
    '-b', project.file.path(  'protein_1_blastout.tbl' ),
    '-f', fastTreeCall, '-h true', '-m true', '-e ALL'
  )
)

checkTrue( file.exists( project.file.path( 'Protein_1', 'go_term_predictions.tbl' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'ml_tree.newick' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta-gb' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs_stats.txt' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa_stats.txt' ) ) )

# clean up:
unlink( "Protein_1", recursive=T )

# Test PhyloFun with HTML reporting:
system( 
  paste( 'Rscript', project.file.path( 'src', 'runPhyloFun.R' ),
    '-q', project.file.path(  'protein_1.fasta' ),
    '-b', project.file.path(  'protein_1_blastout.tbl' ),
    '-f', fastTreeCall, '-h true', '-m true', '-e ALL', '-r true'
  )
)

checkTrue( file.exists( project.file.path( 'Protein_1', 'go_term_predictions.tbl' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'ml_tree.newick' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa.fasta-gb' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'homologs_stats.txt' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'msa_stats.txt' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'report', 'Protein_1_report.html' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'report', 'Protein_1_phylo_fun_tree.newick' ) ) )
checkTrue( file.exists( project.file.path( 'Protein_1', 'report', 'Protein_1_phylo_fun_tree.png' ) ) )

# clean up:
unlink( "Protein_1", recursive=T )