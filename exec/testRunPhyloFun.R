require( PhyloFun )

fastTreeCall <- if( try( system( 'FastTreeMP' ), silent=T ) == 0 ) 'FastTreeMP' else 'FastTree'

system( 
  paste( 'Rscript', project.file.path( 'exec', 'runPhyloFun.R' ),
    '-q', project.file.path(  'protein_1.fasta' ),
    '-p', project.file.path(  'protein_1_jackhmmer_out.tbl' ),
    '-f', fastTreeCall, '-h true', '-m true'
  )
)

checkTrue( file.exists( 'Protein_1/go_term_predictions.tbl' ) ) 
checkTrue( file.exists( 'Protein_1/homologs.fasta' ) ) 
checkTrue( file.exists( 'Protein_1/ml_tree.newick' ) ) 
checkTrue( file.exists( 'Protein_1/msa.fasta' ) ) 
checkTrue( file.exists( 'Protein_1/msa.fasta-gb' ) ) 
checkTrue( file.exists( 'Protein_1/homologs_stats.txt' ) ) 
checkTrue( file.exists( 'Protein_1/msa_stats.txt' ) ) 

# clean up:
unlink( "Protein_1", recursive=T )

# Test PhyloFun with Blast results and not filtering for EVIDENCE CODES:
system( 
  paste( 'Rscript', project.file.path( 'exec', 'runPhyloFun.R' ),
    '-q', project.file.path(  'protein_1.fasta' ),
    '-b', project.file.path(  'protein_1_blastout.tbl' ),
    '-f', fastTreeCall, '-h true', '-m true', '-e ALL'
  )
)

checkTrue( file.exists( 'Protein_1/go_term_predictions.tbl' ) )
checkTrue( file.exists( 'Protein_1/homologs.fasta' ) ) 
checkTrue( file.exists( 'Protein_1/ml_tree.newick' ) ) 
checkTrue( file.exists( 'Protein_1/msa.fasta' ) ) 
checkTrue( file.exists( 'Protein_1/msa.fasta-gb' ) ) 
checkTrue( file.exists( 'Protein_1/homologs_stats.txt' ) ) 
checkTrue( file.exists( 'Protein_1/msa_stats.txt' ) ) 

# clean up:
unlink( "Protein_1", recursive=T )

# Test PhyloFun with HTML reporting:
system( 
  paste( 'Rscript', project.file.path( 'exec', 'runPhyloFun.R' ),
    '-q', project.file.path(  'protein_1.fasta' ),
    '-b', project.file.path(  'protein_1_blastout.tbl' ),
    '-f', fastTreeCall, '-h true', '-m true', '-e ALL', '-r true'
  )
)

checkTrue( file.exists( 'Protein_1/go_term_predictions.tbl' ) ) 
checkTrue( file.exists( 'Protein_1/homologs.fasta' ) ) 
checkTrue( file.exists( 'Protein_1/ml_tree.newick' ) ) 
checkTrue( file.exists( 'Protein_1/msa.fasta' ) ) 
checkTrue( file.exists( 'Protein_1/msa.fasta-gb' ) ) 
checkTrue( file.exists( 'Protein_1/homologs_stats.txt' ) ) 
checkTrue( file.exists( 'Protein_1/msa_stats.txt' ) ) 
checkTrue( file.exists( 'Protein_1/report/Protein_1_report.html' ) ) 
checkTrue( file.exists( 'Protein_1/report/Protein_1_phylo_fun_tree.newick' ) ) 
checkTrue( file.exists( 'Protein_1/report/Protein_1_phylo_fun_tree.png' ) ) 

# clean up:
unlink( "Protein_1", recursive=T )
