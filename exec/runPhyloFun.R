require( PhyloFun )

# Hail User:
print( paste(
  "Usage: Rscript runPhyloFun.R -q path/2/query_proteins.fasta ( -p path/2/phmmer_results.tbl OR -b path/2/blast_results.tbl )",
  "[ -c cores_to_use (default all) ] [ -f FastTree[MP] (default FastTreeMP) ]",
  "[ -e add.evidence.codes ( example: '-e TAS,IC' - Default all experimentally verified. Set to 'ALL', if no filtering for evidence-codes is wanted. ) ]",
  "[ -n n.best.hits Maximum number of best scoring results from sequence similarity search to use for each query protein. (default 1000) ]",
  "[ -h true ( Write out a statistics table 'homologs_stats.txt' for each Query Protein's set of homologs: Set-size and bit.score distribution. ) ]",
  "[ -m true ( Write out a statistics table 'msa_stats.txt' for each Query Protein's multiple sequence alignments: Differences in number of sequences and positions betwen original and filtered MSAs. ) ]",
  "[ -r true ( Generate an HTML report for each Query Protein. Will be stored in an extra folder 'report'. Default: false ) ]"
) )
print( '' )
print(
  paste( "WARNING: The PhyloFun pipeline uses other programs to generate multiple sequence alignments (MAFFT),",
    "filter them for conserved regions (GBlocks), and generate a phylogenetic tree of the MSA (FastTree[MP]).",
    "These programs need to be in your $PATH and require protein accessions of your PHMMER or Blast homolgy searches to be Uniprot accessions.",
    "Finally the accessions of your query proteins should consist only of the following character class [a-zA-Z0-9_-]"
  )
)

# Input
phylo.fun.args <- commandLineArguments( commandArgs(trailingOnly = TRUE), list( 'c'=detectCores(), 'f'='FastTreeMP', 'n'=1000 ) )

# Evidence codes of GO annotations to accept:
go.anno.evdnc.cds <- if ( is.null( phylo.fun.args[[ 'e' ]] ) ) {
    EVIDENCE.CODES
  } else if ( length( phylo.fun.args[[ 'e' ]] ) &&
    phylo.fun.args[[ 'e' ]] == 'ALL' ) {
    'ALL'
  } else {
    c( EVIDENCE.CODES, str_split( phylo.fun.args[[ 'e' ]], ',' )[[ 1 ]] )
  }

# Read fasta:
aa.seqs <- sapply( readAAStringSet( phylo.fun.args[[ 'q' ]] ), function(s) replaceSelenocystein( toString(s) ) )
print( paste("Read", length(aa.seqs), "sequences from", phylo.fun.args[[ 'q' ]] ) )

# Parse sequence similarity search results:
seq.search.rslts <- if ( ! is.null( phylo.fun.args[[ 'b' ]] ) ) {
  parseBlastTable( read.table( phylo.fun.args[[ 'b' ]] ) )
} else if ( ! is.null( phylo.fun.args[[ 'p' ]] ) ) {
  parsePhmmerTable( 
    scan( file=phylo.fun.args[[ 'p' ]], what=character(), sep="\n" )
  )
}
print( paste( "Parsed sequence similarity search results table. Got", nrow(seq.search.rslts), "query-hit pairs" ) )

# Sanitize protein accessions:
accs <- unlist( setNames( lapply( names(aa.seqs), sanitizeUniprotAccession ), names(aa.seqs) ) )
print( "Parsed the query proteins' accessions using function sanitizeUniprotAccession(...) ." )
print( 
  paste( "WARNING: If your query accessions are not matching, PhyloFun will fail to find their",
  "accessions in the sequence similarity search results, nor will Gblocks accept such sequence names!",
  "See function sanitizeUniprotAccession(…) for details." )
)

# For each query protein, do:
for ( prot.acc in accs ) {
  homologs <- sanitizeUniprotAccessions(
    bestHits( seq.search.rslts, prot.acc, n.best.hits=phylo.fun.args[[ 'n' ]] )
  )

  if ( nrow( homologs ) > 0 ) {
    orig.acc <- names( accs[ accs[] == prot.acc ] )
    if ( ! file.exists( prot.acc ) )
      dir.create( prot.acc )

    # Log statistics of the sequence homologs, if requested:
    if ( ! is.null( phylo.fun.args[[ 'h' ]] ) ) {
      write.table( homologsStats( homologs, prot.acc ),
        file=paste( prot.acc, '/homologs_stats.txt', sep='' )
      )
    }

    # Generate multiple sequence alignment ( MSA ) using MAFFT:
    print( "Generating multiple sequence alignment (MSA)" )
    acc.hmlgs.file <- paste( prot.acc, '/homologs.fasta', sep='' )
    hit.accs <- as.character( homologs[ , 'hit.name' ] )
    # Obtain AA sequences from Uniprot's Web Service:
    downloadSequences( hit.accs, acc.hmlgs.file )   
    # Replace all possible occurrences of u or U with x and X respectively:
    replaceSelenocysteinInFasta( acc.hmlgs.file )   
    # Read in the resulting filtered homologous AA Sequences:
    hit.seqs <- readAAStringSet( acc.hmlgs.file )   
    upr.accs <- names( hit.seqs )
    # Append Query's AA-Sequence and use sanitized original name:
    acc.hmlgs <- setNames( c( hit.seqs, AAStringSet( aa.seqs[ orig.acc ] ) ),
      c( upr.accs, prot.acc ) )
    writeXStringSet( acc.hmlgs, file=acc.hmlgs.file )
    # Generate Multiple Sequence Alignment:
    acc.msa.file <- paste( prot.acc, "/msa.fasta", sep="" )
    system( paste( "mafft --auto", acc.hmlgs.file, ">", acc.msa.file ) )

    # Remove duplicated accessions from MSA:
    print( "Removing duplicated accessions -if existing- from MSA" )
    uniqueHomologs( acc.msa.file )

    # Filter the MSA for highly conserved regions using GBlocks:
    print( "Filtering MSA for highly conserved regions" )
    acc.gblocks.msa.file <- paste( acc.msa.file, '-gb', sep='' )
    system( paste( 'Gblocks', acc.msa.file, '-b5=h -t=p -p=n' ) )

    # Phylo-Filter the GBlocks result discarding empty sequences or those
    # consisting of too many gap characters:
    acc.gblocks.msa <- readAAStringSet( acc.gblocks.msa.file ) 
    acc.msa.phylo.filtered <- filterMultipleSequenceAlignment( acc.gblocks.msa )
    # Should we use the original unfiltered MSA or the MSA filtered in two
    # steps (Gblocks followed by phyloFun's filter)?
    acc.msa.chosen.file <- if ( chooseFilteredAlignment(
        readAAStringSet( acc.msa.file ),
        acc.msa.phylo.filtered
      ) ) {
      # Only write out PhyloFun filtered MSA, if it is different from the
      # Gblocks' filtered one:
      if ( ! msaEqual( acc.gblocks.msa, acc.msa.phylo.filtered ) ) {
        acc.msa.phylo.filtered.file <- paste( acc.gblocks.msa.file, '-phylo_filtered', sep='' )
        writeXStringSet( acc.msa.phylo.filtered, acc.msa.phylo.filtered.file )
        acc.msa.phylo.filtered.file
      } else {
        acc.gblocks.msa.file
      }
    } else {
      acc.msa.file
    }
    # Tell the User, which MSA is going to be used:
    print( paste( 'MSA has been filtered with Gblocks and PhyloFun. Going to use',
      acc.msa.chosen.file, 'for the phylogenetic reconstruction' ) )
    print( '' )
    # Report MSA statistics, if requested:
    if ( ! is.null( phylo.fun.args[[ 'm' ]] ) ) {
      write.table( 
        msaStats( readAAStringSet( acc.msa.file ),
          readAAStringSet( acc.msa.chosen.file ),
          prot.acc 
        ),
        file=paste( prot.acc, '/msa_stats.txt', sep='' )
      )
    }

    # Construct the phylogenetic max likelihood tree of the above alignment
    # using FastTree[MP]:
    print( "Constructing maximum likelihood phylogenetic tree" )
    acc.phyl.tree.file <- paste( prot.acc, '/ml_tree.newick', sep='' )
    system( paste( phylo.fun.args[[ 'f' ]], '<', acc.msa.chosen.file, '>', acc.phyl.tree.file ) ) 

    # Compute probability distributions for GO terms of the three different
    # types 'biological_process', 'cellular_component', and
    # 'molecular_function':
    print( "Starting PhyloFun on phylogenetic tree" )
    print( "Note: Branch Lengths of phylogenetic tree will be rounded to two decimal digits!" )
    acc.phyl.tree <- roundBranchLengths( read.tree( acc.phyl.tree.file ) )
    acc.hmlgs.annos <- if ( length( go.anno.evdnc.cds ) == 1 &&
      go.anno.evdnc.cds == 'ALL' ) {
      retrieveUniprotAnnotations( hit.accs )
    } else {
      retrieveExperimentallyVerifiedGOAnnotations(
        hit.accs, evidence.codes=go.anno.evdnc.cds )
    }    
    
    if ( ! is.null( acc.hmlgs.annos ) && ncol( acc.hmlgs.annos ) > 0 ) {
      go.con <- connectToGeneOntology()
      acc.go.type.annos  <- goTypeAnnotationMatrices( acc.hmlgs.annos,
        go.con=go.con
      )
      dbDisconnect( go.con )
      acc.go.anno.spaces <- goAnnotationSpaceList( acc.go.type.annos )
      quoted.acc <- surroundEachWithQuotes( prot.acc )
  
      go.types <- c( 'biological_process', 'cellular_component',
        'molecular_function'
      )
      acc.go.predictions <- setNames(
        lapply( go.types, function( go.type ) {
  
          print( paste( 'Computing function predictions for GO type', go.type ) )
  
          acc.bayes.evdnc <- annotationMatrixForBayesNetwork(
            acc.go.type.annos[[ go.type ]], 
            all.accessions=setdiff( acc.phyl.tree$tip.label, prot.acc )
          )
          if ( ! is.null( acc.bayes.evdnc ) ) {
            acc.bayes.netw <- grain( compileCPT(
              bayesNodes( acc.phyl.tree, acc.go.anno.spaces[[ go.type ]] )
            ) )
            predict.grain( acc.bayes.netw, response=quoted.acc,
              newdata=acc.bayes.evdnc, type='dist'
            )
          } else {
            NULL
          }
        }),
        go.types
      )
      # Finished predicting GO term annotations.
      
      # Write out complete results:
      f <- file( paste( prot.acc, "/phyloFun_R_serialized.txt", sep="" ), "w" )
      serialize( acc.go.predictions, f )
      close( f )
  
      # Human readable results:
      acc.pf.predictions <- goTermPredictionTable( acc.go.predictions, quoted.acc )
      write.table( acc.pf.predictions,
        file=paste( prot.acc, '/go_term_predictions.tbl', sep='' ),
        row.names=F
      )

      # If requested generate an HTML report:
      if ( ! is.null( phylo.fun.args[[ 'r' ]] ) ) {
        print( "Generating HTML report" )
        report.dir <- paste( prot.acc, '/report', sep='' )
        report.tree.fn <- paste( prot.acc, '_phylo_fun_tree', sep='' ) 
        report.tree.path <- paste( report.dir, '/', report.tree.fn, sep='' )
        if ( ! file.exists( report.dir ) )
          dir.create( report.dir, recursive=T )

        plot.rslt <- plotPhyloFunTree( prot.acc, acc.phyl.tree, acc.pf.predictions, acc.go.type.annos,
          paste( report.tree.path, '.png', sep='' ) )
        htmlReport( paste( report.tree.fn, '.png', sep='' ), prot.acc, plot.rslt$caption,
          paste( report.dir, '/', prot.acc, '_report.html', sep='' ) )
        write.tree( plot.rslt$tree.with.abbrevs, file=paste( report.tree.path, '.newick', sep='' ) )
      }
    } else {
      # No PhyloFun annotations can be assigned, because there are no GO term
      # annotations available for the Query's homologs:
      msg <- c(
        "WARNING: There were no GO term annotations available for the found homologous sequences of significant similarity.",
        paste( prot.acc, "will therefore have NO PhyloFun predictions!" ) 
      )
      print( msg )
      f <- file( paste( prot.acc, '/go_term_predictions.tbl', sep='' ), 'w' )
      writeLines( msg, con=f )
      close( f )
    }

    print( paste( "Finished computations for", orig.acc ) )
  } else {
    print( paste(
      "Warning: No sequence homologs were found for query protein",
      prot.acc ) )
  }
}

print( "DONE" )
