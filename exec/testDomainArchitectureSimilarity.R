require( PhyloFun )

# Load domain.weights.table:
dwt <- read.table(
  project.file.path( "data", "uniprotkb_based_domain_weights.tbl" )
  )[ , 'DW', drop=F ]

# Load test annotations:
f <- file( project.file.path(  "test_annotations_2.tbl" ), "r" )
am <- unserialize( f )
close( f )

# Test constructVectorSpaceModel
print( "Testing constructVectorSpaceModel(...)" )
expc.uniq.annos <- sort(
  c("IPR003593", "IPR001957", "IPR020591", "IPR018312",
    "IPR013159", "IPR013317", "IPR024633", "IPR010921")
)
vsm <- constructVectorSpaceModel( am )
checkEquals( vsm, expc.uniq.annos )

# Test generateDomainArchitectureSpaceVectors
print("Testing generateDomainArchitectureSpaceVectors(...)")
# Initialize domain weights database:
dwd <- read.table( project.file.path(  "domain_weights_database.tbl" ),
  header=T )
# Initialize expected result:
f <- file( project.file.path(  "testDomainArchitectureSpaceVectors.txt"),
  "r" )
exp.dasv <- unserialize( f )
close(f)

dasv <- generateDomainArchitectureSpaceVectors( vsm, am, dwd )
# print( dasv )
checkEquals( exp.dasv, dasv )

# Test partial generateDomainArchitectureSpaceVectors
print("Testing partial generateDomainArchitectureSpaceVectors(...)")
test.accs <- c('A0RLX8', 'Protein_1')
part.dasv <- generateDomainArchitectureSpaceVectors( vsm, am, dwd,
  vectors.4.accessions=test.accs )
# print( part.dasv )
checkEquals( exp.dasv[ , colnames(exp.dasv) %in% test.accs ], part.dasv )

# Test pairwiseDomainArchitectureDistance
print("Testing pairwiseDomainArchitectureDistance(...)")
checkEquals(
  round( pairwiseDomainArchitectureDistance(c(1.0,1.0), c(1.0,0.0)), 6 ),
  ( 1.0 - 0.707107 )
)
checkEquals(
  round( pairwiseDomainArchitectureDistance(c(1.0,0.0), c(0.0,1.0)), 6 ),
  1.0
)
checkEquals(
  round( pairwiseDomainArchitectureDistance(c(1.0,0.0), c(2.0,0.0)), 6 ),
  0.0
)
checkEquals(
  pairwiseDomainArchitectureDistance(c(0.0, 2.0, 3.0, 0.0), c(0.0, 2.0, 3.0, 0.0)),
  0.0
)

# Test domainArchitectureDistances
print("Testing domainArchitectureDistances(...)")
# Construct a matrix of test domain architecture space vectors:
dsvs <- matrix(
  c(1.0, 2.0, 3.0, 4.0,  0.0, 2.0, 3.0, 0.0, 
    1.0, 0.0, 0.0, 4.0,  0.0, 0.0, 3.0, 4.0),
  byrow=F, nrow=4,
  dimnames=list( 
    c("IPR000111", "IPR000222", "IPR000333", "IPR000444"),
    c("Protein_A", "Protein_B", "Protein_C", "Protein_D")
  )
)
# expected result:
exp.rslt <- matrix( 
  c(0.00000000, 0.34171941, 0.24722735, 0.08712907,
    0.34171941, 0.00000000, 1.00000000, 0.50076982,
    0.24722735, 1.00000000, 0.00000000, 0.22388600,
    0.08712907, 0.50076982, 0.22388600, 0.00000000),
  nrow=4, ncol=4,
  dimnames=list(
    c('Protein_A', 'Protein_B', 'Protein_C', 'Protein_D'),
    c('Protein_A', 'Protein_B', 'Protein_C', 'Protein_D'))
)
#            Protein_A Protein_B Protein_C  Protein_D
# Protein_A 0.00000000 0.3417194 0.2472273 0.08712907
# Protein_B 0.34171941 0.0000000 1.0000000 0.50076982
# Protein_C 0.24722735 1.0000000 0.0000000 0.22388600
# Protein_D 0.08712907 0.5007698 0.2238860 0.00000000
rslt <- domainArchitectureDistances( dsvs )
checkEquals( as.matrix(rslt), exp.rslt )

# Test partialDomainArchitectureDistances
print("Testing partialDomainArchitectureDistances(...)")
test.prot.pairs <- list( c( 'A0RLX8', 'A0PKB2' ), c( 'A0PKB2', 'Protein_1' ) )
part.das.dists <- partialDomainArchitectureDistances( am, dwd, test.prot.pairs )
# print( part.das.dists )
checkTrue( ! is.na( part.das.dists ) && ! is.null( part.das.dists ) )
checkEquals( class( part.das.dists ), "list" )
checkEquals( length( part.das.dists ), length( test.prot.pairs ) )
checkEquals( round( part.das.dists[[ 1 ]], 6 ), 0.322616 )

# Test pairwiseSequenceDistance
print("Testing pairwiseSequenceDistance(...)")
p1 <- "LALDTKQIWFTTLGTLQNQILRYDYDTWLKTTALVSVANDLAVIGAPNVTTKQVIEDRFMSVLRRALGEVLGYQVNVRVIISSATPAPSEPVAVTPSEPSPTTEVAEPSFASFNQAAPMLNQLPLGDPNRSSVLNPRYTFSSFIVGTSNRLAHAACMAVAEHPAQAYNPLFLYGGVGLGKTHLLQAIGNYALDRNPEVNVLYVSSEKFTNDLINAIRRQQTEEFRIRYRNIDILLIDDIQFIAGKEGTQEEFFHTFNTLHGAGKQIVLSSDRPPKAILTLEERLRSRFEWGLIVDVQNPDLETRTAILRAKGETLQVPVSSEVIDFLAQRIQSNIRELEGCLNRVIAYANLNRTPVTVEVASAALADLLDTSRRKRVTADDIFREVSQHYGIDQRAIRGRGRSRNVVLPRQVVMYLLREETDASLVEIGELLGGRDHTTVMHGYNKITDDLTSDARLRNDITSLRQRLYGENAR"
p2 <- "MQSIEDIWQETLQIVKKNMSKPSYDTWMKSTTAHSLEGNTFIISAPNNFVRDWLEKSYTQFIANILQEITGRLFDVRFIDGEQEENFEYTVIKPNPALDEDGIEIGKHMLNPRYVFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAVGHYVQQHKDNAKVMYLSSEKFTNEFISSIRDNKTEEFRTKYRNVDVLLIDDIQFLAGKEGTQEEFFHTFNTLYDEQKQIIISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKADGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLVNKDITAGLAAEALKDIIPSSKSQVITISGIQETVGEYFHVRLEDFKAKKRTKSIAFPRQIAMYLSRELTDASLPKIGDEFGGRDHTTVIHAHEKISQLLKTDQVLKNDLAEIEKNLRKSQNMF"
p.dist <- round( pairwiseSequenceDistance( p1, p2 ), 6 )
checkTrue( ! is.null(p.dist) )
checkEquals( class(p.dist), 'numeric' )
checkEquals( p.dist, 0.813185 )
# Assert symmetric distance measure:
checkEquals( round( pairwiseSequenceDistance( p2, p1 ), 6 ), 0.813185 )

# Test replaceSelenocystein
print("Testing replaceSelenocystein(...)")
p3 <- "MTASLWQQCLNRLQDELPSAEFSMWIRPLQAELSDNTLTLYAPNRFVLDWVRDKYLIRVNGIINELCGVDGPTLRFDIGNRPHPVAVARAPARGADPVNNSQKSWESKAEAKPEPNHKSNTNVNYTFENFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNAIKERKQDAKVIYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMRKADENQIHLPDEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAINIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKLADLLSKRRSRSVARPRQLAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLKEESHDIKEDYSNLIRTLSS"
p3.san <- "MTASLWQQCLNRLQDELPSAEFSMWIRPLQAELSDNTLTLYAPNRFVLDWVRDKYLIRVNGIINELCGVDGPTLRFDIGNRPHPVAVARAPARGADPVNNSQKSWESKAEAKPEPNHKSNTNVNYTFENFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNAIKERKQDAKVIYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMRKADENQIHLPDEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAINIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKLADLLSKRRSRSVARPRQLAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLKEESHDIKEDYSNLIRTLSS"
checkEquals( replaceSelenocystein( p3 ), p3.san )
# Check on set of AA seqs:
aa.seqs <- readAAStringSet( project.file.path( 
  'non_unique_hmlgs.fasta' ) )
no.u.aa.seqs <- lapply( aa.seqs, replaceSelenocystein )
exp.aa.seq <- 'MAYNHDGMPDQQFLWDVFQRVDKDRSGHISADELQVALSNGTWSAFNPETIRLMIGMFDRENKGTVSFKDFGALWKYVTDWQNCFRSFDRDNSGNIDKTELKTALTSFGYRLSDHLIDVLLRKFDRFGRGTILFDDFIQCCIVLYTLTTAFRQHDTDLDGIITIHYEQFLXSMVFSLKI'
checkEquals( exp.aa.seq, no.u.aa.seqs[[ 'Q8UUX9' ]] )

# Test sequenceDistances
print("Testing sequenceDistances(...)")
prot.lst <- list( "P1"=p1, "P2"=p2, "P3"=p3 )
prot.seq.dists <- sequenceDistances( prot.lst )
checkEquals( class(prot.seq.dists), 'dist' )
prot.seq.dists.mtrx <- round( as.matrix( prot.seq.dists ), 6 )
checkEquals( prot.seq.dists.mtrx[[ 1, 1 ]], 0.0 )
checkEquals( prot.seq.dists.mtrx[[ 1, 2 ]], 0.813185 )
print( prot.seq.dists.mtrx[[ 1, 3 ]] )
checkEquals( prot.seq.dists.mtrx[[ 1, 3 ]], 0.925284 )

# Test partialSequenceDistances
print("Testing partialSequenceDistances(...)")
prot.pairs <- list( c( "P1", "P2" ) )
prot.seq.dists <- partialSequenceDistances( prot.lst, prot.pairs )
# print( prot.seq.dists )
checkEquals( class( prot.seq.dists ), 'list' )
checkEquals( length( prot.seq.dists ), 1 )
checkEquals( round( prot.seq.dists[[ 1 ]], 6 ), 0.813185 )

# Test distanceMatrixIndices
print("Testing distanceMatrixIndices(...)")
accs <- c("A", "B", "C", "D", "E", "F", "G", "H")
exp.inds <- list( c("A","B"), c("A","C"), c("A","D"), c("A","E"), c("A","F"),
  c("A","G"), c("A","H"), c("B","C"), c("B","D"), c("B","E"), c("B","F"),
  c("B","G"), c("B","H"), c("C","D"), c("C","E"), c("C","F"), c("C","G"),
  c("C","H"), c("D","E"), c("D","F"), c("D","G"), c("D","H"), c("E","F"),
  c("E","G"), c("E","H"), c("F","G"), c("F","H"), c("G","H")
)
checkEquals( exp.inds, distanceMatrixIndices( accs ) )

# Test distanceIndices
print("Testing distanceIndices(...)")
batch.size <- 3
dist.ind.1 <- distanceIndices( 1, batch.size, accs )
checkEquals( dist.ind.1, exp.inds[1:3] )

dist.ind.2 <- distanceIndices( 2, batch.size, accs )
checkEquals( dist.ind.2, exp.inds[4:6] )

dist.ind.3 <- distanceIndices( 3, batch.size, accs )
checkEquals( dist.ind.3, exp.inds[7:9] )

dist.ind.4 <- distanceIndices( 4, batch.size, accs )
checkEquals( dist.ind.4, exp.inds[10:12] )

# Test edge cases of function distanceIndices(...)
dist.ind.5 <- distanceIndices( 1, 9, accs )
checkEquals( dist.ind.5, exp.inds[1:9] )

dist.ind.6 <- distanceIndices( 2, 8, accs )
checkEquals( dist.ind.6, exp.inds[9:16] )

dist.ind.7 <- distanceIndices( 4, 8, accs )
checkEquals( dist.ind.7, exp.inds[25:28] )

dist.ind.8 <- distanceIndices( 4, 6, accs )
checkEquals( dist.ind.8, exp.inds[19:24] )

batch.sizes <- 1:29
for ( batch.size in batch.sizes ) {
  all.batches <- unlist(
    lapply( 1:ceiling( 28 / batch.size ),
      function( batch.no ) {
        distanceIndices( batch.no, batch.size, accs )
      }
    ),
    recursive=F
  )
  print( checkEquals( all.batches, exp.inds ) )
}

# Test pairsFromBlastResult
print("Testing pairsFromBlastResult(...)")
blast.result <- matrix(
  c("Query_A", "Query_B", "Query_C", "Hit_A", "Hit_B", "Hit_C" ),
  nrow=3, ncol=2
)
blast.pairs <- pairsFromBlastResult( blast.result )
exp.pairs <- list( c("Query_A", "Hit_A"), c("Query_B", "Hit_B"), c("Query_C", "Hit_C") )
checkEquals( blast.pairs, exp.pairs )

# Test uniquePairs
print("Testing uniquePairs(...)")
pairs.test <- matrix( c( "A", "B", "B", "A", "C", "D", "D", "C", "X", "X" ), byrow=T, nrow=5 )
uniq.pairs.test <- uniquePairs( pairs.test )
checkEquals( uniq.pairs.test, matrix( c( "A", "B", "C", "D" ), byrow=T, nrow=2 ) )
pairs.test <- matrix( character() )
uniq.pairs.test <- uniquePairs( pairs.test )
checkTrue( is.null( uniq.pairs.test ) )
pairs.test <- matrix( character(), nrow=5 )
uniq.pairs.test <- uniquePairs( pairs.test )
checkTrue( is.null( uniq.pairs.test ) )
pairs.test <- matrix( character(), ncol=5 )
uniq.pairs.test <- uniquePairs( pairs.test )
checkTrue( is.null( uniq.pairs.test ) )
