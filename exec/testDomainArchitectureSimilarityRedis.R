require( PhyloFun )

redis.con <- try( redisConnect() )
if( identical( class( redis.con ), "try-error" ) )
  stop( "Could not connect to redis. Did you start it on localhost?" )

# Empty redis for testing:
redisFlushAll()

# Test initializeDomainWeights
dwt <- read.table( project.file.path(  "domain_weights_database_test.tbl" ) )
print("Testing initializeDomainWeights(...)")
no.res <- initializeDomainWeights( dwt )
checkEquals( length( redisKeys() ), nrow( dwt ) )
dom.id <- rownames(dwt)[[1]] 
checkEquals( redisGet( dom.id ), dwt[[ dom.id, 1 ]] )

# Test initializeDomainAnnotations
print("Testing initializeDomainAnnotations(...)")
redisFlushAll()
f <- file( project.file.path(  "test_annotations_2.tbl" ), "r" )
dom.annos <- unserialize( f )
close( f )
no.res <- initializeDomainAnnotations( dom.annos )
checkEquals( length( redisKeys() ), ncol( dom.annos ) )
prot.acc <- colnames( dom.annos )[[1]]
checkEquals( sort( unlist( redisSMembers( prot.acc ) ) ),
  sort( dom.annos[[ "InterPro", prot.acc ]] )
)

# Test pairwiseVectorSpaceModel
print("Testing pairwiseVectorSpaceModel(...)")
prot.acc.a <- colnames( dom.annos )[[1]]
prot.acc.b <- colnames( dom.annos )[[2]]
vsm <- pairwiseVectorSpaceModel( prot.acc.a, prot.acc.b )
vsm.exp <- sort( union( dom.annos[[ "InterPro", prot.acc.a ]],
    dom.annos[[ "InterPro", prot.acc.b ]] ) )
checkEquals( sort( unlist( vsm ) ), vsm.exp )

# Test generateDomainArchitectureSpaceVectorRedis
print("Testing generateDomainArchitectureSpaceVectorRedis(...)")
no.res <- initializeDomainWeights( dwt ) # Initialize domain weights
prot.acc.c <- "SHEEP"
vsm.2 <- c("IPR003350", "IPR022724", "IPR025668", "IPR004305", "IPR727272")
for( d in vsm.2 ) {
  redisSAdd( prot.acc.c, charToRaw(d) )
}
das.vec.c <- generateDomainArchitectureSpaceVectorRedis( prot.acc.c, vsm.2 )
das.vec.c.exp <- c( 3.32700478081248, 14.9922872943264,
  7.19508911404042, 4.26253229157448, 0.0 )
checkEquals( das.vec.c, das.vec.c.exp )

# Test pairwiseDistanceKey
print("Testing pairwiseDistanceKey(...)")
prot.acc.d <- "goat"
prot.acc.e <- "sheep"
key.exp <- "goat_sheep_das_dist"
checkEquals( pairwiseDistanceKey(prot.acc.d, prot.acc.e), key.exp )
checkEquals( pairwiseDistanceKey(prot.acc.e, prot.acc.d), key.exp )
checkEquals( pairwiseDistanceKey(prot.acc.d, prot.acc.e, distance.type="seq_dist"),
  paste( prot.acc.d, prot.acc.e, "seq_dist", sep="_" )
)

# Test pairwiseDomainArchitectureDistanceRedis
print("Testing pairwiseDomainArchitectureDistanceRedis(...)")
redisFlushAll()
dom.annos.2 <- do.call("cbind",
  list(
    "BAAA"=list( "InterPro"=list("IPR001957", "IPR010921", "IPR018312") ),
    "MAAA"=list( "InterPro"=list("IPR001957", "IPR010921", "IPR018312", "IPR024633", "IPR013159", "IPR003593")),
    "MABA"=list( "InterPro"=list("IPR001957", "IPR018312") ),
    "BAMA"=list( "InterPro"=list("IPR001957", "IPR024633", "IPR013159", "IPR003593"))
  )
)
dom.wghts.2 <- matrix(
  c(3.327005, 14.992287, 7.195089, 4.262532, 0.456778),
  ncol=1, nrow=5,
  dimnames=list(
    c("IPR001957", "IPR010921", "IPR018312", "IPR024633", "IPR013159"),
    c("DW")
  )
)
# Initialize domain weights:
no.res <- initializeDomainWeights( dom.wghts.2 )
# Initialize protein domain annotations:
no.res <- initializeDomainAnnotations( dom.annos.2 )
das.dist.a.b <- pairwiseDomainArchitectureDistanceRedis( "BAAA", "MAAA" )
checkEquals( round(das.dist.a.b, 7), 0.0304956 )
checkEquals( round( redisGet( pairwiseDistanceKey("BAAA", "MAAA") ), 7 ),
  0.0304956 )

# Test partialDomainArchitectureDistancesRedis
print("Testing partialDomainArchitectureDistancesRedis(...)")
# Initialization done by last test.
no.res <- redisDelete( pairwiseDistanceKey( "BAAA", "MAAA" ) )
no.res <- partialDomainArchitectureDistancesRedis(
  colnames(dom.annos.2), colnames(dom.annos.2)[1:2]
)
# Check copied from previous test:
checkEquals( round( redisGet( pairwiseDistanceKey("BAAA", "MAAA") ), 7 ),
  0.0304956 )
checkTrue( is.null( redisGet(
  pairwiseDistanceKey( colnames(dom.annos.2)[[ 3 ]], colnames(dom.annos.2)[[ 4 ]] )
) ) )
# Test in parallel mode:
no.res <- redisDelete( pairwiseDistanceKey( "BAAA", "MAAA" ) )
init.thread.funk <- function() { redisConnect() }
close.thread.funk <- function() { redisClose() }
redisClose()
no.res <- partialDomainArchitectureDistancesRedis(
  colnames(dom.annos.2), colnames(dom.annos.2)[1:2],
  lapply.funk=mclapply, init.thread.funk=init.thread.funk,
  close.thread.funk=close.thread.funk, mc.preschedule=T,
  mc.cores=detectCores()
)
redisConnect()
# Check copied from previous test:
checkEquals( round( redisGet( pairwiseDistanceKey("BAAA", "MAAA") ), 7 ),
  0.0304956 )
checkTrue( is.null( redisGet(
  pairwiseDistanceKey( colnames(dom.annos.2)[[ 3 ]], colnames(dom.annos.2)[[ 4 ]] )
) ) )

# Test partialSequenceDistancesRedis
print("Testing partialSequenceDistancesRedis(...)")
redisFlushAll()
p1 <- "LALDTKQIWFTTLGTLQNQILRYDYDTWLKTTALVSVANDLAVIGAPNVTTKQVIEDRFMSVLRRALGEVLGYQVNVRVIISSATPAPSEPVAVTPSEPSPTTEVAEPSFASFNQAAPMLNQLPLGDPNRSSVLNPRYTFSSFIVGTSNRLAHAACMAVAEHPAQAYNPLFLYGGVGLGKTHLLQAIGNYALDRNPEVNVLYVSSEKFTNDLINAIRRQQTEEFRIRYRNIDILLIDDIQFIAGKEGTQEEFFHTFNTLHGAGKQIVLSSDRPPKAILTLEERLRSRFEWGLIVDVQNPDLETRTAILRAKGETLQVPVSSEVIDFLAQRIQSNIRELEGCLNRVIAYANLNRTPVTVEVASAALADLLDTSRRKRVTADDIFREVSQHYGIDQRAIRGRGRSRNVVLPRQVVMYLLREETDASLVEIGELLGGRDHTTVMHGYNKITDDLTSDARLRNDITSLRQRLYGENAR"
p2 <- "MQSIEDIWQETLQIVKKNMSKPSYDTWMKSTTAHSLEGNTFIISAPNNFVRDWLEKSYTQFIANILQEITGRLFDVRFIDGEQEENFEYTVIKPNPALDEDGIEIGKHMLNPRYVFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAVGHYVQQHKDNAKVMYLSSEKFTNEFISSIRDNKTEEFRTKYRNVDVLLIDDIQFLAGKEGTQEEFFHTFNTLYDEQKQIIISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKADGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLVNKDITAGLAAEALKDIIPSSKSQVITISGIQETVGEYFHVRLEDFKAKKRTKSIAFPRQIAMYLSRELTDASLPKIGDEFGGRDHTTVIHAHEKISQLLKTDQVLKNDLAEIEKNLRKSQNMF"
p3 <- "MQSIEDIWQETLQIVKKNMSKPSYDTWMKSTTAHSLEGNTFIISAPNNFVRDWLEKSYTQFIANILQEIT"
p4 <- "ILRYDYDTWLKTTALVSVANDLAVIGAPNVTTKQVIEDRFMSVLRRALGEVLGY"
p1.acc <- "Test_Protein_1"
p2.acc <- "Test_Protein_2"
p3.acc <- "Test_Protein_3"
p4.acc <- "Test_Protein_4"
aa.seqs <- setNames( list( p1, p2, p3, p4 ), c( p1.acc, p2.acc, p3.acc, p4.acc ) )
no.res <- partialSequenceDistancesRedis( aa.seqs, c( p1.acc, p2.acc ) )
p1.p2.dist.key <- pairwiseDistanceKey( p1.acc, p2.acc, distance.type="seq_dist" )
checkTrue( ! is.null( redisGet( p1.p2.dist.key ) ) )
checkTrue( is.null( redisGet( pairwiseDistanceKey( p3.acc, p4.acc, distance.type="seq_dist" ) ) ) )
checkEquals( pairwiseSequenceDistance( p1, p2 ), redisGet( p1.p2.dist.key ) )
# Check in parallel mode:
redisFlushAll()
p1 <- "LALDTKQIWFTTLGTLQNQILRYDYDTWLKTTALVSVANDLAVIGAPNVTTKQVIEDRFMSVLRRALGEVLGYQVNVRVIISSATPAPSEPVAVTPSEPSPTTEVAEPSFASFNQAAPMLNQLPLGDPNRSSVLNPRYTFSSFIVGTSNRLAHAACMAVAEHPAQAYNPLFLYGGVGLGKTHLLQAIGNYALDRNPEVNVLYVSSEKFTNDLINAIRRQQTEEFRIRYRNIDILLIDDIQFIAGKEGTQEEFFHTFNTLHGAGKQIVLSSDRPPKAILTLEERLRSRFEWGLIVDVQNPDLETRTAILRAKGETLQVPVSSEVIDFLAQRIQSNIRELEGCLNRVIAYANLNRTPVTVEVASAALADLLDTSRRKRVTADDIFREVSQHYGIDQRAIRGRGRSRNVVLPRQVVMYLLREETDASLVEIGELLGGRDHTTVMHGYNKITDDLTSDARLRNDITSLRQRLYGENAR"
p2 <- "MQSIEDIWQETLQIVKKNMSKPSYDTWMKSTTAHSLEGNTFIISAPNNFVRDWLEKSYTQFIANILQEITGRLFDVRFIDGEQEENFEYTVIKPNPALDEDGIEIGKHMLNPRYVFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAVGHYVQQHKDNAKVMYLSSEKFTNEFISSIRDNKTEEFRTKYRNVDVLLIDDIQFLAGKEGTQEEFFHTFNTLYDEQKQIIISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKADGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLVNKDITAGLAAEALKDIIPSSKSQVITISGIQETVGEYFHVRLEDFKAKKRTKSIAFPRQIAMYLSRELTDASLPKIGDEFGGRDHTTVIHAHEKISQLLKTDQVLKNDLAEIEKNLRKSQNMF"
p3 <- "MQSIEDIWQETLQIVKKNMSKPSYDTWMKSTTAHSLEGNTFIISAPNNFVRDWLEKSYTQFIANILQEIT"
p4 <- "ILRYDYDTWLKTTALVSVANDLAVIGAPNVTTKQVIEDRFMSVLRRALGEVLGY"
p1.acc <- "Test_Protein_1"
p2.acc <- "Test_Protein_2"
p3.acc <- "Test_Protein_3"
p4.acc <- "Test_Protein_4"
aa.seqs <- setNames( list( p1, p2, p3, p4 ), c( p1.acc, p2.acc, p3.acc, p4.acc ) )
no.res <- partialSequenceDistancesRedis( aa.seqs, c( p1.acc, p2.acc ), lapply.funk=mclapply,
  init.thread.funk=init.thread.funk, close.thread.funk=close.thread.funk,
  mc.cores=detectCores(), mc.preschedule=F
)
p1.p2.dist.key <- pairwiseDistanceKey( p1.acc, p2.acc, distance.type="seq_dist" )
checkTrue( ! is.null( redisGet( p1.p2.dist.key ) ) )
checkTrue( is.null( redisGet( pairwiseDistanceKey( p3.acc, p4.acc, distance.type="seq_dist" ) ) ) )
checkEquals( pairwiseSequenceDistance( p1, p2 ), redisGet( p1.p2.dist.key ) )

# Clean up, boy:
redisFlushAll()
