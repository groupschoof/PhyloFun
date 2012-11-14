library(tools)
library(RUnit)
library(rredis)

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
src.project.file( 'src','domainArchitectureSimilarityRedis.R' )
src.project.file( 'src','domainArchitectureSimilarity.R' )

redis.con <- try( redisConnect() )
if( identical( class( redis.con ), "try-error" ) )
  stop( "Could not connect to redis. Did you start it on localhost?" )

# Empty redis for testing:
redisFlushAll()

# Test initializeDomainWeights
dwt <- read.table( project.file.path( "test", "domain_weights_database_test.tbl" ) )
print("Testing initializeDomainWeights(...)")
no.res <- initializeDomainWeights( dwt )
checkEquals( length( redisKeys() ), nrow( dwt ) )
dom.id <- rownames(dwt)[[1]] 
checkEquals( redisGet( dom.id ), dwt[[ dom.id, 1 ]] )

# Test initializeDomainAnnotations
print("Testing initializeDomainAnnotations(...)")
redisFlushAll()
f <- file( project.file.path( "test", "test_annotations_2.tbl" ), "r" )
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

# Test pairwiseDASDistanceKey
print("Testing pairwiseDASDistanceKey(...)")
prot.acc.d <- "goat"
prot.acc.e <- "sheep"
key.exp <- "goat_sheep"
checkEquals( pairwiseDASDistanceKey(prot.acc.d, prot.acc.e), key.exp )
checkEquals( pairwiseDASDistanceKey(prot.acc.e, prot.acc.d), key.exp )

# Test pairwiseDomainArchitectureDistanceRedis
print("Testing pairwiseDomainArchitectureDistanceRedis(...)")
redisFlushAll()
dom.annos.2 <- do.call("cbind",
  list(
    "BAAA"=list( "InterPro"=list("IPR001957", "IPR010921", "IPR018312") ),
    "MAAA"=list( "InterPro"=list("IPR001957", "IPR010921", "IPR018312", "IPR024633", "IPR013159", "IPR003593"))
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
checkEquals( round( redisGet( pairwiseDASDistanceKey("BAAA", "MAAA") ), 7 ),
  0.0304956 )

# Clean up, boy:
redisFlushAll()
