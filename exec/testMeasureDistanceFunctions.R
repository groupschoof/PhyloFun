require( PhyloFun )

# Test measureDistances
print("Testing measureDistances(...)")
f <- file( project.file.path("test", "test_annotations.tbl"), "r" )
annos <- unserialize(f)
close(f)
aa.seqs <- sapply(
  readAAStringSet( project.file.path("test", "measureDistanceFunctions.fasta") ),
  function(s) toString(s)
)
blast.rslt.tbl <- matrix(
  c( 
    "A0RLX8", "A0PKB2",
    "A0RLX8", "A0K2M8",
    "A0K2M8", "A0KR35",
    "A0Q3U6", "A0KEC3"
  ),
  ncol=2, byrow=T
)
domain.weights.table <- read.table( project.file.path( "test",
    "domain_weights_database.tbl" )
)
go.term <- "GO:0017111"
dists <- measureDistances( go.term, annos, blast.rslt.tbl, aa.seqs,
  domain.weights.table )
# print( dists )
checkTrue( ! is.null( dists ) )
checkEquals( class(dists), 'matrix' )
checkEquals( nrow(dists), 4 )
exp.dists <- matrix(
  c(
    1.13, 0.56, T, 1.26,
    1.13, 0.00, T, 1.13,
    0.67, 0.32, 0, 0.74,
    0.79, 1.00, 0, 1.27
  ),
  byrow=T, nrow=4,
  dimnames=list(
    c( "A0PKB2_A0RLX8", "A0K2M8_A0RLX8", "A0K2M8_A0KR35", "A0KEC3_A0Q3U6" ),
    c( "Sequence.Distance", "Domain.Architecture.Distance", "Share.GO:0017111", "Euclidean.Distance.To.Origin" )
  )
)
# print( exp.dists )
checkEquals( round( dists, 2 ), exp.dists )
no.dists <- measureDistances( "GO:0005737", annos, blast.rslt.tbl, aa.seqs,
  domain.weights.table )
checkTrue( is.null( no.dists ) )

# Test pMutation
print("Testing pMutation(...)")
checkEquals( pMutation( 0, 0 ), 0 )
checkEquals( pMutation( 3, 0 ), 0 )
checkEquals( pMutation( 0, 3 ), 1 )
checkEquals( pMutation( 3, 1, 0.5 ), 0.5 )

# Test mutationProbabilityDistribution
dist.tbl <- read.table( stringsAsFactors=FALSE, text=
"A B 0.1
A C 0.5
A D 1.0")
annot.tbl <- read.table( stringsAsFactors=FALSE, text=
"GO:0001234 IEA A
GO:0001234 IEA B 
GO:0001234 IEA D")
res.mutationProbabilityDistribution <- mutationProbabilityDistribution(
  dist.tbl, annot.tbl, 'GO:0001234' )


# Test gridPMutation
print("Testing gridPMutation(...)")
p.mut.quants <- gridPMutation( p.mut.das.seq )
exp.p.mut.quants.1 <- matrix( c( 0.0, NA, NA, 0.33, NA, NA, NA, NA, NA, NA ),
  ncol=1, dimnames=list(
    c(),
    "p.mutation|Euclidean.Distance.To.Origin"
  )
)
checkEquals( p.mut.quants[ , 2, drop=F ], exp.p.mut.quants.1 )

# Test pMutationMinMaxParentValues
print("Testing pMutationMinMaxParentValues(...)")
go_0036057.p.mut.dists <- matrix(
  c(0.0, 0.1, 0.1, 0.4, 0.4, 0.6, 0.7, 1.0,
    0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.9, 1.2,
    0.0, 0.95, 0.1, 0.9, 0.15, 0.8, 0.25, 1.0,
    1, 1, 1, 1, 0, 0, 0, 0,
    0.5, 0.75, 0.7, 1.0, 1.25, 1.15, 1.5, 1.75
  ), nrow=8, ncol=5,
  dimnames=list(
    c(),c("p.mutation|Sequence.Distance", "Sequence.Distance",
      "Domain.Architecture.Distance", "Share.GO:0036057",
      "Euclidean.Distance.To.Origin")
  )
)
go_0036057.p.mut.min.max <- pMutationMinMaxParentValues( go_0036057.p.mut.dists,
  "p.mutation|Sequence.Distance" )
exp.go_0036057.p.mut.min.max <- matrix( 
  c( 0.0, 0.5, 0.00, 0.50, 0.5, 0.00, 0.50,
     0.1, 0.6, 0.10, 0.70, 0.6, 0.95, 0.75,
     0.4, 0.7, 0.15, 1.00, 0.7, 0.90, 1.25,
     0.6, 0.8, 0.80, 1.15, 0.8, 0.80, 1.15,
     0.7, 0.9, 0.25, 1.50, 0.9, 0.25, 1.50,
     1.0, 1.2, 1.00, 1.75, 1.2, 1.00, 1.75
  ), ncol=7, byrow=T, dimnames=list( c(),
  c( "p.mutation|Sequence.Distance", "min.Sequence.Distance",
      "min.Domain.Architecture.Distance", "min.Euclidean.Distance.To.Origin",
      "max.Sequence.Distance", "max.Domain.Architecture.Distance",
      "max.Euclidean.Distance.To.Origin" )
  )
)
checkEquals( go_0036057.p.mut.min.max, exp.go_0036057.p.mut.min.max )
