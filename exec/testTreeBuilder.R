require( PhyloFun )

# Test buildVectorSpaceModel
print("Testing buildVectorSpaceModel(...)")
f <- file( project.file.path( "test_annotations.tbl"), "r" )
annos <- unserialize(f)
close(f)
ua <- uniq.annotations(annos, type='InterPro')
vsm <- buildVectorSpaceModel(ua)
checkEquals( vsm,
  c( "IPR001957", "IPR003593", "IPR010921", "IPR013159", "IPR013317",
    "IPR018312", "IPR020591", "IPR024633" ))

# Read test domain weights 
domain.weights <- read.table(project.file.path( "domainWeights.tbl"))
