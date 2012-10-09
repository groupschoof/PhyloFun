library(RUnit)
library(tools)
library(RCurl)

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
src.project.file( 'src','domainArchitectureSimilarity.R' )

# Load domain.weights.table:
dwt <- read.table(
  project.file.path( "data", "uniprotkb_based_domain_weights.tbl" )
  )[ , 'DW', drop=F ]

# Load test annotations:
f <- file( project.file.path( "test", "test_annotations_2.tbl" ), "r" )
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
dwd <- read.table( project.file.path( "test", "domain_weights_database.tbl" ),
  header=T )
# Initialize expected result:
f <- file( project.file.path( "test", "testDomainArchitectureSpaceVectors.txt"),
  "r" )
exp.dasv <- unserialize( f )
close(f)

dasv <- generateDomainArchitectureSpaceVectors( vsm, am, dwd )
checkEquals( exp.dasv, dasv )

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

# Test pairwiseSequenceDistance
print("Testing pairwiseSequenceDistance(...)")
p1 <- "LALDTKQIWFTTLGTLQNQILRYDYDTWLKTTALVSVANDLAVIGAPNVTTKQVIEDRFMSVLRRALGEVLGYQVNVRVIISSATPAPSEPVAVTPSEPSPTTEVAEPSFASFNQAAPMLNQLPLGDPNRSSVLNPRYTFSSFIVGTSNRLAHAACMAVAEHPAQAYNPLFLYGGVGLGKTHLLQAIGNYALDRNPEVNVLYVSSEKFTNDLINAIRRQQTEEFRIRYRNIDILLIDDIQFIAGKEGTQEEFFHTFNTLHGAGKQIVLSSDRPPKAILTLEERLRSRFEWGLIVDVQNPDLETRTAILRAKGETLQVPVSSEVIDFLAQRIQSNIRELEGCLNRVIAYANLNRTPVTVEVASAALADLLDTSRRKRVTADDIFREVSQHYGIDQRAIRGRGRSRNVVLPRQVVMYLLREETDASLVEIGELLGGRDHTTVMHGYNKITDDLTSDARLRNDITSLRQRLYGENAR"
p2 <- "MQSIEDIWQETLQIVKKNMSKPSYDTWMKSTTAHSLEGNTFIISAPNNFVRDWLEKSYTQFIANILQEITGRLFDVRFIDGEQEENFEYTVIKPNPALDEDGIEIGKHMLNPRYVFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAVGHYVQQHKDNAKVMYLSSEKFTNEFISSIRDNKTEEFRTKYRNVDVLLIDDIQFLAGKEGTQEEFFHTFNTLYDEQKQIIISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKADGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLVNKDITAGLAAEALKDIIPSSKSQVITISGIQETVGEYFHVRLEDFKAKKRTKSIAFPRQIAMYLSRELTDASLPKIGDEFGGRDHTTVIHAHEKISQLLKTDQVLKNDLAEIEKNLRKSQNMF"
p.dist <- round( pairwiseSequenceDistance( p1, p2 ), 6 )
checkTrue( ! is.null(p.dist) )
checkEquals( class(p.dist), 'numeric' )
checkEquals( p.dist, 0.813185 )
