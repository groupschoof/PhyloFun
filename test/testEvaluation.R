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
src.project.file <- function(...) {
  source(project.file.path(...))
}
src.project.file( 'src', 'evaluation.R' )

# Test precision
print("Testing precision(...)")
true.gos <- c( 'A', 'B', 'C' )
res.precision <- precision( true.gos, true.gos )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( true.gos, 'D' ), true.gos )
exp.precision <- 0.75
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( 'A', 'B' ), true.gos )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( 'D' ), true.gos )
exp.precision <- 0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( 'A' ), c() )
exp.precision <- 0
checkEquals( res.precision, exp.precision ) 

res.precision <- precision( c( ), true.gos )
exp.precision <- 1.0
checkEquals( res.precision, exp.precision ) 

# Test recall
print("Testing recall(...)")
true.gos <- c( 'A', 'B', 'C' )
res.recall <- recall( true.gos, true.gos )
exp.recall <- 1.0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( true.gos, 'D' ), true.gos )
exp.recall <- 1.0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( 'A', 'B' ), true.gos )
exp.recall <- 2 / 3
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( 'D' ), true.gos )
exp.recall <- 0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( 'A' ), c() )
exp.recall <- 1.0
checkEquals( res.recall, exp.recall ) 

res.recall <- recall( c( ), true.gos )
exp.recall <- 0
checkEquals( res.recall, exp.recall ) 

# Test fScore
print("Testing fScore(...)")
res.fScore <- fScore( true.gos, true.gos )
exp.fScore <- 1.0
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c( 'A' ), true.gos )
exp.fScore <- 0.5
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c(), c() )
exp.fScore <- 1.0
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c(), true.gos )
exp.fScore <- 0
checkEquals( res.fScore, exp.fScore ) 

res.fScore <- fScore( c( 'A' ), c() )
exp.fScore <- 0
checkEquals( res.fScore, exp.fScore ) 
