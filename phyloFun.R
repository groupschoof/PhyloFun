# PhyloFun v0.1
# Authors: Asis Hallab, Heiko Schoof
library(gRain)
library(Matrix)
# library(gdata)

# Test-Tree
# (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
#  
# GO-Annotations for Tips:
# A -> GO:1
# B -> GO:2
# C -> GO:3
# D -> QUERY-PROTEIN

go_terms <- c("GO:1","GO:2","GO:3")
no_evidence <- c(1,1,1)
spec <- 0.03

# Use unity matrix with entries 1.0
base.transition.probs <- matrix(1,
  dimnames=list(go_terms,go_terms),
  nrow=3,ncol=3)

# Compute conditional probability tables:
# NOTE: 'values' is parsed row by row
# yn  <- c("yes","no")
# a   <- cptable(~asia, values=c(1,99),levels=yn)
# t.a <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)
# t.a's conditional probability table depending on its parent a 
# then is:
#         a=yes a=no
# t.a=yes     5   95
# t.a=no      1   99
# Having been in Asia the probability of having conceived tuberculosis is 5
# times higher than normal.
 
node.f <- cptable(~F, values=no_evidence, levels=go_terms)
node.a <- cptable(~A | F, values=expm(0.1*spec*base.transition.probs), levels=go_terms)
node.b <- cptable(~B | F, values=expm(0.2*spec*base.transition.probs), levels=go_terms)
node.e <- cptable(~E | F, values=expm(0.5*spec*base.transition.probs), levels=go_terms)
node.c <- cptable(~C | E, values=expm(0.3*spec*base.transition.probs), levels=go_terms)
node.d <- cptable(~D | E, values=expm(0.4*spec*base.transition.probs), levels=go_terms)

# Compile Bayesian Network
plist <- compileCPT(list(node.f,node.a,node.b,node.c,node.d,node.e))
b.n   <- grain(plist)

# Enter evidence
# 1.
# phylo.fun <- setFinding(b.n,
#   nodes=c("A","B","C"),
#   states=c("GO:1","GO:1","GO:1")
#   )
#
# querygrain(phylo.fun,nodes=c("D"))

# 2.
# phylo.fun <- setFinding(b.n,
#   nodes=c("A","B","C"),
#   states=c("GO:1","GO:2","GO:3")
#   )
#
# querygrain(phylo.fun,nodes=c("D"))

# Returns the same as setFinding and subsequent querygrain
# predict's arguments are.
# 1. Network without any evidence
# 2. Nodes whose probability distributions are to be predicted
# 3. Nodes having the following evidence
# 4. The evidence in a data.frame
res <- predict(b.n, 
  c("F","E","D"),       
  c("A","B","C"),   
  data.frame(list(A="GO:1",B="GO:2",C="GO:3")),
  type="distribution"
  )
