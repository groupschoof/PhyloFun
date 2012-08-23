library(gRain)

yn   <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub|asia, values=c(5,95,1,99),levels=yn)
s    <- cptable(~smoke, values=c(5,5), levels=yn)
l.s  <- cptable(~lung|smoke, values=c(1,9,1,99), levels=yn)
b.s  <- cptable(~bronc|smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either|lung:tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
x.e  <- cptable(~xray|either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp|bronc:either, values=c(9,1,7,3,8,2,1,9), levels=yn)

plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
in1 <- grain(plist)
querygrain(in1,nodes=c("lung","bronc"), type="marginal")
querygrain(in1,nodes=c("lung","bronc"), type="joint")
in12 <- setFinding(in1,nodes=c("asia","dysp"),states=c("yes","yes"))
querygrain(in12,nodes=c("lung","bronc"))
querygrain(in12,nodes=c("lung","bronc"), type="joint")
in13 <- setFinding(in1,nodes=c("either","tub"),states=c("no","yes"))
pFinding(in13)
querygrain(in13,nodes=c("lung","bronc"), type="joint")


