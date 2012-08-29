print(paste("Usage: Rscript runPhyloFun.R input=myQueryProteins.fasta",
"[jackhmmer.path=/path/to/hmmer3/binaries/jackhmmer]",
"[fast.tree.path=/path/to/FastTreeMP] [jackhmmer.args='-E 1e-10 -N 2']",
"[fast.tree.args='']"))

argsToList <- function(args) {
  sapply(args, function(a) {
      if (grepl("=",a,fixed=T)) {
        kv <- strsplit(a,'=')
        setNames(list(kv[[1]][[2]]),kv[[1]][[1]])
      }   
    }, USE.NAMES=F)
}

phyloFunArgs <- function(args.list, 
  defaults=list('jackhmmer.path'='jackhmmer',
    'fast.tree.path'='FastTreeMP',
    'jackhmmer.args'='-E 1e-10 -N 2',
    'fast.tree.args'='')) {
  sapply(names(defaults),function(n) {
      if(is.null(args.list[[n]]))
          setNames(list(defaults[[n]]),n)
      else
          setNames(list(args.list[[n]]),n)
    }, USE.NAMES=F)
}

input.args <- phyloFunArgs(argsToList(commandArgs(trailingOnly = TRUE)))
print(input.args)
