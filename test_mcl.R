sheep <- function(x) {
  out <- vector(mode='numeric')
  out <- append(out, as.numeric(x))
  for(i in 1:10){
    print(out[i])
    out <- append(out, (out[i]^2))
  }
  # return
  out
}

goat <- function(x, sums) {
  out <- vector(mode='numeric')
  for(i in 1:length(x)) {
    out <- append(out, (x[i] / sums[i]))
  }
  # return
  out
}



