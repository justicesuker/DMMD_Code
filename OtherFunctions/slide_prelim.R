# Generate S matrix in SLIDE model corresponding to known total and joint ranks.
getS <- function(r1, r2, r_c, ndata){
  S = matrix(rep(0, ndata * (r1 + r2 - r_c)),nrow = ndata)
  if (r_c > 0){
      S[ , 1:r_c] = 1
  }
  if (r1-r_c>0){
    S[1, (1 + r_c):r1] = 1
  }
  if (r2-r_c>0){
    S[2, (1+r1):(r2 + r1 - r_c)] = 1
  }
  return(S)
}

# getS(8,6,6,2)
# getS(8,6,0,2)
# getS(0,0,0,2)
# getS(6,5,5,2)
# getS(4,3,1,2)