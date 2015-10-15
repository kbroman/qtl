get.additive.coef.from.2.genoprobs <- function(pa, pb) {
  
  if (missing(pb) & ncol(pa) == 2 & is.list(pa)) {
    pb <- pa[[2]]
    pa <- pa[[1]]
  } else if (missing(pb) & ncol(pa) == 2 & !is.list(pa)) {
    pb <- pa[,2]
    pa <- pa[,1]
  }  
  
  return(-pa + pb)
}