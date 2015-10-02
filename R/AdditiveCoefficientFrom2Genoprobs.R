AdditiveCoefficientFrom2Genoprobs <- function(pa, pb) {
  
  if (missing(pb) & ncol(pa) == 2 & is.list(pa)) {
    pb <- pa[[2]]
    pa <- pa[[1]]
  } else if (missing(pb) & ncol(pa) == 2 & !is.list(pa)) {
    pb <- pa[,2]
    pa <- pa[,1]
  }  
  
  return(-pa + pb)
}