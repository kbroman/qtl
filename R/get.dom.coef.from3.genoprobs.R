get.dom.coef.from.3.genoprobs <- function(paa, pab, pbb) {
  
  if (missing(pab) & missing(pbb) & ncol(paa) == 3 & is.list(paa)) {
    pab <- paa[[2]]
    pbb <- paa[[3]]
    paa <- paa[[1]]
  } else if (missing(pab) & missing(pbb) & ncol(paa) == 3 & !is.list(paa)) {
    pab <- paa[,2]
    pbb <- paa[,3]
    paa <- paa[,1]
  }  
  
  return(pab)
}