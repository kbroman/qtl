#'  @title Compute Additive Coefficient From Two Genotype Probabilities
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{get.additive.coef.from.2.genoprobs} should not typically be called by a user.
#'    It is used to reliably set up model coefficients from genotype probabilities.
#'  
#'  @param pa A vector of length n giving the probability of A-type for each individual, 
#'    or a matrix or data.frame of dimension n-by-2 giving probability of A-type and 
#'    the probability of B-type for each individual.
#'  @param pb A vector of length n giving the probability of B-type for each individual,
#'    or missing if \code{pa} is a data.frame or matrix.
#'    
#'  @return Additive coefficient vector.
#'  
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