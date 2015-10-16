#'  @title Compute Additive Coefficient From Three Genotype Probabilities
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{get.add.coef.from.3.genoprobs} should not typically be called by a user.
#'    It is used to reliably set up model coefficients from genotype probabilities.
#'  
#'  @param paa A vector of length n giving the probability of AA-type for each individual, 
#'    or a matrix or data.frame of dimension n-by-3 giving probability of AA-type and 
#'    the probability of AB-type and the probability of BB-type for each individual.
#'  @param pab A vector of length n giving the probability of AB-type for each individual,
#'    or missing if \code{pa} is a data.frame or matrix.
#'  @param pbb A vector of length n giving the probability of BB-type for each individual,
#'    or missing if \code{pa} is a data.frame or matrix.
#'    
#'  @return Additive coefficient vector.
#'  
#'  
get.additive.coef.from.3.genoprobs <- function(paa, pab, pbb) {

  if (missing(pab) & missing(pbb) & ncol(paa) == 3 & is.list(paa)) {
    pab <- paa[[2]]
    pbb <- paa[[3]]
    paa <- paa[[1]]
  } else if (missing(pab) & missing(pbb) & ncol(paa) == 3 & !is.list(paa)) {
    pab <- paa[,2]
    pbb <- paa[,3]
    paa <- paa[,1]
  }  
  
  return(-paa + pbb)
}