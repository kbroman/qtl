#'  @title Get Genotypes From a Cross Object By Marker Name
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{get.genotypes.by.marker.name} is an accessor function that 
#'    returns the most likely genotype of each individual in the cross at the given marker.
#'  
#'  @param cross The cross from which the genetic information will be extracted.
#'  @param marker.name The name of the marker where we want to know each individuals most likely genotype.
#'  @param use.genoprobs Defaults to TRUE.  Should we look at the genoprobs to figure out the 
#'    most likely genotype?  This ensures that there will be no NA.  But in some cases there may be 
#'    significant uncertainty, so this may oversimplify the true situation.
#'  @param as.matrix Defaults to FALSE.  Should the resulting genotypes be returns as a vector of numeric
#'    values (default) or a matrix?
#'    
#'  @return Most likely genotype (or NA) for all individuals in the cross at the specified locus.
#'  
get.genotypes.by.marker.name <- function(cross, marker.name, use.genoprobs = TRUE, as.matrix = FALSE) {
  
  chr.of.interest <- which(sapply(X = cross$geno, FUN = function(chr) { marker.name %in% colnames(chr$data)}))

  if (use.genoprobs) {
    
    vec.max <- function(v) {
      ret <- rep(0, length(v))
      ret[which.max(v)] <- 1
      return(ret)
    }
    
    genoprobs <- cross$geno[[chr.of.interest]]$prob[,marker.name,]

    if (as.matrix) {
      genoprob.mat <- t(apply(X = genoprobs, MARGIN = 1, FUN = vec.max))
      colnames(genoprob.mat) <- c('AA', 'AB', 'BB')
      return(genoprob.mat)
    }
    if (!as.matrix) {
      return(apply(X = genoprobs, MARGIN = 1, FUN = which.max))
    }
  }
  
  if (!use.genoprobs) {
    
    if (as.matrix) {
      stop('Not sure it makes sense to use genotypes as matrix...what should the design matrix be for an individual with a missing genotype?')
    }
    
    if (!as.matrix) {
      return(cross$geno[[chr.of.interest]][['data']][,marker.name])
    }
  }
  
}