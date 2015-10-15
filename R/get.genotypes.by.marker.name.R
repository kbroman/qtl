get.genotypes.by.marker.name <- function(cross, marker.name, use.genoprobs = TRUE, as.matrix = TRUE) {
  
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