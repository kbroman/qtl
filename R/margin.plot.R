margin.plot <- function(cross, 
                        focal.phenotype.name,
                        marginal.phen.names = NULL,
                        marginal.marker.names = NULL,
                        ...) {
  
  if (any(missing(cross), missing(focal.phenotype.name), !(focal.phenotype.name %in% names(cross$pheno)))) { stop('Must provide a cross and a focal phenotype in that cross.')}
  
  num.plots <- sum(length(marginal.phen.names), length(marginal.marker.names))
  if (num.plots == 0) { stop('Must provide a marginal phenotype or marker.')}
  
  par(mfrow = c(1, num.plots))
  
  focal.phen <- cross$pheno[[focal.phenotype.name]]
  
  for (marginal.phen.name in marginal.phen.names) {
    
    plot(x = cross$pheno[[marginal.phen.name]], y = focal.phen, xlab = marginal.phen.name, ...)
  }
  
  for (marginal.marker.name in marginal.marker.names) {
    
    chr.of.interest <- which(sapply(X = cross$geno, FUN = function(chr) { marginal.marker.name %in% colnames(chr$data)}))
    
    plot(x = jitter(cross$geno[[chr.of.interest]]$data[,marginal.marker.name]),
         y = focal.phen, 
         xlab = marginal.marker.name,
         ...)
  }
  
  
}