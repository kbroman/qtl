margin.plot <- function(cross, 
                        focal.phenotype.name,
                        marginal.phen.names = NULL,
                        marginal.marker.names = NULL,
                        genotype.plotting.names = c('A', 'H', 'B'),
                        ...) {
  
  if (any(missing(cross), missing(focal.phenotype.name), !(focal.phenotype.name %in% names(cross$pheno)))) { 
    stop('Must provide a cross and a focal phenotype in that cross.')
  }
  
  num.plots <- sum(length(marginal.phen.names), length(marginal.marker.names))
  if (num.plots == 0) { stop('Must provide a marginal phenotype or marker.')}
  
  # store current graphical parameters and customize them for this plot
  start.pars <- par(no.readonly = TRUE)
  
  par(mfrow = c(1, num.plots))
  
  focal.phen <- cross$pheno[[focal.phenotype.name]]
  
  for (marginal.phen.name in marginal.phen.names) {
    
    marginal.phen <- cross$pheno[[marginal.phen.name]]
    if (is.factor(marginal.phen)) {
      lev.names <- levels(marginal.phen)
      plotting.phen <- as.numeric(marginal.phen)
    } else {
      plotting.phen <- marginal.phen
    }
    
    plot(x = jitter(plotting.phen), 
         y = focal.phen, 
         xlab = marginal.phen.name, 
         xaxt = 'n',
         main = paste(focal.phenotype.name, 'by', marginal.phen.name),
         ...)
    
    if (is.factor(marginal.phen)) {
      axis(side = 1, at = unique(marginal.phen), labels = lev.names, tick = TRUE)
    }
  
  }
  
  for (marginal.marker.name in marginal.marker.names) {
    
    chr.of.interest <- which(sapply(X = cross$geno, FUN = function(chr) { marginal.marker.name %in% colnames(chr$data)}))
    
    genotypes <- cross$geno[[chr.of.interest]]$data[,marginal.marker.name]
    plot(x = jitter(genotypes),
         y = focal.phen,
         xaxt = 'n',
         xlab = marginal.marker.name,
         main = paste(focal.phenotype.name, 'by', marginal.marker.name),
         ...)
    axis(side = 1, at = 1:3, labels = genotype.plotting.names, tick = TRUE)
    
    means <- aggregate(x = focal.phen, by = list(genotypes), FUN = mean)[,2]
    sds <- aggregate(x = focal.phen, by = list(genotypes), FUN = sd)[,2]
    
    x.start.wide <- 1:3 - 0.2
    x.end.wide <- 1:3 + 0.2
    x.start.narrow <- 1:3 - 0.1
    x.end.narrow <- 1:3 + 0.1
    
    # horizontal lines at means and means +/- 1SD
    segments(x0 = c(x.start.narrow, x.start.wide, x.start.wide),
             y0 = c(means, means - sds, means + sds),
             x1 = c(x.end.narrow, x.end.wide, x.end.wide),
             y1 = c(means, means - sds, means + sds),
             lwd = rep(c(4, 2, 2), each = 3))
    
    # vertical line down the middle of each genotype group
    segments(x0 = 1:3,
             y0 = means - sds,
             x1 = 1:3,
             y1 = means + sds,
             lwd = 2)
    
    # dotted horizontal lines connecting means and SD's
    segments(x0 = rep(1:2, 3),
             y0 = c(means[1:2], means[1:2] + sds[1:2], means[1:2] - sds[1:2]),
             x1 = rep(2:3, 3),
             y1 = c(means[2:3], means[2:3] + sds[2:3], means[2:3] - sds[2:3]),
             lty = 2)
  }
  
  # reset graphical parameteers to how they were on start
  par(start.pars)
  
  # return nothing
  invisible()
}