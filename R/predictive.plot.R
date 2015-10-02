predictive.plot <- function(varscan,
                            cross,
                            phen.name,
                            marker.name,
                            plot.by.genoprobs = FALSE, 
                            ...) {
  
  if (any(length(phen.name) != 1, length(marker.name) != 1)) {
    stop('Must input exactly one phenotype name and one marker name')
  }
  
  means <- vars <- mean.ses <- var.ses <- rep(0, nind(cross))
  
  # get design matrix of phenotype and genotype
  phen <- as.numeric(cross$pheno[[phen.name]]) - 1
  chr.of.interest <- which(sapply(X = cross$geno, FUN = function(chr) { marker.name %in% colnames(chr$data)}))
  genoprobs <- cross$geno[[chr.of.interest]]$prob[,marker.name,]
  vec.max <- function(v) {
    ret <- rep(0, length(v))
    ret[which.max(v)] <- 1
    return(ret)
  }
  genotypes <- t(apply(X = genoprobs, MARGIN = 1, FUN = vec.max))
  
  if (!plot.by.genoprobs) {
    add.coef <- AdditiveCoefficientFrom3Genoprobs(genotypes)
    dom.coef <- DominantCoefficientFrom3Genoprobs(genotypes)
  }
  
  if (plot.by.genoprobs) {
    add.coef <- AdditiveCoefficientFrom3Genoprobs(genoprobs)
    dom.coef <- DominantCoefficientFrom3Genoprobs(genoprobs)
  }
  
  # get the right row out of the varscan
  myrow <- varscan[which(attr(varscan$pos, 'names') == marker.name),]
  if (nrow(myrow) != 1) { stop('Marker must match exactly one row in varscan.')}
  
  means <- myrow[['mean.intercept']] + phen * myrow[[paste0('mean.', phen.name)]] + add.coef * myrow[['mean.mean.QTL.add']] + dom.coef * myrow[['mean.mean.QTL.dom']]
  vars <- myrow[['var.intercept']] + phen * myrow[[paste0('var.', phen.name)]] + add.coef * myrow[['var.var.QTL.add']] + dom.coef * myrow[['var.var.QTL.dom']]
  
  mean.ses <- phen * myrow[[paste0('se.mean.', phen.name)]] + add.coef * myrow[['se.mean.mean.QTL.add']] + dom.coef * myrow[['se.mean.mean.QTL.dom']]
  var.ses <- phen * myrow[[paste0('se.var.', phen.name)]] + add.coef * myrow[['se.var.var.QTL.add']] + dom.coef * myrow[['se.var.var.QTL.dom']]
  
  plot(x = means, y = vars, 
       xlim = range(c(means + mean.ses, means - mean.ses)), 
       ylim = range(c(vars + var.ses, vars - var.ses)))
  segments(x0 = means, y0 = vars, x1 = means + mean.ses, y1 = vars)
  segments(x0 = means, y0 = vars, x1 = means - mean.ses, y1 = vars)
  segments(x0 = means, y0 = vars, x1 = means           , y1 = vars + var.ses)
  segments(x0 = means, y0 = vars, x1 = means           , y1 = vars - var.ses)
  
  genotypes <- factor(apply(X = genotypes, MARGIN = 1, FUN = which.max))
  levels(genotypes) <- c('AA', 'AB', 'BB')
  
  to.plot <- !duplicated(means)
  text(x = means[to.plot], y = vars[to.plot], labels = cross$pheno[[phen.name]][to.plot], pos = 2, cex = 3)
  text(x = means[to.plot], y = vars[to.plot], labels = genotypes[to.plot], pos = 1, cex = 3)
}