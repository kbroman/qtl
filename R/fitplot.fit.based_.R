
fitplot.fit.based_ <- function(cross, 
                               name.of.marker.to.plot,
                               varscan,
                               main,
                               col,
                               ...)	
{
  
  # pull out phenotype values
  phen.name <- attr(x = varscan, which = 'pheno')
  phen <- cross$pheno[,phen.name]
  if (is.null(phen)) { stop('Phenotype of scan not found in cross')}
  
  # pull out genotypes at the marker to plot
  orig.name <- strsplit(x = name.of.marker.to.plot, split = '_', fixed = TRUE)[[1]][-1]
  genoprobs <- GetGenoprobsByMarkerName(cross = cross,
                                        marker.name = orig.name)
  
  genotypes <- apply(X = genoprobs, MARGIN = 1, FUN = which.max)
  idealized.genoprobs <- t(apply(X = genoprobs, MARGIN = 1, FUN = function(v) { 
    ret <- c(0, 0, 0)
    ret[which.max(v)] <- 1
    return(ret)
  }))
  add.design <- AdditiveCoefficientFrom3Genoprobs(data.frame(idealized.genoprobs))
  dom.design <- DominantCoefficientFrom3Genoprobs(data.frame(idealized.genoprobs))
  
  # check that varscan included chr where marker is
  mles <- dplyr::filter(varscan, marker.name == name.of.marker.to.plot)
  if (nrow(mles) == 0) { stop('Desired marker is not in scan.')}
  
  # count number of covariates
  mean.covars <- as.character(attr(varscan, 'mean.null.formula')[[3]])
  var.covars <- as.character(attr(varscan, 'var.null.formula')[[2]])
  unique.covars <- unique(c(mean.covars, var.covars))
  
  # set up plotting areas
  par(mfrow = c(1, 1 + length(unique.covars)))
  
  # plot phenotypes broken up by genotype
  plot(x = genotypes + runif(n = nind(cross), min = -0.2, max = 0.2),
       y = phen,
       xlab = NA,
       ylab = NA,
       xaxt = 'n',
       xlim = c(0.5, 3.5),
       pch = 20,
       col = col, 
       main = paste0(name.of.marker.to.plot, main),
       ...)
  mtext(text = 'Genotype', side = 1, line = 2.2, cex = 1.2)
  mtext(text = phen.name, side = 2, line = 2.2, cex = 1.2, las = 0)
  
  axis(side = 1, at = 1:3,
       labels = c('AA', 'AB', 'BB'),
       line = 0)
  
  means <- mles$mean.intercept + mles$mean.sex*(as.numeric(cross$pheno$sex) - 1) + mles$mean.mean.QTL.add*add.design + mles$mean.mean.QTL.dom*dom.design
  vars <- mles$var.intercept + mles$var.sex*(as.numeric(cross$pheno$sex) - 1) + mles$var.var.QTL.add*add.design + mles$var.var.QTL.dom*dom.design

  group.means <- aggregate(x = means, by = list(genotypes), FUN = mean)[,2]
  group.vars <- aggregate(x = vars, by = list(genotypes), FUN = mean)[,2]
  

  x.start.wide <- 1:3 - 0.2
  x.end.wide <- 1:3 + 0.2
  x.start.narrow <- 1:3 - 0.1
  x.end.narrow <- 1:3 + 0.1
  
  # horizontal lines at means and means +/- 1SD
  segments(x0 = c(x.start.narrow, x.start.wide, x.start.wide),
           y0 = c(group.means, group.means - group.vars, group.means + group.vars),
           x1 = c(x.end.narrow, x.end.wide, x.end.wide),
           y1 = c(group.means, group.means - group.vars, group.means + group.vars),
           lwd = rep(c(4, 2, 2), each = 3))
  
  # vertical line down the middle of each genotype group
  segments(x0 = 1:3,
           y0 = group.means - group.vars,
           x1 = 1:3,
           y1 = group.means + group.vars,
           lwd = 2)
  
  # dotted horizontal lines connecting means and SD's
  segments(x0 = rep(1:2, 3),
           y0 = c(group.means[1:2], group.means[1:2] + group.vars[1:2], group.means[1:2] - group.vars[1:2]),
           x1 = rep(2:3, 3),
           y1 = c(group.means[2:3], group.means[2:3] + group.vars[2:3], group.means[2:3] - group.vars[2:3]),
           lty = 2)
  
  
  
  # plot phenotype broken up by other covariates
  for (covar in unique.covars) {
    
    plot(x = as.numeric(cross$pheno[[covar]]) + runif(n = nind(cross), min = -0.3, max = 0.3),
         y = phen,
         xlab = NA,
         ylab = NA,
         xaxt = 'n',
         pch = 20,
         main = covar,
         ...)
    mtext(text = covar, side = 1, line = 2.2, cex = 1.2)
    mtext(text = phen.name, side = 2, line = 2.2, cex = 1.2, las = 0)
    
    if (covar == 'sex') {
      axis(side = 1, at = 1:2,
           labels = c('female', 'male'),
           line = 0)
    }
    
    group.means <- aggregate(x = means, by = list(cross$pheno[[covar]]), FUN = mean)[,2]
    group.vars <- aggregate(x = vars, by = list(cross$pheno[[covar]]), FUN = mean)[,2]
    
    x.start.wide <- 1:length(group.means) - 0.2
    x.end.wide <- 1:length(group.means) + 0.2
    x.start.narrow <- 1:length(group.means) - 0.1
    x.end.narrow <- 1:length(group.means) + 0.1
    
    # horizontal lines at means and means +/- 1SD
    segments(x0 = c(x.start.narrow, x.start.wide, x.start.wide),
             y0 = c(group.means, group.means - group.vars, group.means + group.vars),
             x1 = c(x.end.narrow, x.end.wide, x.end.wide),
             y1 = c(group.means, group.means - group.vars, group.means + group.vars),
             lwd = rep(c(4, 2, 2), each = length(group.means)))
    
    # vertical line down the middle of each genotype group
    segments(x0 = 1:length(group.means),
             y0 = group.means - group.vars,
             x1 = 1:length(group.means),
             y1 = group.means + group.vars,
             lwd = 2)
    
    # dotted horizontal lines connecting means and SD's
    segments(x0 = rep(1:(length(group.means) - 1), 3),
             y0 = c(group.means[1:(length(group.means) - 1)], 
                    group.means[1:(length(group.means) - 1)] + group.vars[1:(length(group.means) - 1)], 
                    group.means[1:(length(group.means) - 1)] - group.vars[1:(length(group.means) - 1)]),
             x1 = rep(length(group.means), 3),
             y1 = c(group.means[length(group.means)],
                    group.means[length(group.means)] + group.vars[length(group.means)], 
                    group.means[length(group.means)] - group.vars[length(group.means)]),
             lty = 2)
  }
  
  
  
  
  
  
  
  
#   
#   means <- mles$mean.intercept + mles$mean.sex*(as.numeric(cross$pheno$sex) - 1) + mles$mean.mean.QTL.add*add.design + mles$mean.mean.QTL.dom*dom.design
#   vars <- mles$var.intercept + mles$var.sex*(as.numeric(cross$pheno$sex) - 1) + mles$var.var.QTL.add*add.design + mles$var.var.QTL.dom*dom.design
#   
  # 		# compute predictions for each group based on
  # 		# parameters fitted in scanonevar
  # 		aa.mean <- mles$mean.baseline 
  # 		ab.mean <- mles$mean.baseline
  # 		bb.mean <- mles$mean.baseline
  # 
  # 		if ('mean.add' %in% names(mles)) {
  # 		  aa.mean <- aa.mean - mles$mean.add
  # 		  bb.mean <- bb.mean + mles$mean.add
  # 		}
  # 		
  # 		if ('mean.dom' %in% names(mles)) {
  # 		  ab.mean <- ab.mean + mles$mean.dom
  # 		}
  # 		
  # 		aa.var <- mles$var.baseline 
  # 		ab.var <- mles$var.baseline
  # 		bb.var <- mles$var.baseline
  
  
  # 		if ('var.add' %in% names(mles)) {
# 		  aa.var <- aa.var - mles$var.add
# 		  bb.var <- bb.var + mles$var.add
# 		}
# 		
# 		if ('var.dom' %in% names(mles)) {
# 		  ab.var <- ab.var + mles$var.dom
# 		}
# }

# 	# X chromosome-specific stuff
# 	if (mles$chrtype == 'X') {
# 
# 		axis(side = 1, at = 1:3,
# 				 labels = c('A- or AA', 'AB', 'B- or BB'),
# 				 line = 0)
# 
# 		# compute predictions for each group based on
# 		# parameters fitted in scanonevar (additive only here)
# 	  aa.mean <- mles$mean.baseline 
# 	  ab.mean <- mles$mean.baseline
# 	  bb.mean <- mles$mean.baseline
# 	  
# 	  if ('mean.add' %in% names(mles)) {
# 	    ab.mean <- ab.mean + 0.5*mles$mean.add.effect
# 	    bb.mean <- bb.mean + mles$mean.add.effect
# 	  }
# 	  
# 	  aa.var <- mles$var.baseline 
# 	  ab.var <- mles$var.baseline
# 	  bb.var <- mles$var.baseline
# 	  
# 	  if ('var.add' %in% names(mles)) {
# 	    ab.var <- ab.var + 0.5*mles$var.add.effect
# 	    bb.var <- bb.var + mles$var.add.effect
# 	  }
# 	}


# 	aa.sd <- sqrt(exp(aa.var))
# 	ab.sd <- sqrt(exp(ab.var))
# 	bb.sd <- sqrt(exp(bb.var))
# 
# 	means <- c(aa.mean, ab.mean, bb.mean)
# 	sds <- c(aa.sd, ab.sd, bb.sd)
# 
# 	# vertical lines
# 	segments(x0 = (1:3),
# 					 y0 = means - sds,
# 					 x1 = (1:3),
# 					 y1 = means + sds,
# 					 lwd = 3,
# 					 col = 'black')
# 
# 	# fitted mean lines
# 	segments(x0 = (1:3) - 0.05,
# 					 x1 = (1:3) + 0.05,
# 					 y0 = means,
# 					 y1 = means,
# 					 lwd = 3,
# 					 col = 'black')
# 
# 	# 1 SD away from mean lines
# 	segments(x0 = rep((1:3) - 0.05, each = 2),
# 					 x1 = rep((1:3) + 0.05, each = 2),
# 					 y0 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
# 					 y1 = rep(means, each = 2) + c(aa.sd, -aa.sd, ab.sd, -ab.sd, bb.sd, -bb.sd),
# 					 lwd = 3,
# 					 col = 'black')
# 
# 	# horizontalish connectors
# 	segments(x0 = rep(1:2, 3),
# 					 x1 = rep(2:3, 3),
# 					 y0 = c(means[-3] + sds[-3], means[-3], means[-3] - sds[-3]),
# 					 y1 = c(means[-1] + sds[-1], means[-1], means[-1] - sds[-1]),
# 					 col = c(rep('red', 2), rep('blue', 2), rep('red', 2)))
  
}
