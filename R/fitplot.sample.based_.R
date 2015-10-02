fitplot.sample.based_ <- function(cross, 
                                   name.of.marker.to.plot,
                                   phenotype.name,
                                   main,
                                   col,
                                   ...) {
  
  # pull out phenotype values
  phen <- cross$pheno[,phenotype.name]
  if (is.null(phen)) { stop('Phenotype of scan not found in cross')}
  
  # pull out genotypes at the marker to plot
  orig.name <- strsplit(x = name.of.marker.to.plot, split = '_', fixed = TRUE)[[1]][-1]
  genotypes <- GetGenoprobsByMarkerName(cross = cross,
                                        marker.name = orig.name)
  genotypes <- apply(X = genotypes, MARGIN = 1, FUN = which.max)
  
 
  # find all the peaks on this chromosome
  # name.of.chr.with.marker <- mles$chr
  # 	chr.peaks <- GetPeaksFromVarscan(vs = dplyr::filter(varscan, chr == name.of.chr.with.marker))
  # 
  # 	# distance of desired marker to each peak
  # 	full.peaks <- dplyr::filter(chr.peaks, full.peak == TRUE)
  # 	dists <- abs(mles$pos - (full.peaks)$pos)
  # 	best.full.locus <- slice(full.peaks, which.min(dists))
  # 
  # 	mean.peaks <- dplyr::filter(chr.peaks, mean.peak == TRUE)
  # 	dists <- abs(mles$pos - (mean.peaks)$pos)
  # 	best.mean.locus <- slice(mean.peaks, which.min(dists))
  # 
  # 	var.peaks <- dplyr::filter(chr.peaks, var.peak == TRUE)
  # 	dists <- abs(mles$pos - (var.peaks)$pos)
  # 	best.var.locus <- slice(var.peaks, which.min(dists))
  
  
  # plot phenotypes, broken out by genotype at this marker
  plot(x = genotypes + runif(n = nind(cross), min = -0.2, max = 0.2),
       y = phen,
       xlab = NA,
       ylab = NA,
       xaxt = 'n',
       xlim = c(0.5, 3.5),
       ylim = range(phen)*c(0.95, 1.05),
       pch = 20,
       col = col, 
       main = paste0(name.of.marker.to.plot, main),
       ...)
  mtext(text = 'Genotype', side = 1, line = 2.2, cex = 1.2)
  mtext(text = phenotype.name, side = 2, line = 2.2, cex = 1.2, las = 0)
  
  # 	if (attr(varscan, which = 'units') == 'lods') {
  # 		title(paste('Full Model: LOD here',
  # 								round(mles$lod.full, 2),
  # 								'LOD at closest full peak:',
  # 								round(best.full.locus$lod.full, 2),
  # 								'\nMean Model: LOD here',
  # 								round(mles$lod.mean, 2),
  # 								'LOD at closest mean peak:',
  # 								round(best.mean.locus$lod.mean, 2),
  # 								'\nVar Model: LOD here',
  # 								round(mles$lod.var, 2),
  # 								'LOD at closest var peak:',
  # 								round(best.var.locus$lod.var, 2)))
  # 	}
  # 	if (units(varscan) == 'emp.ps') {
  # 		title(paste('Full Model: -log10p here',
  # 								round(-log10(mles$emp.p.full), 2),
  # 								'-log10p at closest full peak:',
  # 								round(-log10(best.full.locus$emp.p.full), 2),
  # 								'\nMean Model: -log10p here',
  # 								round(-log10(mles$emp.p.mean), 2),
  # 								'-log10p at closest mean peak:',
  # 								round(-log10(best.mean.locus$emp.p.mean), 2),
  # 								'\nVar Model: -log10p here',
  # 								round(-log10(mles$emp.p.var), 2),
  # 								'-log10p at closest var peak:',
  # 								round(-log10(best.var.locus$emp.p.var), 2)))
  # 	}
  
  # autosome-specific stuff
  # if (mles$chrtype == 'A') {
  
  axis(side = 1, at = 1:3,
       labels = c('AA', 'AB', 'BB'),
       line = 0)
  
  means <- aggregate(x = phen, by = list(genotypes), FUN = mean)[,2]
  sds <- aggregate(x = phen, by = list(genotypes), FUN = sd)[,2]
  
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