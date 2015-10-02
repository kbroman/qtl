simplex.plot.scanonevar <- function(cross, 
                                    name.of.marker.to.plot,
                                    varscan,
                                    marker.or.pseudomarker = 'marker',
                                    show.sample.sd = TRUE,
                                    horizontal.connectors = TRUE,
                                    main = NULL,
                                    col = 'gray',
                                    ...)
{
  
  
  # store current graphical parameters and customize them for this plot
  set.seed(27599)
  library(scatterplot3d)
  library(ggtern)
 
  # pull out phenotype values
  phen.name <- attr(x = varscan, which = 'pheno')
  phen <- cross$pheno[,phen.name]
  if (is.null(phen)) { stop('Phenotype of scan not found in cross')}
  
  # pull out genotypes at the marker to plot
  chr.of.genoprobs <- as.character(varscan$chr[varscan$marker.name == name.of.marker.to.plot])
  genoprobs <- cross$geno[[chr.of.genoprobs]]$prob[,strsplit(x = name.of.marker.to.plot, split = '_')[[1]][2],]
  
  tern.df <- data.frame(cbind(phen, genoprobs))
  
  simplex.y  <- function(x) {
    return( sqrt(0.75) *  x[3] / sum(x) )
  } 
  simplex.x  <- function(x) {
    return( (x[2] + 0.5 * x[3]) / sum(x) )
  }
  to.simplex <- function(v) {
    return(c(simplex.x(v), simplex.y(v)))
  }
  
  genoprobs.on.simplex <- t(apply(X = genoprobs, MARGIN = 1, FUN = to.simplex))
  
  scatterplot3d(x = genoprobs.on.simplex[,1],
                y = genoprobs.on.simplex[,2],
                z = phen)
  
  ggtern(data = tern.df, aes(x = AA, y = AB, z = BB)) + stat_density_tern(contour = FALSE)
  
  
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
  mtext(text = phen.name, side = 2, line = 2.2, cex = 1.2, las = 0)
  
  
  axis(side = 1, at = 1:3,
       labels = c('AA', 'AB', 'BB'),
       line = 0)
  
  if (show.sample.sd) {
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
}