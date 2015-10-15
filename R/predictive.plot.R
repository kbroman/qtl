predictive.plot <- function(cross,
                             mean.formula, var.formula, 
                             marker.name, phen.name,
                             genotype.plotting.names = NULL) {
  
  phenotypes <- cross$pheno[[phen.name]]
  genoprobs <- GetGenoprobsByMarkerName(cross = cross, marker.name = marker.name)
  genotypes <- GetGenotypesByMarker(cross = cross, marker.name = marker.name, as.matrix = FALSE)
  response.phen <- mean.formula[[2]]
  
  if (!is.null(genotype.plotting.names)) {
    if (length(genotype.plotting.names) != length(unique(genotypes))) {
      stop('length of genotype.plotting.names must be equal to the number of genotypes at marker.name')
    }
    
    plotting.genotypes <- mapvalues(x = genotypes, from = unique(genotypes), to = genotype.plotting.names)
  }
#   if (!is.null(phenotype.plotting.names)) {
#     if (length(phenotype.plotting.names) != length(unique(phenotypes))) {
#       stop('length of phenotype.plotting.names must be equal to the number of phenotypes')
#     }
#     
#     phenotypes <- mapvalues(x = phenotypes, from = unique(phenotypes), to = phenotype.plotting.names)
#   }
  
  
  
  model.df <- data.frame(mean.QTL.add = AdditiveCoefficientFrom3Genoprobs(genoprobs),
                         mean.QTL.dom = DominantCoefficientFrom3Genoprobs(genoprobs),
                         var.QTL.add = AdditiveCoefficientFrom3Genoprobs(genoprobs),
                         var.QTL.dom = DominantCoefficientFrom3Genoprobs(genoprobs))
  
  model.df <- cbind(model.df, cross$pheno)
  
  dglm.fit <- dglm(formula = mean.formula,
                   dformula = var.formula,
                   data = model.df)
  
  
  mean.pred <- predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  
  var.pred <- predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  var.estim <- var.pred$fit/var.pred$residual.scale
  var.se <- var.pred$se.fit/var.pred$residual.scale
  
  prediction.tbl <- data_frame(genotype = genotypes,
                               plotting.genotype = plotting.genotypes,
                               phen = phenotypes,
                               indiv.mean.estim = mean.estim,
                               indiv.mean.lb = mean.estim - mean.se,
                               indiv.mean.ub = mean.estim + mean.se,
                               indiv.var.estim = exp(var.estim),
                               indiv.var.lb = exp(var.estim - var.se),
                               indiv.var.ub = exp(var.estim + var.se))
  
  plotting.tbl <- prediction.tbl %>% 
    group_by(phen, genotype, plotting.genotype) %>%
    summarise(group.mean.estim = mean(indiv.mean.estim),
              group.mean.lb = mean(indiv.mean.lb),
              group.mean.ub = mean(indiv.mean.ub),
              group.var.estim = mean(indiv.var.estim),
              group.var.lb = mean(indiv.var.lb),
              group.var.ub = mean(indiv.var.ub)) %>%
    arrange(phen, genotype)
  
  # set up colors -- special red/blue is phen is 'sex'
  unique.phens <- sort(unique(plotting.tbl$phen))
  num.unique.phens <- length(unique.phens)
  if (phen.name == 'sex') {
    phen.group.colors <- c('red', 'blue')
  } else {
    phen.group.colors <- brewer.pal(n = num.unique.phens, name = 'Set1')
  }
  plotting.tbl$phen.col <- mapvalues(x = plotting.tbl$phen, from = unique.phens, to = phen.group.colors)
 
  # blank plot of correct size 
  with(plotting.tbl,
       plot(x = 1, y = 1,
            type = 'n',
            xlab = NA,
            ylab = NA,
            xlim = range(c(group.mean.lb, group.mean.ub)), 
            ylim = range(c(group.var.lb, group.var.ub)),
            main = paste('Predictive of', response.phen, 'from', phen.name, 'and', marker.name)))
  mtext(side = 1, text = paste(response.phen, 'mean'), line = 2)
  mtext(side = 2, text = paste(response.phen, 'SD'), line = 2)
  
  # horizontal lines
  with(plotting.tbl,
       segments(x0 = group.mean.lb, 
                y0 = group.var.estim,
                x1 = group.mean.ub, 
                y1 = group.var.estim,
                col = phen.col, 
                lwd = 3))
  
  # vertical lines
  with(plotting.tbl,
       segments(x0 = group.mean.estim, 
                y0 = group.var.lb,
                x1 = group.mean.estim, 
                y1 = group.var.ub,
                col = phen.col,
                lwd = 3))
  
  # draw light lines connecting same-phenotype groups
  # e.g. connect male AA, male AB, and male BB together (same for female)
  for (phen.group.idx in 1:num.unique.phens) {
    phen <- unique.phens[phen.group.idx]
    col <- col2rgb(phen.group.colors[phen.group.idx])
    
    this.group.idxs <- plotting.tbl$phen == phen
    
    lines(x = plotting.tbl$group.mean.estim[this.group.idxs],
          y = plotting.tbl$group.var.estim[this.group.idxs],
          col = rgb(red = col[1], green = col[2], blue = col[3], maxColorValue = 255, alpha = 50),
          lwd = 20)
  }
  
  # white circles clear space to write genotype names
  with(plotting.tbl,
       points(x = group.mean.estim, y = group.var.estim,
         pch = 19, col = rgb(1, 1, 1), cex = 5))
  
  # write genotype names
  with(plotting.tbl,
       text(x = group.mean.estim, y = group.var.estim,
            labels = plotting.genotype, cex = 1.5,
            col = phen.col))
  
  
  
#   p <- ggplot(plotting.tbl, mapping = aes(x = group.mean.estim, y = group.var.estim))
#   p <- p + geom_segment(mapping = aes(color = rep(factor(phen), 2),
#                                       x = c(group.mean.estim, group.mean.lb), xend = c(group.mean.estim, group.mean.ub),
#                                       y = c(group.var.lb, group.var.estim), yend = c(group.var.ub, group.var.estim)),
#                         size = 2)
#   p <- p + geom_point(color = 'white', size = 15)
#   p <- p + geom_text(aes(label = genotype), size = 10)
#   p <- p + theme_bw()
#   print(p)
}