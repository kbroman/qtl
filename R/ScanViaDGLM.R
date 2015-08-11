ScanViaDGLM <- function(mean.alt.formula,
                        var.alt.formula,
                        genoprobs,
                        mapping.df,
                        chr.by.marker,
                        pos.by.marker,
                        marker.names) { 
  
  # figure out which tests to do based on where the QTL term appears
  # also make null formulae programmatically
  # todo: allow user to specity null formulae optionally
  {
    test.mean.effect <- test.var.effect <- test.meanvar.effect <- FALSE
    mean.terms <- terms(mean.alt.formula)
    mean.qtl.terms <- grep(pattern = 'QTL', x = attr(mean.terms, 'term.labels'))
    if (any(mean.qtl.terms)) {
      mean.null.terms <- drop.terms(termobj = mean.terms,
                                    dropx = mean.qtl.terms)
      mean.null.formula <- reformulate(termlabels = attr(mean.null.terms, 'term.labels'),
                                       response = mean.alt.formula[[2]])
      test.mean.effect <- TRUE
    }
    
    var.terms <- terms(var.alt.formula)
    var.qtl.terms <- grep(pattern = 'QTL', x = attr(var.terms, 'term.labels'))
    if (any(var.qtl.terms)) {
      var.null.terms <- drop.terms(termobj = var.terms,
                                   dropx = var.qtl.terms)
      var.null.formula <- reformulate(termlabels = attr(var.null.terms, 'term.labels'))
      test.var.effect <- TRUE
    }
    
    if (!any(test.mean.effect, test.var.effect)) {
      stop('No genetic effects specified (use QTL.add and QTL.dom in mean.formula or var.formula).')
    }
    test.meanvar.effect <- test.mean.effect & test.var.effect
  }
  
  # calculate necessary LODs and their differences
  full.lod <- mean.lod <- var.lod <- rep(0, length(marker.names))
  log10lik.full <- log10lik.mean <- log10lik.var <- log10lik.null <- NA
  
  # calculate null for double-test
  if (test.meanvar.effect) {
    null.fit <- dglm(formula = mean.null.formula,
                     dformula = var.null.formula,
                     data = mapping.df)
    log10lik.null <- -0.5*null.fit$m2loglik / log(10)  
  }
  
  # loop through markers (including pseudomarkers)
  for (marker.idx in 1:length(marker.names)) {
    
    # select relevant markers from genoprob df
    marker.name <- marker.names[marker.idx]
    marker.genoprobs <- select(genoprobs, starts_with(marker.name))
    
    # for chr 1 - 19 in F2
    if (ncol(marker.genoprobs) == 3) {
      mapping.df$QTL.add <- AdditiveCoefficientFrom3Genoprobs(marker.genoprobs)
      mapping.df$QTL.dom <- DominantCoefficientFrom3Genoprobs(marker.genoprobs)
    }
    
    # for chr X in F2...maybe also for all of backcross?
    if (ncol(marker.genoprobs) == 2) {
      next
      mapping.df$QTL.add <- AdditiveCoefficientFrom2Genoprobs(marker.genoprobs)
      mapping.df$QTL.dom <- NA
    }
    
    full.fit <- dglm(formula = mean.alt.formula,
                     dformula = var.alt.formula,
                     data = mapping.df)
    log10lik.full <- -0.5*full.fit$m2loglik / log(10)
    
    if (test.var.effect) {
      mean.fit <- dglm(formula = mean.alt.formula,
                       dformula = var.null.formula,
                       data = mapping.df)
      log10lik.mean <- -0.5*mean.fit$m2loglik / log(10)
    }
    
    if (test.mean.effect) {
      var.fit <- dglm(formula = mean.null.formula,
                      dformula = var.alt.formula,
                      data = mapping.df)
      log10lik.var <- -0.5*var.fit$m2loglik / log(10)
    }
    
    
    full.lod[marker.idx] <- log10lik.full - log10lik.null
    mean.lod[marker.idx] <- log10lik.full - log10lik.var
    var.lod[marker.idx] <- log10lik.full - log10lik.mean
    
  }
 
  varscan <- data_frame(chr = factor(chr.by.marker,
                                     levels = mixedsort(unique(chr.by.marker))),
                        pos = pos.by.marker,
                        marker.name = marker.names,
                        full.lod,
                        mean.lod,
                        var.lod)  
  
  class(varscan) <- c('scanonevar', class(varscan))
  attr(varscan, 'units') <- 'lods'
  if (test.mean.effect) { attr(varscan, 'mean.null.formula') <- mean.null.formula }
  if (test.var.effect) { attr(varscan, 'var.null.formula') <- var.null.formula }
  attr(varscan, 'mean.alt.formula') <- mean.alt.formula
  attr(varscan, 'var.alt.formula') <- var.alt.formula
  
  return(varscan)
}