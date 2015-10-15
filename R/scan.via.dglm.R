scan.via.dglm <- function(mean.alt.formula,
                        var.alt.formula,
                        genoprobs,
                        mapping.df,
                        chr.by.marker,
                        pos.by.marker,
                        marker.names,
                        return.effects = FALSE,
                        return.effect.ses = FALSE,
                        return.effect.ps = FALSE,
                        cor.threshold = 0.8,
                        perm = 1:nrow(genoprobs))
{ 
  
#### figure out which tests to do based on where the QTL term appears ####
# also make null formulae programmatically
# todo: allow user to specity null formulae optionally
  test.mean.effect <- test.var.effect <- test.meanvar.effect <- FALSE
  mean.terms <- terms(mean.alt.formula)
  mean.qtl.terms <- grep(pattern = 'mean.QTL', x = attr(mean.terms, 'term.labels'))
  if (any(mean.qtl.terms)) {
    mean.null.terms <- drop.terms(termobj = mean.terms,
                                  dropx = mean.qtl.terms)
    mean.null.formula <- reformulate(termlabels = attr(mean.null.terms, 'term.labels'),
                                     response = mean.alt.formula[[2]])
    test.mean.effect <- TRUE
  }
  
  var.terms <- terms(var.alt.formula)
  var.qtl.terms <- grep(pattern = 'var.QTL', x = attr(var.terms, 'term.labels'))
  if (any(var.qtl.terms)) {
    var.null.terms <- drop.terms(termobj = var.terms,
                                 dropx = var.qtl.terms)
    var.null.formula <- reformulate(termlabels = attr(var.null.terms, 'term.labels'))
    test.var.effect <- TRUE
  }
  
  if (!any(test.mean.effect, test.var.effect)) {
    stop('No genetic effects specified (use (mean|var).QTL.(add|dom) in mean.formula or var.formula).')
  }
  test.meanvar.effect <- test.mean.effect & test.var.effect
  
#### set up containers for the scan ####
  full.lod <- mean.lod <- var.lod <- rep(NA, length(marker.names))
  log10lik.bothalt <- log10lik.bothnull <- log10lik.meanalt <- log10lik.meannull <- log10lik.varalt <- log10lik.varnull <- NA
  
  if (any(return.effects, return.effect.ses, return.effect.ps)) { 
    num.mean.effects <- length(all.vars(mean.alt.formula)) # minus 1 for response, plus one for intercept term
    num.var.effects <- length(all.vars(var.alt.formula)) + 1 # add one for the intercept
    num.markers <- length(marker.names)
    mean.indices <- 1:(num.mean.effects)
    var.indices <- (num.mean.effects + 1):(num.mean.effects + num.var.effects)
  }
  if (return.effects) {
    fitted.effects <- array(data = NA, dim = c(num.markers, num.mean.effects + num.var.effects),              # add one for null marker
                            dimnames = list(marker.names,
                                            c('mean.intercept', paste0('mean.', all.vars(mean.alt.formula)[-1]),
                                              'var.intercept', paste0('var.', all.vars(var.alt.formula)))))
  }
  if (return.effect.ses) {
    effect.ses <- array(data = NA, dim = c(num.markers, num.mean.effects + num.var.effects),              # add one for null marker
                        dimnames = list(marker.names,
                                        c('se.mean.intercept', paste0('se.mean.', all.vars(mean.alt.formula)[-1]),
                                          'se.var.intercept', paste0('se.var.', all.vars(var.alt.formula)))))
  }
  if (return.effect.ps) {
    effect.ps <- array(data = NA, dim = c(num.markers, num.mean.effects + num.var.effects),              # add one for null marker
                       dimnames = list(marker.names,
                                       c('p.mean.intercept', paste0('p.mean.', all.vars(mean.alt.formula)[-1]),
                                         'p.var.intercept', paste0('p.var.', all.vars(var.alt.formula)))))
  }
  
#### calculate null for double-test ####
  if (test.meanvar.effect) {
    both.null.fit <- dglm(formula = mean.null.formula,
                     dformula = var.null.formula,
                     data = mapping.df)
    log10lik.bothnull <- -0.5*both.null.fit$m2loglik / log(10)  
  }
  
#### loop through markers (including pseudomarkers) ####
  for (marker.idx in 1:length(marker.names)) {
    
    # select relevant markers from genoprob df
    marker.name <- marker.names[marker.idx]
    marker.genoprobs <- select(genoprobs, starts_with(paste0(marker.name, '_')))
    
    # need to set QTL.dom to something that won't throw a "sd = 0" error from cor()
    # just temporarily, it gets set back to 0 in a few lines
    mapping.df$mean.QTL.add <- mapping.df$mean.QTL.dom <- NA
    mapping.df$var.QTL.add  <- mapping.df$var.QTL.dom  <- NA
    
    # skip this locus if it is highly correlated with a marker covariate
    # todo: could skip only the test that has this marker as a covariate
    this.marker.vec <- marker.genoprobs[[2]]
    cors <- sapply(X = mapping.df, FUN = function(v) ifelse(test = is.numeric(v), cor(v, this.marker.vec, use = 'complete.obs'), FALSE))
    if (any(cors > cor.threshold)) { next }
    
    # BOTH TESTING
    if (ncol(marker.genoprobs) == 3) {
      mapping.df$mean.QTL.add <- mapping.df$var.QTL.add <- get.additive.coef.from.3.genoprobs(marker.genoprobs)[perm]
      mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- get.dom.coef.from.3.genoprobs(marker.genoprobs)[perm]
    }
    if (ncol(marker.genoprobs) == 2) {
      mapping.df$mean.QTL.add <- mapping.df$var.QTL.add <- get.additive.coef.from.2.genoprobs(marker.genoprobs)[perm]
      mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- 0
    }
    
    both.alt.fit <- dglm(formula = mean.alt.formula,
                     dformula = var.alt.formula,
                     data = mapping.df)
    log10lik.bothalt <- -0.5*both.alt.fit$m2loglik / log(10)
    
    if (any(return.effects, return.effect.ses, return.effect.ps)) { 
      if (ncol(marker.genoprobs == 3)) {
        
        mean.indices <- 1:(num.mean.effects)
        var.indices <- (num.mean.effects + 1):(num.mean.effects + num.var.effects)
      }
      if (ncol(marker.genoprobs) == 2) {
        
        mean.idx.to.rm <- grep(pattern = 'dom', x = colnames(effect.ps[,mean.indices]))
        mean.indices <- mean.indices[-mean.idx.to.rm] 
        
        var.idx.to.rm <- grep(pattern = 'dom', x = colnames(effect.ps[,var.indices]))
        var.indices <- var.indices[-var.idx.to.rm]
      } 
    }
    if (return.effects) {
      fitted.effects[marker.idx, mean.indices] <- summary(both.alt.fit)$coef[,'Estimate']
      fitted.effects[marker.idx, var.indices] <- summary(both.alt.fit$dispersion.fit)$coef[,'Estimate']
    }
    if (return.effect.ses) {
        effect.ses[marker.idx , mean.indices] <- summary(both.alt.fit)$coef[,'Std. Error']
        effect.ses[marker.idx, var.indices] <- summary(both.alt.fit$dispersion.fit)$coef[,'Std. Error']
    }
    if (return.effect.ps) {
        effect.ps[marker.idx, mean.indices] <- summary(both.alt.fit)$coef[,'Pr(>|t|)']
        effect.ps[marker.idx, var.indices] <- summary(both.alt.fit$dispersion.fit)$coef[,'Pr(>|t|)']
    }
    
    # store results
    full.lod[marker.idx] <- log10lik.bothalt - log10lik.bothnull
    
    # MEAN TESTING
    if (test.mean.effect) {
      
      # if we are not doing a permutation (doing a real scan) the mean alt is the same as the both alt
      if (identical(perm, 1:nrow(genoprobs))) {
        log10lik.meanalt <- log10lik.bothalt
        
      # if we are doing a permutation, we calculate a different mean alt
      } else {
        
        if (ncol(marker.genoprobs) == 3) {
          mapping.df$mean.QTL.add <- AdditiveCoefficientFrom3Genoprobs(marker.genoprobs)[perm]
          mapping.df$var.QTL.add <- AdditiveCoefficientFrom3Genoprobs(marker.genoprobs)
          mapping.df$mean.QTL.dom <- DominantCoefficientFrom3Genoprobs(marker.genoprobs)[perm]
          mapping.df$var.QTL.dom <- DominantCoefficientFrom3Genoprobs(marker.genoprobs)
        }
        if (ncol(marker.genoprobs) == 2) {
          mapping.df$mean.QTL.add <- AdditiveCoefficientFrom2Genoprobs(marker.genoprobs)[perm]
          mapping.df$var.QTL.add <- AdditiveCoefficientFrom2Genoprobs(marker.genoprobs)
          mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- 0
        }
        
        mean.alt.fit <- dglm(formula = mean.alt.formula,
                             dformula = var.alt.formula,
                             data = mapping.df)
        log10lik.meanalt <-  -0.5*mean.alt.fit$m2loglik / log(10)
      }
      
      mean.null.fit <- dglm(formula = mean.null.formula,
                            dformula = var.alt.formula,
                            data = mapping.df)
      log10lik.meannull <- -0.5*mean.null.fit$m2loglik / log(10)
      
      # store results
      mean.lod[marker.idx] <- log10lik.meanalt - log10lik.meannull
    }
    
    # VAR TESTING
    if (test.var.effect) {
      
      # if we are not doing a permutation (doing a real scan) the mean alt is the same as the both alt
      if (identical(perm, 1:nrow(genoprobs))) {
        log10lik.varalt <- log10lik.bothalt
        
        # if we are doing a permutation, we calculate a different mean alt
      } else {
        
        if (ncol(marker.genoprobs) == 3) {
          mapping.df$mean.QTL.add <- AdditiveCoefficientFrom3Genoprobs(marker.genoprobs)
          mapping.df$var.QTL.add <- AdditiveCoefficientFrom3Genoprobs(marker.genoprobs)[perm]
          mapping.df$mean.QTL.dom <- DominantCoefficientFrom3Genoprobs(marker.genoprobs)
          mapping.df$var.QTL.dom <- DominantCoefficientFrom3Genoprobs(marker.genoprobs)[perm]
        }
        if (ncol(marker.genoprobs) == 2) {
          mapping.df$mean.QTL.add <- AdditiveCoefficientFrom2Genoprobs(marker.genoprobs)
          mapping.df$var.QTL.add <- AdditiveCoefficientFrom2Genoprobs(marker.genoprobs)[perm]
          mapping.df$mean.QTL.dom <- mapping.df$var.QTL.dom <- 0
        }
        
        var.alt.fit <- dglm(formula = mean.alt.formula,
                             dformula = var.alt.formula,
                             data = mapping.df)
        log10lik.varalt <-  -0.5*var.alt.fit$m2loglik / log(10)
      }
      
      var.null.fit <- dglm(formula = mean.alt.formula,
                            dformula = var.null.formula,
                            data = mapping.df)
      
      log10lik.varnull <- -0.5*var.null.fit$m2loglik / log(10)
      
      # store results
      var.lod[marker.idx] <- log10lik.varalt - log10lik.varnull
    }
    
  }
 
#### compile data into return format ####
  varscan <- data_frame(chr = factor(chr.by.marker,
                                     levels = mixedsort(unique(chr.by.marker))),
                          pos = pos.by.marker,
                          marker.name = marker.names,
                          full.lod,
                          mean.lod,
                          var.lod) 
  if (return.effects) {
    varscan <- bind_cols(varscan, data.frame(fitted.effects))
  }
  if (return.effect.ses) {
    varscan <- bind_cols(varscan, data.frame(effect.ses))
  }
  if (return.effect.ps) {
    varscan <- bind_cols(varscan, data.frame(effect.ps))
  }
  
  class(varscan) <- c('scanonevar', class(varscan))
  attr(varscan, 'pheno') <- as.character(mean.null.formula[[2]])
  attr(varscan, 'units') <- 'lods'
  if (test.mean.effect) { attr(varscan, 'mean.null.formula') <- mean.null.formula }
  if (test.var.effect) { attr(varscan, 'var.null.formula') <- var.null.formula }
  attr(varscan, 'mean.alt.formula') <- mean.alt.formula
  attr(varscan, 'var.alt.formula') <- var.alt.formula
  attr(varscan, 'null.fit') <- both.null.fit
  
  return(varscan)
}