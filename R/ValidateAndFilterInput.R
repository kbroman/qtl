ValidateAndFilterInput <- function(cross, mean.formula, var.formula)
{
  
  # calc genoprobs if needed
  if (!('prob' %in% names(cross$geno[[1]]))) {
    cross <- calc.genoprob(cross = cross, step = 2)
  }
  
  # turn genoprobs into one big tbl_df
  for (this.chr.name in names(cross$geno))
  {
    # pull out genoprobs and store in tbl_df
    this.chr.probs <- cross$geno[[this.chr.name]]$prob

    # add chr name to pseudomarker names so they don't collide
    chr.specific.names <- paste0('chr', 
                                 this.chr.name,
                                 '_',
                                 colnames(this.chr.probs))
    
    colnames(this.chr.probs) <- chr.specific.names
    
    this.chr.tbl <- tbl_df(data.frame(this.chr.probs))

    # aggregate
    if (!exists('all.chr.genoprobs')) {
      all.chr.genoprobs <- this.chr.tbl
      chr.by.marker <- rep(this.chr.name, dim(this.chr.probs)[2])
      pos.by.marker <- attr(this.chr.probs, 'map')
      marker.names <- colnames(this.chr.probs)
    } else { 
      all.chr.genoprobs <- bind_cols(all.chr.genoprobs, this.chr.tbl)
      chr.by.marker <- c(chr.by.marker, rep(this.chr.name, dim(this.chr.probs)[2]))
      pos.by.marker <- c(pos.by.marker, attr(this.chr.probs, 'map'))
      marker.names <- c(marker.names, colnames(this.chr.probs))
    }
  }
  
  # get the names of all the variables needed (markers and phenotypes)
  mean.var.names <- all.vars(mean.formula, unique = TRUE)
  var.var.names <- all.vars(var.formula, unique = TRUE)

  # split var names into those that refer to phenotypes and those that refer to markers
  mean.pheno.vars <- mean.marker.vars <- NULL
  for (var in mean.var.names) {
    if (var == 'QTL.add' | var == 'QTL.dom') {
      next
    } else if (var %in% names(cross$pheno)) {
      mean.pheno.vars <- c(mean.pheno.vars, var)
    } else if (var %in% names(pos.by.marker)) {
      mean.marker.vars <- c(mean.marker.vars, var)
    } else  {
      stop(paste("Variable", var, "not found in phenotypes or genotypes."))
    }
  }
  
  var.pheno.vars <- var.marker.vars <- NULL
  for (var in var.var.names) {
    if (var == 'QTL.add' | var == 'QTL.dom') {
      next
    } else if (var %in% names(cross$pheno)) {
      var.pheno.vars <- c(var.pheno.vars, var)
    } else if (var %in% names(pos.by.marker)) {
      var.marker.vars <- c(var.marker.vars, var)
    } else {
      stop(paste("Variable", var, "not found in phenotypes or genotypes."))
    }
  }
  
  mapping.df <- data.frame(QTL.add = rep(NA, nrow(all.chr.genoprobs)),
                           QTL.dom = rep(NA, nrow(all.chr.genoprobs)))
  
  pheno.vars <- c(mean.pheno.vars, var.pheno.vars)
  if (!is.null(pheno.vars))  {
    pheno.vars <- cross$pheno %>%
      select(matches(paste(pheno.vars,
                           collapse = '|')))
    mapping.df <- bind_cols(mapping.df, pheno.vars)
  }
  
  marker.vars <- unique(c(mean.marker.vars, var.marker.vars))
  if (!is.null(marker.vars))  {
    
    # treat each marker var one at a time
    for (marker.var in marker.vars) {
      
      marker.df <- all.chr.genoprobs %>%
        select(-matches('loc')) %>%
        select(matches(marker.var))
      
      marker.df <- marker.df[,-1]
      
      if (ncol(marker.df) == 1) {
        names(marker.df) <- marker.var
        mapping.df <- bind_cols(mapping.df, marker.df)
      } else {
        if (marker.var %in% mean.var.names) {
          mean.terms <- attr(terms(mean.formula), 'term.labels')
          mean.formula <- reformulate(termlabels = gsub(pattern = marker.var, 
                                                        replacement = paste0('(',
                                                                             paste(names(marker.df),
                                                                                   collapse = ' + '),
                                                                             ')'),
                                                        x = mean.terms),
                                      response = mean.formula[[2]])
        }
        if (marker.var %in% var.var.names) {
          var.terms <- attr(terms(var.formula), 'term.labels')
          var.formula <- reformulate(termlabels = gsub(pattern = marker.var, 
                                                        replacement = paste0('(',
                                                                             paste(names(marker.df),
                                                                                   collapse = ' + '),
                                                                             ')'),
                                                        x = var.terms))
        }
      }
      
      mapping.df <- cbind(mapping.df, marker.df)
    }
  }
  
  return(list(genoprobs = all.chr.genoprobs,
              chr.by.marker = chr.by.marker,
              pos.by.marker = pos.by.marker,
              marker.names = marker.names,
              mapping.df = mapping.df,
              mean.formula = mean.formula,
              var.formula = var.formula))
}