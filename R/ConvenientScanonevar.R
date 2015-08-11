ConvenientScanoneVar <- function(cross,
                                 pheno.name,
                                 chrs = 1:length(cross$geno),
                                 mean.covar.names = NULL,
                                 var.covar.names = NULL,
                                 num.perms = 5,
                                 num.processes = 1,
                                 mean.add = TRUE,
                                 mean.dom = TRUE,
                                 var.add = TRUE,
                                 var.dom = TRUE,
                                 return.effects = TRUE,
                                 quiet = FALSE,
                                 seed = 27599, # ZIP code of UNC Chapel Hill
                                 mean.perm = TRUE,
                                 var.perm = TRUE,
                                 meanvar.perm = TRUE) {
  
  set.seed(seed)
  
  # check that pheno.name and covar.names exist in the cross
  # and that cross is of an acceptable type
  # return the phenotype and covariate data (not names)
  validated.input <- ValidateScanonevarInput(cross = cross,
                                             pheno.name = pheno.name,
                                             mean.add = mean.add,
                                             mean.dom = mean.dom,
                                             var.add = var.add,
                                             var.dom = var.dom,
                                             mean.covar.names = mean.covar.names,
                                             var.covar.names = var.covar.names)
  cross <- validated.input$cross
  pheno <- validated.input$pheno
  mean.covars <- validated.input$mean.covars
  var.covars <- validated.input$var.covars
  
  # set up design matrix (X) and formulae for use in dglm
  X.and.formulae <- AssembleDesignMatAndFormulae(pheno = pheno,
                                                 mean.add = mean.add,
                                                 mean.dom = mean.dom,
                                                 var.add = var.add,
                                                 var.dom = var.dom,
                                                 mean.covars = mean.covars,
                                                 var.covars = var.covars)
  X <- X.and.formulae$X
  mean.null.formula <- X.and.formulae$mean.null.formula
  var.null.formula <- X.and.formulae$var.null.formula
  mean.alt.formula <- X.and.formulae$mean.alt.formula
  var.alt.formula <- X.and.formulae$var.alt.formula
  
  # workhorse of scanonevar()
  # applies dglm to each locus and compares to various null fits to return LOD for each locus
  scan <- scanonevar_(cross = cross,
                      chrs = chrs,
                      mean.add = mean.add,
                      mean.dom = mean.dom,
                      var.add = var.add,
                      var.dom = var.dom,
                      X = X,
                      mean.null.formula = mean.null.formula,
                      var.null.formula = var.null.formula,
                      mean.alt.formula = mean.alt.formula,
                      var.alt.formula = var.alt.formula,
                      return.effects = return.effects,
                      quiet = quiet)
  
  # iterate through num.perms, calling the workhorse of scanonevar.perm() each time
  # do this in parallel if num.cores > 1, and in series otherwise
  if (num.processes == 1) {
    
    perms <- vector(mode = 'list', length = num.perms)
    for (perm.num in 1:num.perms) {
      
      if (!quiet) { print(paste('Starting permutation number', perm.num, '...')) }
      perm <- sample(nrow(pheno))
      
      meanvarperm <- scanonevar_(cross = cross, 
                                 chrs = chrs, 
                                 mean.add = mean.add,
                                 mean.dom = mean.dom,
                                 var.add = var.add,
                                 var.dom = var.dom,
                                 X = X,
                                 mean.null.formula = mean.null.formula,
                                 var.null.formula = var.null.formula,
                                 mean.alt.formula = mean.alt.formula,
                                 var.alt.formula = var.alt.formula,
                                 mean.perm = perm, 
                                 var.perm = perm,
                                 return.effects = FALSE,
                                 calc.lod.full = TRUE,
                                 calc.lod.mean = FALSE, 
                                 calc.lod.var = FALSE,
                                 quiet = TRUE)
      
      meanperm <- scanonevar_(cross = cross, 
                              chrs = chrs, 
                              mean.add = mean.add,
                              mean.dom = mean.dom,
                              var.add = var.add,
                              var.dom = var.dom,
                              X = X,
                              mean.null.formula = mean.null.formula,
                              var.null.formula = var.alt.formula,
                              mean.alt.formula = mean.alt.formula,
                              var.alt.formula = var.alt.formula,
                              mean.perm = perm, 
                              return.effects = FALSE,
                              calc.lod.full = FALSE,
                              calc.lod.mean = TRUE, 
                              calc.lod.var = FALSE,
                              quiet = TRUE)
      
      
      varperm <- scanonevar_(cross = cross, 
                             chrs = chrs, 
                             mean.add = mean.add,
                             mean.dom = mean.dom,
                             var.add = var.add,
                             var.dom = var.dom,
                             X = X,
                             mean.null.formula = mean.alt.formula,
                             var.null.formula = var.null.formula,
                             mean.alt.formula = mean.alt.formula,
                             var.alt.formula = var.alt.formula,
                             var.perm = perm,
                             return.effects = FALSE,
                             calc.lod.full = FALSE,
                             calc.lod.mean = FALSE, 
                             calc.lod.var = TRUE,
                             quiet = TRUE)
      
      n.loci <- nrow(meanvarperm)
      result <- cbind(perm.idx = rep(perm.num, n.loci),
                      as.data.frame(meanvarperm[,c('chrtype', 'chr', 'pos', 'lod.full')]),
                      as.data.frame(meanperm[, 'lod.mean', drop = FALSE]),
                      as.data.frame(varperm[, 'lod.var', drop = FALSE]))
      
      result <- tbl_df(result) %>%
        group_by(perm.idx, chrtype) %>%
        summarize(max.lod.full = max(lod.full),
                  max.lod.mean = max(lod.mean),
                  max.lod.var = max(lod.var))
      perms[[perm.num]] <- result
    }
  }
  
  if (num.processes > 1) {
    
    # set up parllel backend for permutations
    cl <- makeCluster(num.processes, outfile = '')
    registerDoParallel(cl)
    
    perms <- foreach(perm.num = 1:num.perms, .errorhandling = 'remove', .verbose = quiet) %dorng% {
      
      if (!quiet) { print(paste('Starting permutation number', perm.num, '...')) }
      perm <- sample(nrow(pheno))
      
      meanvarperm <- scanonevar_(cross = cross, 
                                 chrs = chrs, 
                                 mean.add = mean.add,
                                 mean.dom = mean.dom,
                                 var.add = var.add,
                                 var.dom = var.dom,
                                 X = X,
                                 mean.null.formula = mean.null.formula,
                                 var.null.formula = var.null.formula,
                                 mean.alt.formula = mean.alt.formula,
                                 var.alt.formula = var.alt.formula,
                                 mean.perm = perm, 
                                 var.perm = perm,
                                 return.effects = FALSE,
                                 calc.lod.full = TRUE,
                                 calc.lod.mean = FALSE, 
                                 calc.lod.var = FALSE)
      
      meanperm <- scanonevar_(cross = cross, 
                              chrs = chrs, 
                              mean.add = mean.add,
                              mean.dom = mean.dom,
                              var.add = var.add,
                              var.dom = var.dom,
                              X = X,
                              mean.null.formula = mean.null.formula,
                              var.null.formula = var.alt.formula,
                              mean.alt.formula = mean.alt.formula,
                              var.alt.formula = var.alt.formula,
                              mean.perm = perm, 
                              return.effects = FALSE,
                              calc.lod.full = FALSE,
                              calc.lod.mean = TRUE, 
                              calc.lod.var = FALSE)
      
      
      varperm <- scanonevar_(cross = cross, 
                             chrs = chrs, 
                             mean.add = mean.add,
                             mean.dom = mean.dom,
                             var.add = var.add,
                             var.dom = var.dom,
                             X = X,
                             mean.null.formula = mean.alt.formula,
                             var.null.formula = var.null.formula,
                             mean.alt.formula = mean.alt.formula,
                             var.alt.formula = var.alt.formula,
                             var.perm = perm,
                             return.effects = FALSE,
                             calc.lod.full = FALSE,
                             calc.lod.mean = FALSE, 
                             calc.lod.var = TRUE)
      
      n.loci <- nrow(meanvarperm)
      result <- cbind(perm.idx = rep(perm.num, n.loci),
                      as.data.frame(meanvarperm[,c('chrtype', 'chr', 'pos', 'lod.full')]),
                      as.data.frame(meanperm[, 'lod.mean', drop = FALSE]),
                      as.data.frame(varperm[, 'lod.var', drop = FALSE]))
      
      result <- tbl_df(result) %>%
        group_by(perm.idx, chrtype) %>%
        summarize(max.lod.full = max(lod.full),
                  max.lod.mean = max(lod.mean),
                  max.lod.var = max(lod.var))
      
      return(result)
    }
    
    # close down parallel backend
    stopCluster(cl)
  }
  
  perms <- rbind_all(perms)
  
  # convert to empirical p values and return that
  return(vqtl_ConvertLODsToEmpPs(scan, perms))
  
}
