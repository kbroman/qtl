ConvenientScanoneVar <- function(cross,
                                 pheno.name,
                                 chrs = 1:length(cross$geno),
                                 mean.covar.names = NULL,
                                 var.covar.names = NULL,
                                 dom = TRUE,
                                 return.effects = TRUE,
                                 quiet = TRUE) {
  
  # check that pheno.name and covar.names exist in the cross
  # and that cross is of an acceptable type
  # return the phenotype and covariate data (not names)
  validated.input <- ValidateScanonevarInput(cross = cross,
                                             pheno.name = pheno.name,
                                             mean.covar.names = mean.covar.names,
                                             var.covar.names = var.covar.names)
  cross <- validated.input$cross
  pheno <- validated.input$pheno
  mean.covars <- validated.input$mean.covars
  var.covars <- validated.input$var.covars
  
  # set up design matrix (X) and formulae for use in dglm
  X.and.formulae <- assemble.design.mat.and.formulae(pheno = pheno,
                                                     dom = dom,
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
                      dom = dom,
                      X = X,
                      mean.null.formula = mean.null.formula,
                      var.null.formula = var.null.formula,
                      mean.alt.formula = mean.alt.formula,
                      var.alt.formula = var.alt.formula,
                      return.effects = return.effects)
}
