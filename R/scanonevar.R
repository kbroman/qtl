#'  @title Conduct a Scanonevar.
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{scanonevar} conducts a scanonevar, a genome scan that takes into account 
#'    genetic and non-genetic effects on residual (environmental) variation in the trait of interest.
#'  
#'  @param cross The cross on which the scanonevar will be conducted.
#'  @param mean.formula The formula that describes the response, and the covariates and genetic effects that influence it.
#'    The left hand side of the ~ must be a single phenotype that is in the cross.
#'    The right hand side must use only phenotypes thata are in the cross, markers that are in the cross,
#'    and the special terms: mean.QTL.add (additive effect on the mean) and mean.QTL.dom (dominance deviation from additive on the mean).
#'  @param var.formula The formula that describes the covariates and the genetic effects that influence residual (environmental) variation.
#'    There should be nothing on the left of the ~ (Inferred to be residual variation).
#'    The right hand side must use only phenotypes thata are in the cross, markers that are in the cross,
#'    and the special terms: var.QTL.add (additive effect on the variance) and var.QTL.dom (dominance deviation from additive on the variance).
#'  @param return.effects Logical indicating whether the estimated effects should be returned.
#'  @param return.effect.ses Logical indicating whether the standard errors of the estimated effects should be returned.
#'  @param return.effect.ps Logical indicating whether the p-value of the estimated effects should be returned.
#'  @param chrs The subset of chromosomes to scan (defaults to all chromosomes).
#'  @param exclusion.window Numeric between 0 and 1 indicating how tightly a locus must be correlated with a covariate to be skipped.
#'    e.g. if cor.threshold is 0.8 (it's default) any locus with \code{cor(locus, covariate) > 0.8} will be skipped.
#'    
#'  @return A scanonevar object.
#'  
#'  
#'    
scanonevar <- function(cross,
                       mean.formula,
                       var.formula,
                       return.effects = FALSE,
                       return.effect.ses = FALSE,
                       return.effect.ps = FALSE,
                       chrs = unique(names(cross$geno)),
                       exclusion.window = 0.8)
{
  
  validated.input <- validate.input.scanonevar(cross = cross, 
                                               mean.formula = mean.formula,
                                               var.formula = var.formula, 
                                               chrs = chrs)
  
  scan <- scan.via.dglm(mean.alt.formula = validated.input$mean.formula,
                        var.alt.formula = validated.input$var.formula,
                        genoprobs = validated.input$genoprobs,
                        mapping.df = validated.input$mapping.df,
                        chr.by.marker = validated.input$chr.by.marker,
                        pos.by.marker = validated.input$pos.by.marker,
                        marker.names = validated.input$marker.names,
                        return.effects = return.effects,
                        return.effect.ses = return.effect.ses,
                        return.effect.ps = return.effect.ps,
                        cor.threshold = exclusion.window)
  
  return(scan)
}