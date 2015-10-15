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