newScanonevar.perm <- function(cross,
                               mean.formula,
                               var.formula,
                               n.perms,
                               chrs = unique(names(cross$geno))) 
{
  
  validated.input <- ValidateAndFilterInput(cross = cross, 
                                            mean.formula = mean.formula, 
                                            var.formula = var.formula,
                                            chrs = chrs)
  genoprobs <- validated.input$genoprobs
  mapping.df <- validated.input$mapping.df
  chr.by.marker <- validated.input$chr.by.marker
  pos.by.marker <- validated.input$pos.by.marker
  marker.names <- validated.input$marker.names
  mean.formula <- validated.input$mean.formula
  var.formula <- validated.input$var.formula
  
  all.perms <- NULL
  for (perm.idx in 1:n.perms) {
    perm.scan <- ScanViaDGLM(mean.alt.formula = mean.formula,
                             var.alt.formula = var.formula,
                             genoprobs = genoprobs,
                             mapping.df = mapping.df,
                             chr.by.marker = chr.by.marker,
                             pos.by.marker = pos.by.marker,
                             marker.names = marker.names,
                             perm = sample(nrow(genoprobs)))
    
    this.perm <- perm.scan %>%
      select(full.lod, mean.lod, var.lod) %>%
      summarise_each(funs(max(., na.rm = TRUE)))
    
    if (is.null(all.perms)) {
      all.perms <- this.perm
    } else {
      all.perms <- rbind(all.perms, this.perm)
    }
  
  }

  return(all.perms)
  
}