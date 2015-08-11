newConvenientScanonevar <- function(cross,
                                    mean.formula,
                                    var.formula,
                                    n.perms = 100,
                                    n.cores = 8) 
{
  
  validated.input <- ValidateAndFilterInput(cross, mean.formula, var.formula)
  genoprobs <- validated.input$genoprobs
  mapping.df <- validated.input$mapping.df
  chr.by.marker <- validated.input$chr.by.marker
  pos.by.marker <- validated.input$pos.by.marker
  marker.names <- validated.input$marker.names
  mean.formula <- validated.input$mean.formula
  var.formula <- validated.input$var.formula
  
  # do the actual scan
  scan <- ScanViaDGLM(mean.alt.formula = mean.formula,
                      var.alt.formula = var.formula,
                      genoprobs = genoprobs,
                      mapping.df = mapping.df,
                      chr.by.marker = chr.by.marker,
                      pos.by.marker = pos.by.marker,
                      marker.names = marker.names)
  
  # set up parallel permutations
  cl <- makeCluster(n.perms)
  registerDoParallel(cl)
  
  perms <- foreach(perm = 1:n.perms) %dorng% {
    
    perm.scan <- ScanViaDGLM(mean.alt.formula = mean.formula,
                             var.alt.formula = var.formula,
                             genoprobs = genoprobs[sample(nrow(genoprobs)),],
                             mapping.df = mapping.df,
                             chr.by.marker = chr.by.marker,
                             pos.by.marker = pos.by.marker,
                             marker.names = marker.names)
    
    perm.scan %>%
      select(full.lod, mean.lod, var.lod) %>%
      summarise_each(funs(max))
    
  }
  
  stopCluster(cl)
  
  # evaluate actual scan as quantiles of EVD from permutations scans
  perms <- bind_rows(perms)
  for (lod.type in c('full.lod', 'mean.lod', 'var.lod')) {

    scan[[paste0('emp.p.', lod.type)]] <- 1
    
    if (sum(is.na(scan[[lod.type]])) < 0.5*length(scan[[lod.type]])) {
      
      evd <- fgev(perms[[lod.type]])
      
      scan[[paste0('emp.p.', lod.type)]] <- pgev(q = scan[[lod.type]],
                                                 loc = fitted(evd)[1],
                                                 scale = fitted(evd)[2],
                                                 shape = fitted(evd)[3],
                                                 lower.tail = FALSE)
    }    
  }
  
  attr(scan, 'units') <- 'emp.ps'
  return(scan)
}