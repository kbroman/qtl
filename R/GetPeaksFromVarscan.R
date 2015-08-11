GetPeaksFromVarscan <- function(vs, thresh) {
  
  if (units(vs) == 'lods') {
    
    if (missing(thresh)) { thresh <- 3 }
    
    peaks <- vs %>%
      mutate(full.peak = and(lod.full > lag(lod.full, default = 0),
                             lod.full > lead(lod.full, default = 0))) %>%
      mutate(mean.peak = and(lod.mean > lag(lod.mean, default = 0),
                             lod.mean > lead(lod.mean, default = 0))) %>%
      mutate(var.peak = and(lod.var > lag(lod.var, default = 0),
                            lod.var > lead(lod.var, default = 0)))
  }
  
  if (units(vs) == 'emp.ps') {
    
    if (missing(thresh)) { thresh <- 0.05 }
    
    peaks <- vs %>%
      mutate(full.peak = and(emp.p.full < lag(emp.p.full, default = 1),
                             emp.p.full < lead(emp.p.full, default = 1))) %>%
      mutate(mean.peak = and(emp.p.mean < lag(emp.p.mean, default = 1),
                             emp.p.mean < lead(emp.p.mean, default = 1))) %>%
      mutate(var.peak = and(emp.p.var < lag(emp.p.var, default = 1),
                            emp.p.var < lead(emp.p.var, default = 1)))
  }
  
  return(peaks)
}