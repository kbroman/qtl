# function under development
# idea is to find local maxima


get.peaks.from.scanonevar <- function(vs, thresh) {
  
  if (units(vs) == 'lods') {
    
    if (missing(thresh)) { thresh <- 3 }
    
    peaks <- vs %>%
      mutate(full.peak = and(full.lod >= lag(full.lod, default = 0),
                             full.lod >= lead(full.lod, default = 0))) %>%
      mutate(mean.peak = and(mean.lod >= lag(mean.lod, default = 0),
                             mean.lod >= lead(mean.lod, default = 0))) %>%
      mutate(var.peak = and(var.lod >= lag(var.lod, default = 0),
                            var.lod >= lead(var.lod, default = 0)))
  }
  
  if (units(vs) == 'emp.ps') {
    
    if (missing(thresh)) { thresh <- 0.05 }
    
    peaks <- vs %>%
      mutate(full.peak = and(emp.p.full.lod <= lag(emp.p.full.lod, default = 1),
                             emp.p.full.lod <= lead(emp.p.full.lod, default = 1))) %>%
      mutate(mean.peak = and(emp.p.mean.lod <= lag(emp.p.mean.lod, default = 1),
                             emp.p.mean.lod <= lead(emp.p.mean.lod, default = 1))) %>%
      mutate(var.peak = and(emp.p.var.lod <= lag(emp.p.var.lod, default = 1),
                            emp.p.var.lod <= lead(emp.p.var.lod, default = 1)))
  }
  
  return(peaks)
}