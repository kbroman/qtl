validate.convert.scanonevar.to.emp.ps <- function(scan, null.scan.maxes) {
  
  if (!identical(unique(scan$chrtype), unique(null.scan.maxes$chrtype))) {
    stop('chrtypes of scan and null.scan.maxes dont match')
  }
  
  if ('lod.full' %in% names(scan) & !('max.lod.full' %in% names(null.scan.maxes))) {
    stop('no lod.full in null.scan.maxes, though its in the scan')
  }
  
  if ('lod.mean' %in% names(scan) & !('max.lod.mean' %in% names(null.scan.maxes))) {
    stop('no lod.mean in null.scan.maxes, though its in the scan')
  }
  
  if ('lod.var' %in% names(scan) & !('max.lod.var' %in% names(null.scan.maxes))) {
    stop('no lod.var in null.scan.maxes, though its in the scan')
  }
  
}