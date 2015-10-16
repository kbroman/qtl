#'  @title Check the Compatibility of the Scanonevar to be Converted with the Permutations to be Used in the Conversion
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{validate.convert.scanonevar.to.emp.ps} should not typically be called by a user.
#'  This function is used by \code{convert.scanonevar.to.emp.ps}
#'    
#'  @param scan the scanonevar to be converted
#'  @param null.scan.maxes the maximum LODs observed in permutation (null) scans to be used in the conversion.
#'  
#'  @return Returns TRUE if the two arguments are compatible and FALSE otherwise.
#'    
#'  @seealso  \link{\code{convert.scanonevar.to.emp.ps}}, \code{\link{scanonevar}}, \code{\link{scanonevar.perm}}
#'  
#'  
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