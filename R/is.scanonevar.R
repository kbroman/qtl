#'  @title is.scanonevar
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description Tests whether an object is a valid scanonevar object
#'  
#'  @param sov the object which is tested for being a scanonevar object
#'  
#'  @return Returns TRUE if 'sov' is a valid scanonevar object, with attribute 'why.not' an empty list
#'    Returns FALSE if 'sov' is not a valid scanonevar object with attribute 'why.not' a list of
#'    error messages for failed tests
#'    
#'  @seealso \code{\link{scanonevar}}, \code{\link{convert.scanonevar.to.empirical.ps}}

is.scanonevar <- function(sov) {
  
  ret <- TRUE
  attr(ret, 'why.not') <- list()
  
  # check classes
  if (!(all(class(sov) == c('scanonevar', 'tbl_df', 'tbl', 'data.frame')))) {
    ret <- FALSE
    attr(ret, 'why.not') <- c(attr(ret, 'why.not'), "classes must be c('scanonevar', 'tbl_df', 'tbl', 'data.frame')")
  }
  
  
  # check that required columns are present
  check.for.col.name <- function(col.name) {
    if (!(col.name %in% names(sov))) {
      ret <- FALSE
      attr(ret, 'why.not') <- c(attr(ret, 'why.not'), paste("no column named '", col.name,"'"))
    }
  }
  required.col.names <- c('chr', 'chrtype', 'pos', 'marker.name')
  for (col.name in required.col.names) {
    check.for.col.name(col.name = col.name)
  }
  
  
  # check that required attributes are present
  check.for.attr <- function(attr.name) {
    if (is.null(attr(sov, attr.name))) {
      ret <- FALSE
      attr(ret, 'why.not') <- c(attr(ret, 'why.not'), paste("no attr named '", attr.name,"'"))
    }
  }
  required.attr.names <- c('pheno', 'units', 'mean.alt.formula', 'var.alt.formula', 'null.fit')
  for (attr.name in required.attr.names) {
    check.for.attr(attr.name = attr.name)
  }
  
  # check that data values and ranges are valid
  if (attr(sov, 'units') == 'lods') {
    lod.columns <- grep(pattern = 'lod', names(scan), value = TRUE)
    
    for (lod.column in lod.columns) {
      
      if (any(sov[[lod.column]] < 0)) {
        ret <- FALSE
        attr(ret, 'why.not') <- c(attr(ret, 'why.not'), paste("Some LOD scores in", lod.column," are less than zero"))
      }
    }
  }
  
  if (attr(sov, 'units') == 'emp.ps') {
    lod.columns <- grep(pattern = 'emp.p', names(scan), value = TRUE)
    
    for (lod.column in lod.columns) {
      
      if (any(sov[[lod.column]] < 0 | sov[[lod.column]] > 1)) {
        ret <- FALSE
        attr(ret, 'why.not') <- c(attr(ret, 'why.not'), paste("Some empirical p values in", lod.column," are outside [0, 1]"))
      }
    }
  }
  
  return(ret)
}