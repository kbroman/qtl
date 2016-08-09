is.cross.w.genoprobs <- function(x) {

  # check that there is a 'prob' element in each chromosome
  if (!all(sapply(X = x[['geno']], FUN = function(x) 'prob' %in% names(x))))
    return(FALSE)

  # check that the dimension of each prob element matches the number of individuals
  is.valid.prob.array <- function(a, cross, chr) {
    array.dim <- dim(a)
    if (array.dim[1] != qtl::nind(cross))
      return(FALSE)
    if (array.dim[2] < ncol(chr[['data']]))
      return(FALSE)
    if (array.dim[3] != 3)
      return(FALSE)
    return(TRUE)
  }
  if (!all(sapply(X = x[['geno']], FUN = function(chr) is.valid.prob.array(a = chr[['prob']], cross = x, chr = chr))))
    return(FALSE)

  # all values in the prob array must be between 0 and 1
  sapply(X = x[['geno']], FUN = function(x) all(x[['prob']] > 0 & x[['prob']] < 1))

  return(is.cross(x))
}
