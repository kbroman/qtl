#' @title is.cross
#' @name is.cross
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @param x The object being tested for whether or not is is a cross.
#'
#' @return TRUE if x is a cross object, FALSE otherwise.
#' @export
#'
#' @examples
#' is.cross(3)
#' is.cross(qtl::sim.cross(map = qtl::sim.map()))
#'
is.cross <- function(x) {

  # class of x
  if(!('cross' %in% class(x)))
    return(FALSE)

  # pheno, an element of x
  if (!('pheno' %in% names(x)))
    return(FALSE)

  if (!('data.frame' %in% class(x[['pheno']])))
    return(FALSE)


  # geno, an element of x
  if (!('geno' %in% names(x)))
    return(FALSE)

  if (!('list' %in% class(x[['geno']])))
    return(FALSE)


  # chromosomes, elements of geno
  if (!all(sapply(X = x[['geno']], FUN = class) %in% c('A', 'X')))
    return(FALSE)

  if (!all(sapply(X = x[['geno']], FUN = function(x) ('data' %in% names(x)))))
      return(FALSE)


  # consistent number of individuals across pheno and all genos
  pheno.n <- nrow(x[['pheno']])
  geno.ns <- sapply(X = x[['geno']], FUN = function(x) nrow(x[['data']]))
  if (any(pheno.n != geno.ns))
    return(FALSE)

  # length of map matches number of genotypes in each chromosome
  if (!all(sapply(X = x[['geno']], FUN = function(x) ncol(x[['data']]) == length(x[['map']]))))
    return(FALSE)


  return(TRUE)
}



#' @title is.f2.cross
#' @name is.f2.cross
#' @author Robert W. Corty \email{rcorty@@gmail.com}
#'
#' @inheritParams is.cross
#'
#' @return TRUE if x is a cross object of type F2, FALSE otherwise
#' @export
#'
#' @examples
#' is.cross(3)
#' is.cross(qtl::sim.cross(map = qtl::sim.map()))
#'
is.f2.cross <- function(x) {

  if (!('f2' %in% class(x)))
    return(FALSE)

  return(is.cross(x))
}


