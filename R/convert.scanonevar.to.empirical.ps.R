#'  @title Convert Scanonevar from LODs to Empirical p-values
#'  
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'  
#'  @description \code{convert.scanonevar.to.empirical.ps} takes a scanonevar with LODs as units
#'    and maxes from permutation scans, estimates an extreme value distribution for the maxes,
#'    and returns the probability of observing the LOD scores in those EVDs.
#'  
#'  @param scanonevar the \code{scanonevar} in LODs to be converted to emprirical p values
#'  @param perm.scan.maxes the tbl_df object returned by scanonevar.perm, the maximum LOD score
#'    observed on a per-scan, per-chromosome-type basis in permutation scans.
#'  
#'  @return Returns a scanonevar object in terms of p-values, with \code{attr(x, 'units') = 'emp.ps'}.
#'    
#'  @seealso  \code{\link{scanonevar}}, \code{\link{scanonevar.perm}}
#'  
convert.scanonevar.to.empirical.ps <- function(scanonevar, perm.scan.maxes) {

  validate.convert.scanonevar.to.emp.ps(scanonevar, perm.scan.maxes)

	lod.columns <- grep(pattern = 'lod', names(scanonevar), value = TRUE)
	chr.types <- unique(scanonevar$chrtype)
	scan.as.emp.ps <- scanonevar

  # do each chromosome type and lod column separately
	# currently: empirical p values are put in place of the LOD scores
	# todo: put empirical p values in a new column and keep the LOD scores
	# this is easy to do here, but need to make sure the rest of the package deals with this change properly
	for (chr.type in chr.types) {

		null.lods <- dplyr::filter(perm.scan.maxes, chrtype == chr.type)
		obs.lods <- dplyr::filter(scanonevar, chrtype == chr.type)

		for (lod.column in lod.columns) {

			evd <- fgev(null.lods[[lod.column]])

			emp.ps <- pgev(q = obs.lods[[lod.column]],
										 loc = fitted(evd)[1],
										 scale = fitted(evd)[2],
										 shape = fitted(evd)[3],
										 lower.tail = FALSE)

			scan.as.emp.ps[scanonevar$chrtype == chr.type, lod.column] <- emp.ps
		}
	}

	# change names to reflect that we now have empirical p values, not LOD scores
  for (lod.column in lod.columns) {
    col.idx <- which(names(scan.as.emp.ps) == lod.column)
    new.name <- paste0('emp.p.', lod.column)
    names(scan.as.emp.ps)[col.idx] <- new.name
  }
	attr(scan.as.emp.ps, 'units') <- 'emp.ps'

	return(scan.as.emp.ps)
}
