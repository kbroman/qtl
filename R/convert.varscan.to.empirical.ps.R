convert.scanonevar.to.empirical.ps <- function(scan, null.scan.maxes) {

  validate.convert.scanonevar.to.emp.ps(scan, null.scan.maxes)

	lod.columns <- grep(pattern = 'lod', names(scan), value = TRUE)
	chr.types <- unique(scan$chrtype)
	scan.as.emp.ps <- scan

  # do each chromosome type and lod column separately
	# empirical p values are put in place of the LOD scores for now
	for (chr.type in chr.types) {

		null.lods <- dplyr::filter(null.scan.maxes, chrtype == chr.type)
		obs.lods <- dplyr::filter(scan, chrtype == chr.type)

		for (lod.column in lod.columns) {

			evd <- fgev(null.lods[[lod.column]])

			emp.ps <- pgev(q = obs.lods[[lod.column]],
										 loc = fitted(evd)[1],
										 scale = fitted(evd)[2],
										 shape = fitted(evd)[3],
										 lower.tail = FALSE)

			scan.as.emp.ps[scan$chrtype == chr.type, lod.column] <- emp.ps
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
