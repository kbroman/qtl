#####################################################################
#
# summary.scanonevar.R
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/qtl package
# Contains: summary.scanonevar,
#
######################################################################

######################################################################
#
# summary.scanonevar:
# Summarize the key features of a scanonevar
#
######################################################################

summary.scanonevar <- function(s, digits = 3) {

	if(!any(class(s) == "scanonevar")) {
		stop("Input should have class \"scanonevar\".")
	}

	# Print Effects Fitted by Null Model
	null.effects <- attr(s, 'null.effects')
	message('Effect estimates when no marker present')
	message('mean covariate effects:')
	print(round(null.effects[[1]], digits))
	message('variance covariate effects:')
	print(round(null.effects[[2]], digits))


	# Print Peak Info
	if (units(s) == 'lods') {
		peaks <- s %>%
			mutate(full.peak = and(lod.full > lag(lod.full, default = 0),
														 lod.full > lead(lod.full, default = 0))) %>%
			mutate(mean.peak = and(lod.mean > lag(lod.mean, default = 0),
														 lod.mean > lead(lod.mean, default = 0))) %>%
			mutate(var.peak = and(lod.var > lag(lod.var, default = 0),
														lod.var > lead(lod.var, default = 0)))

		message('Full Model Peaks:')
		print(dplyr::filter(peaks, full.peak == TRUE, lod.full > 2*mean(s$lod.full)))
		message('Mean Model Peaks:')
		print(dplyr::filter(peaks, mean.peak == TRUE, lod.mean > 2*mean(s$lod.var)))
		message('Var Model Peaks:')
		print(dplyr::filter(peaks, var.peak == TRUE, lod.var > 2*mean(s$lod.var)))
	}

	if (units(s) == 'emp.ps') {
		peaks <- s %>%
			mutate(full.peak = and(emp.p.full < lag(emp.p.full, default = 1),
														 emp.p.full < lead(emp.p.full, default = 1))) %>%
			mutate(mean.peak = and(emp.p.mean < lag(emp.p.mean, default = 1),
														 emp.p.mean < lead(emp.p.mean, default = 1))) %>%
			mutate(var.peak = and(emp.p.var < lag(emp.p.var, default = 1),
														emp.p.var < lead(emp.p.var, default = 1)))

		message('Full Model Peaks:')
		print(dplyr::filter(peaks, full.peak == TRUE, emp.p.full < 0.5*mean(s$emp.p.full)))
		message('Mean Model Peaks:')
		print(dplyr::filter(peaks, mean.peak == TRUE, emp.p.mean < 0.5*mean(s$emp.p.mean)))
		message('Var Model Peaks:')
		print(dplyr::filter(peaks, var.peak == TRUE, emp.p.var < 0.5*mean(s$emp.p.var)))
	}



}
