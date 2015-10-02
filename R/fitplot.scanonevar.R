#####################################################################
#
# fitplot.scanonevar.R
#
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
# Contains: fitplot.scanonevar,
#
######################################################################

######################################################################
#
# fitplot.scanonevar:
# Plot to show the fitted model that results from scanonevar() at
# one specific locus of interest.
# Useful to illustrate what data pattern drives a significant LOD score.
#
######################################################################

fitplot.scanonevar <- function(cross, 
                               name.of.marker.to.plot,
                               varscan,
                               sample.based = FALSE,
                               fit.based = !sample.based,
                               phenotype.name,
                               main = NULL,
                               col = 'gray',
                               ...) {

	# store current graphical parameters and customize them for this plot
# 	start.pars <- par(no.readonly = TRUE)
# 	par(mar = c(4,4,4,2))
	set.seed(27599)

  if (sample.based) {
    
    if (missing(phenotype.name) & missing(varscan)) { stop('Need either a varscan or a phenotype name to know what to plot')}
    if (!missing(phenotype.name) & !missing(varscan) & (attr(x = varscan, which = 'pheno') != phenotype.name)) {
      stop('Provided phenotype name and varscan, but they point to different phenotypes')
    }
    
    fitplot.sample.based_(cross = cross, 
                          name.of.marker.to.plot = name.of.marker.to.plot,
                          phenotype.name = phenotype.name,
                          main = main,
                          col = col,
                          ...)
  }
  
  if (fit.based) {
    
    fitplot.fit.based_(cross = cross, 
                       name.of.marker.to.plot = name.of.marker.to.plot,
                       varscan = varscan,
                       main = main,
                       col = col,
                       ...)
  }
	# par(start.pars)
}
