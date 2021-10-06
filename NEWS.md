Revision history for the R/qtl package
----------------------------------------------------------------------
copyright (c) 2001-2021, Karl W Broman
<https://rqtl.org>

    The R/qtl package is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License,
    version 3, as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but without any warranty; without even the implied warranty of
    merchantability or fitness for a particular purpose.  See the GNU
    General Public License for more details.

    A copy of the GNU General Public License, version 3, is available
    at https://www.r-project.org/Licenses/GPL-3
----------------------------------------------------------------------

## Verison 1.50, 2021-10-06

## New features

- `cim()` now includes an `addcovar` argument for including additional
  covariates in the model.

### Minor changes

- Revised qtlversion() to handle a case like "1.50".

- Added `#define USE_FC_LEN_T` in C code that calls Fortran, because
  of a change in R 3.6.5 that's going to be required in R 4.2.0.


## Version 1.48, 2021-03-24

### Bug fixes

- Fixed bug in `addqtl()` and `addpair()` in which the argument
  `forceXcovar` wasn't getting passed to `scanqtl()`.

- Fixed bug in `stepwiseqtl()` regarding the way the null LOD score is
  calculated.


## Version 1.47, 2021-01-07

### Minor changes

- Added function `find_large_intervals()` for finding inter-marker
  intervals in a map with length greater than some value.

- Fixed potential problem in documentation, since `plot()` has moved
  from the graphics package to base.

- Acknowledge R Core Team among contributors, as zeroin function
  (in C) had been taken from R version 2.15.1. Also add a Copyrights
  field to the DESCRIPTION file.

- Allow `rescalemap()`, `shiftmap()`, `summaryMap()`, and
  `jittermap()` to work with plain lists.

- Fixed Issue #91 where pull.rf() gives a cryptic error if marker
  names are not all distinct.

### Bug fixes

- Fix a problem in `inferredpartitions()` that occurs in the devel
  version of R.

- Small change to `read.cross()` to avoid warning about length of
  `alleles` argument for `crosstype="4way"`. (Fixes Issue #90.)

- Small change to `read.cross()` to avoid messing with X chromosome
  genotypes when `crosstype=4way"`. (Fixes Issue #88.)


## Version 1.46, 2020-02-28

### Minor changes

- In preparation for R 4.0, when the default for `stringsAsFactors`
  will become `FALSE` rather than `TRUE` in `read.table()` and
  `data.frame()`, we needed to add `stringsAsFactors=TRUE` to calls to
  `read.table()` and `data.frame()` in various places. Also there was
  some ugliness regarding `addpair()`.

### Bug fixes

- Fix bug in `mqmpermutation()` by removing the uses of batches and
  the `batchsize` argument.


## Version 1.45, 2020-02-01

### Minor changes

- Added `plot()` and `summary()` functions for the output of
  `comparegeno()`.

- Added internal functions `crosstype()` and `chrtype()`.

- Added argument `crosstype` to internal function `fitqtlengine()`
  rather than taking it from the cross attributes.

- Went through full package to replace use of `class(blah)=="blah"`
  with `inherits(blah, "blah")`.

- mqmscan() and mqmscanall() now give a warning about omitting the X
  chromosome, and give a more meaningful error if there are no autosomes.

- Improved warning and error messages in several places: rather than
  "Chromosome misspecified" say "Chromosome __ not found"

### Bug fixes

- Fix bug regarding missing phenotypes in `stepwiseqtl()`.

- Fix bug in `addpair()` re converting map to data frame (getting
  error like `cannot coerce class "A" to data.frame`).

- Fix bug related to reading 4-way cross data, to ensure that the
  genetic map for each chromosome is a 2-row matrix.

- Fix bug in `refineqtl()` that gave a warning about `min(diff(a))`
  when there was a single marker on a chromosome.
  ([Issue #78](https://github.com/kbroman/qtl/issues/78))

- Added explanations of a couple of arguments for mqmscan() that had
  previously not been explained.


## Version 1.44, 2019-01-22

### Minor changes

- Added col argument to geno.image() for custom colors.

- Revision to plot.scanone() to handle +/- Inf in the LOD score column.

- Add better error message in read.cross with format="mapqtl"

- Revision to summary.map to just give warning if input does not have
  class "cross" or "map".

- Improved error message in scanone() when phenotypes are not numeric.

- Small change in getgenonames() to avoid a problem if there aren't
  enough allele codes.

- plotLodProfile() now gives a more informative error message if
  called with a null QTL model.

- Revised mqmplot.circle() so that chromosome IDs don't need to be
  numbers.

- Fix small bugs in c.cross() and checkcovar().


## Version 1.42, 2018-02-19

### Major changes:

- Removed the functions plot.errorlod, plot.geno, plot.info,
  plot.missing, plot.pheno, plot.pxg, and plot.rf. Changes in the
  process of submitting packages to CRAN have made this necessary.
  Each of the packages has alternative names that have been used in
  the tutorials for some years:

  ```
  plot.errorlod  plotErrorlod
  plot.geno      plotGeno
  plot.info      plotInfo
  plot.missing   plotMissing
  plot.pheno     plotPheno
  plot.pxg       plotPXG
  plot.rf        plotRF
  ```

  If you have old scripts that use these function names, add the
  following code at the top:

  ```
  source("https://rqtl.org/dotfunc.R")
  ```

### Minor changes:

- Fixed a bug in stepwiseqtl() regarding model="binary"; a number of
  cases where I wasn't passing model to other functions, like
  refineqtl().

- Dealt with possibility NA LOD scores in refineqtl().

- Removed questionable use of try() from subset.scanoneperm() and
  subset.scantwoperm().

- Fixed a bug in mqmscanfdr (to permute all chromosomes not just
  chromosome 1).

- Made estimate.map=FALSE the default in read.cross.

- Revised read.cross for the "bcsft" cross type so that it should be
  able to handle sep=";".

- Removed the file vignettes/fancyheadings.sty

- Fixed a couple of minor problems that showed up in byte-compiling
  the package (breaks outside of loops and use of <<-)


## Version 1.41, 2017-06-11

### Major changes:

- None.

### Minor changes:

- Fixed bugs in formMarkerCovar to give better error messages if
  markers aren't found.

- Fixed a bug in reviseXdata in case of intercross with both sexes but
  backward direction only.

- Fixed a bug in mqmscan that gave an error about duplicate row names.

- Fixed a bug in summary.scantwo for the case of autosome/X chromosome
  specific permutations.

- Fix bugs in pull.argmaxgeno and pull.genoprob for 4-way crosses with
  include.pos.info=TRUE. (Need to account for map being a matrix.)


## Version 1.40, 2016-10-31

### Major changes:

- cim() now halts with an error for cross type "4way". The method has
  not been implemented for 4-way crosses, and the results are not
  meaningful.

### Minor changes:

- Small change to the way Bayesian credible intervals are calculated
  by bayesint(), concerning the treatment of widths of intervals
  between positions.

- Fix bug in switchAlleles() so that it works with cross type "bcsft"
  (and will give an appropriate error message for unsupported cross
  types).

- sim.cross gives a warning if model is specified but not used (this
  is the case for RILs, where we've not implemented the simulation of
  QTL effects)

- plot.pxg (aka plotPXG) passes ... to plot(), so now you can control
  the y-axis limits via ylim.

- Fixed a problem with column names of output of scantwopermhk.


## Version 1.39, 2016-03-04

### Major changes:

- Fixed a bug in scantwo() when using weights with multiple phenotypes
  when some of the phenotypes contain missing values. (This is the
  same as the bug fix in scanone in version 1.38. Reminder of
  important principle: when you find a bug, look for other possible
  instances of the same bug.)

- Changed the color scheme for plot.scantwo and plot.rf from red/blue
  rainbow to Viridis (see Option D at https://bids.github.io/colormap/)
  The plot.rf function has a new argument col.scheme; if you want to
  use the old red/blue scheme, use col.scheme="redblue".

- Fixed a bug in cim() that was causing a segmentation fault.

### Minor changes:

- Added a function table2map() for converting a data frame with marker
  positions (chr, position) into a map object.

- Added a bit more detail in the help file for readMWril().


## Version 1.38, 2015-12-03

### Major changes:

- Get rid of everything to do with degrees of freedom in scanone and
  scantwo. The checks seem to offer little value but rather produce
  cryptic warnings that confuse many users.

- Fixed a bug in scanone() when using weights with multiple phenotypes
  when some of the phenotypes contain missing values.

### Minor changes:

- In plotLodProfile, if col has length > 1 but no equal to the number
  of QTL, give a warning, but either repeat or subset the vector
  rather than just using the first value.

- Modified read.cross with format "mm" to handle files with "bc" where
  we usually see "backcross".


## Version 1.37, 2015-08-24

### Major changes:

- None.

### Minor changes:

- In read.cross with format="csv", "mm", or "tidy", don't let it
  reorder the chromosomes (which it would do if there were chromosomes
  named other than numbers < 1000, "X", or "x").

- drop.markers now gives an error if you try to drop *all* of the
  markers.

- If cross object contains no genotype data, totmar() and nmar() now
  give more meaningful errors.

- Fixed a bug in scantwo and scantwopermhk in which an error would
  occur if reduce2grid had been called and assumeCondIndep=FALSE. Now
  forcing assumeCondIndep=TRUE in this case.

- Fixed a bug in plot.pxg, in the case that not all possible genotypes
  were observed at a marker.

- Fixed a bug in stepwiseqtl, where if covar is not a data frame, they
  don't get considered in the model.

- Fixed a bug in fill.geno for method="maxmarginal" (wasn't putting NAs in
  genotypes with probability < min.prob)

- Fixed a bug in refineqtl that arises when multiple QTL are at exactly
  the same position, which can arise in stepwiseqtl.


## Version 1.36, 2015-03-05

### Major changes:

- None.

### Minor changes:

- Added a function flip.order() for flipping the order of markers on
  selected chromosomes.

- Added scanonevar.meanperm and scanonevar.varperm (from Robert Corty)
  for permutation tests with scanonevar().

- Revised plotPheno (aka plot.pheno) so that one can control the
  x-axis label and title (also, in a histogram, the breaks).

- plotPXG: if infer=FALSE and there are no fully-informative genotypes
  (e.g., in a 4-way cross), give a more informative error.

- geno.image: allow control of x- and y-axis labels; allow suppression
  of axes.

- Removed some warnings about missing end-of-line characters, in
  read.cross with MapQTL format.

- Fixed a bug in scanonevar; was failing with an error about coercing
  class "A" to a data.frame

- Dropped the name summary.scantwo.old(); still available as
  summaryScantwoOld().


## Version 1.35, 2014-12-15:

### Major changes:

- None.

### Minor changes:

- Fix an important bug in summary.cross.

- Change a couple of abs() to fabs() in C code.


## Version 1.34, 2014-10-30:

### Major changes:

- Added ability to do X-chr-specific permutations in scantwo (argument
  perm.Xsp, as in scanone). Separate thresholds are obtained for the
  regions A:A, A:X, and X:X regions, maintaining control of the
  overall false positive rates.

- Added a function scantwopermhk that just performs scantwo
  permutations with Haley-Knott regression; faster and with lower
  memory usage than scantwo.

- With X-chr-specific scantwo permutations, calc.penalties will give
  separate main effect for autosome and X chromosome, and separate
  interaction penalties for A:A, A:X, X:X. For A:A interactions, we
  still use "light" and "heavy" penalties; for A:X and X:X
  interactions, only the "heavy" penalty is used. These penalties may
  be used in stepwiseqtl for better treatment of the X chromosome in
  automated inference of multi-QTL models.

- Added scanonevar() function, for a single-QTL genome scan for QTL
  affecting not just the mean phenotype but also the variance. (Code
  from Lars Ronnegard; method in Ronnegard and Valdar Genetics
  188:435-447, 2011.)

- Added "tidy" format to read.cross and write.cross. This separates
  the data into three comma-delimited files, for genotypes,
  phenotypes, and the marker map. Separating the data in this way
  allows each file to be in a simpler format.

### Minor changes:

- Add another option to fill.geno: impute using maximum marginal
  probability.

- Add function map2table (output like pull.map with as.table=TRUE, but
  starting with a map rather than with a cross).

- Fixed a bug in est_map_ri8self.c (thanks to Rohan Shah)

- Fixed a bug in scanone/scantwo stratified permutations in batch,
  with multiple phenotypes, some with missing values (thanks to
  John Lovell).

- Fixed some bugs in read.cross with format="mapqtl".


## Version 1.33, 2014-08-12:

### Major changes:

- Can read/write 4-way cross data in MapQTL format (thanks to Timothee
  Flutre).

### Minor changes:

- read.cross with format="qtlcart" can now read doubled haploids
  (dh/Ri0).

- Fixed potential problem in read.cross with format="csv" or "csvs" when
  there are many empty cells in the phenotype data.

- Fixed bug in read.cross for format="csv" that shows up in some rare
  cases: markers not ordered by linkage group, no positions provided,
  and chromosome IDs non-numeric. It was a pretty bad bug, as marker
  genotypes got scrambled.

- Fixed some memory leaks in MQM code.


## Version 1.32, 2014-05-28:

### Major changes:

- None.

### Minor changes:

- fitqtl with model="normal" now returns residuals as an attribute.

- Added an additional argument to plot.scanone, bgrect, for making the
  background of the plotting region a different color.

- Revised cleanGeno to work with any cross having two possible
  genotypes (i.e., not just bc but also riself, risib, dh, haploid).

- Revised summary.cross so that overall genotype frequencies are given
  separately for autosomes and the X chromosome.

- Fixed typo in a warning in add.threshold.

- Fixed a bug in reduce2grid, regarding format of attributes

- Fixed a bug in MQM: in some circumstances, the last marker was
  always included as cofactor; other cleanup in MQM code.


## Version 1.31, 3/19/2014:

### Major changes:

- Added a function reduce2grid() for subsetting the genotype
  probabilities (from calc.genoprob) or imputations (from sim.geno) to
  the evenly spaced grid of pseudomarkers, for use in the case of very
  dense markers when it is too computationally intensive to perform a
  genome scan at both the markers and the grid of pseudomarkers.

### Minor changes:

- Fixed a problem in write.cross with format="csv" (for cross types
  other than BC and F2, it would use incorrect genotype codes if there
  was no "allele" attribute for the cross).

- Add a couple of arguments to plotLodProfile: showallchr, to show all
  chromosomes (and not just those with QTL), and textsep, to control
  the separation between the QTL labels and the LOD curves.

- Fix a bug in bcsft.c, regarding potentially over-running an array

- Fix problem with plotLodProfile when it's maximized at multiple
  locations.

- Fix problem in refineqtl and stepwiseqtl; map attribute in qtl
  object would get unintentionally subsetted if the qtl object needed
  to be re-created.

- In addqtl, addint, and addcovarint, have require.fullrank=FALSE be
  the default; require.fullrank=TRUE remains the default in
  stepwiseqtl.

- Fixed bug in which summary.scantwo was re-ordering the chromosome
  factor levels.


## Version 1.30, 2/10/2014:

### Major changes:

- Revised parallel code (in scanone, scantwo, mqmpermutation, and
  mqmscanall) to use a different function for Windows, as mclapply()
  doesn't work there.

### Minor changes:

- Fixed problem in formLinkageGroups when used with results of
  markerlrt(): no recombination fractions so use max.rf=Inf

- Added arguments type, cex, pch, and bg to plot.scanone, to be passed
  to lines() in making the plot. Thus, you can use type="p" to get
  points only.

- Fixed a bug in scanqtl that showed up if there were missing
  phenotypes and no covariates.


## Version 1.29, 12/9/2013:

### Major changes:

- None.

### Minor changes:

- Slight change to summary.qtl to deal with QTL objects with no QTL.

- For the QTLRel package, now "export" reviseXdata() in NAMESPACE.


## Version 1.28, 9/23/2013:

### Major changes:

- Added cross type "haploid". Like backcross ("bc") or doubled
  haploids ("dh") but with genotype labels like "A" and "B" instead of
  "AA" and "AB"/"BB".

- Added crosstype argument to read.cross, to force a particular cross
  type (such as "riself").

- For parallel processing, replaced reliance on the snow library with
  the use of the parallel library.

### Minor changes:

- Added formMarkerCovar(), to facilitate use of markers as covariates
  in QTL analysis.

- Added function addmarker() for adding genotypes for a marker to a
  cross object.

- Added function nqtl() for counting number of QTL in QTL object.

- Added examples to help files for plotPXG and effectplot on getting the
  output.

- Slight change to way to handle random number generation for
  cluster-based computing.

- Fix bug in fitqtl-link functions for model="binary" (re matrix rank)

- Fixed write.cross to allow use for BCsFt crosses.

- Fixed an out-of-bounds error in the C++ code mqmscan.cpp.

- In stepwiseqtl with verbose=FALSE, the initial LOD score is no
  longer printed.


## Version 1.27, 4/11/2013:

### Major changes:

- Implemented HMM algorithms for advanced backcross/intercross.

### Minor changes:

- Improvement in addqtl, addint, scanqtl, stepwiseqtl
  to better handle collinearity in the design matrix that could give
  spurious evidence for QTL in large QTL models.  This can still be a
  problem with addpair (and stepwiseqtl with scan.pairs=TRUE).

- Slight change to help file for read.cross, to be more explicit about
  the "csvs" format.

- Fixed a bug in read.cross for format="qtlcart" with RIL data.

- Made cross type in "qtlcart" files case insensitive in read.cross.
  For example, any of Ri1, RI1, rI1, or ri1 will be treated the same.


## Version 1.26, 11/27/2012:

### Major changes:

- We changed the treatment of the ind argument in subset.cross to
  avoid problems that occurred when the cross contained numeric
  individual IDs.  Now, we match the values in ind against the
  individual IDs only if ind values are character strings.  If the ind
  values are numeric, we treat them as numeric indices and they are
  not compared against the individual IDs.

### Minor changes:

- In plot.rfmatrix, include marker name in title (unless main is
  provided as an argument).

- Add warning message to find.marker if there's no match for a given
  chromosome name.

- Slight change in find.markerpos(), to speed it up.

- Slight change to locateXO to save genotypes to left and right of
  each crossover.

- Small addition to the "A shorter tutorial of R/qtl" (rqtltour2.pdf).

- Fixed slight bug in summary.scanone for format="tabByChr" or
  format="tabByCol".

- Fixed bug in fitqtl for case that individuals have missing
  phenotypes or covariates and there's a QTL on the X chromosome

- Fixed bug in effectplot, regarding coding of genotypes on the X
  chromosome.

- Fixed bug in mqmscan for case when estimate.map=TRUE but plot=FALSE.

- Fixed bug in c.cross for case that there are different sets of markers.

- Fixed bug in xaxisloc.scanone for case that chromosomes don't start at 0.

- Fixed bug in locateXO; gave core dump if there was just one marker on
  the chromosome.

- Fixed bug in scanone and scantwo for case that weights are used but
  there are individuals with missing phenotype; the weights weren't
  being subsetted appropriately.

- Revised example help file for multitrait data to use mqmscanall
  rather than scanall, since the latter function has no help file.

- scanPhyloQTL now gives warning if there are different marker maps.

- Fixed bug in mqmaugment; with one phenotype, its name was getting
  changed to "pheno".


## Version 1.25, 8/13/2012:

### Major changes:

- None.

### Minor changes:

- Revised write.cross to output data in "qtab" format; see
  https://github.com/qtlHD/qtlHD/blob/master/doc/input/qtab.md

- Revised summary.scanone to be more consistent about only picking one
  row per chromosome even when there are multiple positions sharing
  the maximum LOD on that chromosome.

- Fixed bug in stepwiseqtl; backward deletion steps were not dealing
  with the drop-one-qtl results from fitqtl appropriately.

- Revised fitqtl to include formulas and lod scores as attribute in
  drop-one-qtl analysis.

- Fixed a bug in fitqtl regarding the adjustment for the X
  chromosome.  The problem shows up in stepwiseqtl; if X chr enters
  model and is then removed, the covariates adjusting for sex and pgm
  continue to be used.  Added argument to refineqtl, fitqtl, scanqtl,
  addqtl, addpair, to force X-related covariates into model.

- In stepwiseqtl, include penalties as attribute in the output.

- Slight change to checkcovar(): omit individuals with
  phenotypes/covariates that are +/- Inf, as well as those that are
  missing.

- Handle missing values in mf.stahl and imf.stahl.

- Fixed rare bug in fitqtl regarding X chr loci with interactions.

- Fixed bug in sim.cross for type="bc" that resulted in loss of "X"
  chromosome type.

- Speed up some of the examples in the help file, so that R CMD check
  doesn't take so long.


## Version 1.24, 5/25/2012:

### Major changes:

- Fixed a major bug in checkcovar, used by scanone and scantwo to omit
  individuals with missing phenotypes.  If there is an "ID" column
  that is numeric, the wrong individuals could be omitted, and
  genotypes and phenotypes would be misaligned.

- Changed the names of a number of functions, in order to avoid the
  "Note" in R CMD check, and because Prof. Brian Ripley asked me to.

   ```
   plot.map ->      plotMap
   plot.missing ->  plotMissing
   plot.errorlod -> plotErrorlod
   plot.geno ->     plotGeno
   plot.info ->     plotInfo
   plot.pheno ->    plotPheno
   plot.pxg ->      plotPXG
   plot.rf ->       plotRF
   summary.map ->   summaryMap
   summary.scantwo.old -> summaryScantwoOld
   ```

- Revised tutorials to use the new naming scheme.

- Revised the emission probabilities for dominant markers in an F2,
  for the HMM calculations.  Previously, we had
      Pr(O = not A | g = AA) = Pr(O = not B | g = BB) = epsilon/2
  these have been changed to
      Pr(O = not A | g = AA) = Pr(O = not B | g = BB) = epsilon
  This corresponds to "not A" being "H or B"; similarly for not B.
  Results will change only for an intercross with dominant markers,
  and generally only slightly.

### Minor changes:

- Changes to the format of the output of summary.scanPhyloQTL for
  format="lod".  The final column is now the maximum LOD score across
  partitions; the difference between the maximum and the
  second-highest is now third-to-last; the threshold argument is
  applied to the overall maximum rather than to that difference.

- Revisions to read.cross (for the csv formats) from Steffen Moeller,
  to give some more informative warnings/errors.

- Fixed a bug in scantwo permutations in case that a chromosome has
  multiple markers but they span < step from calc.genoprob.

- Fixed a bug in interpPositions; problems if the input had missing
  rownames.

- Renamed the README.txt file as INSTALL_ME.txt; added a new
  README.txt that provides a brief description of the package.


## Version 1.23, 3/6/2012:

### Major changes:

- None.

### Minor changes:

- Added functions pull.genoprob, pull.argmaxgeno, and pull.draws, to
  pull out those bits from a cross object as a single big matrix or
  array.

- Added a function inferFounderHap() for crudely inferring founder
  haplotypes in multi-parent RIL, using groups of adjacent markers.

- Added function nullmarkers() for identifying markers with no
  genotype data.

- Revised sim.cross so that founder genotypes are included in output,
  for 4- and 8-way RIL.

- Added HMM functions to handle a special design for MAGIC lines from
  BioGemma (http://www.biogemma.fr/indexuk.php).

- Revised subset.cross and clean.cross so that cross information in 4-
  and 8-way RIL don't get lost.

- Revised plot.pxg, effectplot and effectscan to give a more
  informative error if the selected phenotype is not numeric.

- Revised qtlversion() to use packageVersion().

- Fixed bug in summary.map: class included the function data.frame;
  not just the character string "data.frame".

- Revised various utility functions to retain the "onlylod" attribute
  in cross$rf, if it's there.


## Version 1.22, 11/28/2011:

### Major changes:

- Revised plot.map to deal with a pair of maps with markers in
  different orders (or with some markers appearing in one map and not
  the other).  We still require that the two maps have the same
  chromosomes and chromosome names (with chromosomes in the same
  order).

- Revised scantwo to allow analysis of individual chromosome pairs,
  and reorganized the way that scantwo permutations are done (first
  summarizing each chromosome pair and then overall); this should
  eliminate the memory problems we've had with scantwo permutations.

### Minor changes:

- Added warning to help file for fitqtl() regarding the estimated
  percent variance explained in the case of linked loci: the values
  are misleading.

- Fixed problem in comparecrosses() regarding Inf/-Inf phenotypes.

- Fixed a bug in scantwo, method="hk", for multiple phenotypes or
  permutation tests in batch.  This showed up only when there were
  missing genotypes at one or both putative QTL.

- write.cross with format="csv" now exports genotypes as AA/BB for
  RIL (previously, genotypes were written as AA/AB).

- Fixed bug in addint() and addcovarint()

- Revised typingGap to have an argument 'terminal'; if TRUE, just look
  at the gap from the terminal markers to the first typed interior
  marker, giving 0 if the terminal markers are both typed.

- Fixed bug in calc.genoprob with stepwidth="max" in the case that no
  pseudomarkers are to be added.

- Fixed bug in geno.table for the case that X chromosome has just one
  observed genotype.

- No longer allow "" in na.strings in read.cross for csv files.

- Revised all calls to data.frame() and as.data.frame() to override
  global option of stringsAsFactors, so that we know what's going to
  happen.

- Revised scantwo so that if the 'chr' argument is a list, we do just
  the scans of the chr in the first component against those in the
  second component.

- Slight change in plot.info to deal with inclusion of 'main'
  argument.

- Added NAMESPACE file.

- Slight changes to avoid some R warning messages.

- Slight change to imf.cf to give more accurate results.

- Fixed warning message in replacemap.

- Fixed warning message in +.scanone.

- Fixed bug in mqmfind.marker.

- Fixed slight bugs in print.addint and print.addcovarint.

- Removed a bunch of unused variables from C code.


## Version 1.21, 3/21/2011:

### Major changes:

- None.

### Minor changes:

- Add "addchr" argument to find.pseudomarker.  The default is TRUE,
  and returned non-marker locations have names like "c5.loc25" (as in
  the output of scanone).  If FALSE, that initial "c5." part is left
  off, to return just strings like "loc25" (as in the genotype
  probabilities from calc.genoprob).

- Revise calls to rainbow() in plot.rf and plot.scantwo so that they
  no longer use the 'gamma' argument, which is being removed from
  future versions of R.

- Slight change to format of verbose output in est.map with m>0 (that
  is, under interference).


## Version 1.20, 2/18/2011:

### Major changes:

- Enabled est.map to use multiple processors via snow; added argument
  n.cluster to indicate the number of cluster nodes to use.

- Added the option stepwidth="max" in calc.genoprob, sim.geno, and
  argmax.geno.  This inserts the minimal number of pseudomarkers so
  that the maximum step between points is as indicated by the "step"
  argument.

### Minor changes:

- Fixed a bug in refineqtl() that kept it from working for a 4-way
  cross.  (The bug also broke stepwiseqtl().)

- Fixed a problem in the internal function dropXcol() that led to a
  crash in scantwo() for 4-way crosses with an X chromosome.

- Fixed a bug in mqmscan() regarding the chromosome names in the
  output.

- Trap cases of X chromosome for crosses other than bc/f2 in
  stepwiseqtl, makeqtl, addqtl, and addpair.

- Added a function phenames() for pulling out the names of the
  phenotypes.

- Small revisions and enhancements to some of the MQM plots.

- Revised subset.cross so that if the cross contains QTL genotypes
  (from sim.cross), these are also subsetted.

- Revised replacemap.cross (aka replace.map) so that it will also
  replace the maps in results from calc.genoprob, sim.geno and
  argmax.geno, using interpolation if necessary.

- Fixed a couple of minor bugs in mqmscan: one giving duplicate row
  names, another resulting in pseudomarkers outside the terminal
  markers even when off.end=0.

- Fixed bugs in scanone and scantwo regarding batch mode for
  model != "normal".

- Fixed a bug in refineqtl; it was including all possible covariates
  refered to in the data frame 'covar', even if they weren't referred
  to in the formula.

- Fixed a bug in plot.geno, introduced in version 1.19, that made it
  not work with horizontal=FALSE.


## Version 1.19, 11/29/2010:

### Major changes:

- Added a tutorial on genetic map construction; find it within the
  package at docs/geneticmaps.pdf, or (more simply, probably), find it
  on the web at https://rqtl.org/tutorials/geneticmaps.pdf

- Added two additional formats to summary.scanone(), "tabByCol" and
  "tabByChr".  These produce tables of LOD peaks organized by LOD
  score column or by chromosome.  The tables include approximate
  confidence intervals for QTL location (as calculated by the lodint()
  or bayesint() function; which one is indicated by the new argument
  ci.function).  Here's an example:

  ```
  bp:
             chr  pos ci.low ci.high  lod  pval
  c1.loc44.5   1 47.8   35.5    85.0 3.56 0.007
  D4Mit164     4 29.5   18.8    30.6 8.09 0.000

  sqrt:
             chr  pos ci.low ci.high  lod  pval
  c1.loc44.5   1 47.8   35.5    84.8 3.63 0.007
  D4Mit164     4 29.5   17.2    30.6 8.08 0.000
  ```

### Minor changes:

- Revised summary.scanoneperm to include a new argument,
  controlAcrossCol.  If TRUE, LOD thresholds will control error rate
  not just across the genome but also across the LOD score columns.

- Added function droponemarker, with the aim of identifying
  problematic markers by dropping one marker at a time and calculating
  a LOD score and a change in the estimated genetic length of the
  respective chromosome.

- Added a function pull.rf for pulling out either the estimated
  recombination fractions or the lod scores, as calculated by
  est.rf(), from a cross object.  Also added a function plot.rfmatrix,
  for plotting a slice through these.

- Added a function cleanGeno for removing genotypes that are possibly
  in error (as indicated by apparent tight double-crossovers).  Use
  this function with caution.

- Added a function typingGap, which calculates, for each individual
  and each chromosome, the maximum distance between typed markers.

- Revised MQMscan so that the output contains covariate information,
  to be plotted with add.cim.covar().

- Revised c.scanone(), c.scanoneperm, c.scantwo, and c.scantwoperm so
  that the input "..." can be a list of
  scanone/scanoneperm/scantwo/scantwoperm objects.

- Revised plot.geno() so that, for the X chromosome in a backcross or
  intercross, the genotypes appear appropriately, with females being
  homozygous or heterozygous and the males hemizygous, though we plot
  the hemizygous genotypes as if they were homozygotes.

- Revised subset.scanoneperm and subset.scantwoperm so that one may
  pull out a subset of replicates (not just columns).  Added functions
  [.scanoneperm and [.scantwoperm so that one can use [] to subset.

- Revised locateXO() so that the output contains a column with the
  number of typed markers between adjacent crossovers.

- In orderMarkers(), if verbose is numeric and > 1, even more
  information on the progress of the calculations is provided.

- Fixed a problem in subset.cross where, in the case that a cross
  contained numeric IDs, subsetting the individuals resulted in them
  being sorted according to their IDs.

- Added a manual page for the function getid(), which was not previously
  documented.  (It is used internally a great deal, and it may be
  useful more generally.)

- Revised effectplot so that if mname2 or mark2 are given but mname1
  and mark1 are not, the arguments get switched.

- Fixed a bug in shiftmap().

- Added an argument, force, to reviseXdata(), to force a change in the
  genotypes; this is for use within plot.geno().

- Slight change to top.errorlod, so that the output columns are not
  factors but character strings.

- Slight change in plot.geno() so that we use filled circles (pch=23) rather
  than calling points() twice for each point.

- Small changes to mqmplot.circle().

- Subtle changes to tutorials new_multiqtl.pdf,
  new_summary_scanone.pdf, new_summary_scantwo.pdf; added .R files
  with the code for these tutorials.

- Fixed slight bug in getid().


## Version 1.18, 8/18/2010:

### Major changes:

- None.

### Minor changes:

- Revised geno.table so that, for 4-way crosses, it gives P-values in
  most cases.  (Previously, it just did so for fully informative
  markers.)

- Changed the default format for max.scanPhyloQTL and
  summary.scanPhyloQTL from format="lod" to format="postprob".

- Added function inferredpartitions for pulling out the inferred
  partitions for a specified chromosome from the output of
  scanPhyloQTL.

- Fixed a problem in ripple with method="likelihood".  (The problem
  arose in revisions in version 1.17.)

- Fixed a problem in est.map that resulted in NAs.  (The problem arose
  in revisions in version 1.17.)

- Slight changes to "inverse" map functions imf.k, imf.h, imf.cf, so
  that recombination fractions >= 0.5 return large map distances
  rather than NAs.

- Fixed a bug in summary.scanPhyloQTL so that it works when the input
  has just one chromosome and so the colnames remain like "AB|CD"
  rather than getting converted to "AB.CD".


## Version 1.17, 7/28/2010:

### Major changes:

- Implemented the fit of models for binary traits in fitqtl(), by
  Haley-Knott regression and multiple imputation.  This model can also
  be used in refineqtl, scanqtl, addqtl, addpair, addint, stepwiseqtl,
  and addcovarint.

- Implemented Haley-Knott regression for binary traits in scanone and
  scantwo.

- In mqmscan(), replaced the arguments step.min and step.max with
  a single argument, off.end (to be more like scanone).

- Added functions for the joint analysis of multiple crosses, in order
  to map QTL to a phylogenetic tree.  (A paper describing the methods
  is in preparation.)  The key function is scanPhyloQTL.  simPhyloQTL
  is used to simulate data.  plot.scanPhyloQTL, max.scanPhyloQTL, and
  summary.scanPhyloQTL are the plot, max, and summary functions for
  the output from scanPhyloQTL.

### Minor changes:

- Added 'offset' argument to est.map, which defines the starting position
  for each chromosome.  If missing, we use the starting positions that
  are currently present in the input cross object.

- Added function shiftmap, for shifting the starting points of a
  genetic map (in a map or cross object).

- Added function switchAlleles, for switching the alleles at selected
  markers in a cross object.  (For example, in a backcross, switching
  AA and AB at a marker; in an intercross, exchanging AA for BB.)

- Revised geno.table to have an additional argument, scanone.output.
  If scanone.output=TRUE, the output is as produced by scanone(), so
  that one may use plot.scanone() to plot the results.

- Added the ability to have ripple() run in parallel, if the snow
  package is installed.  The added argument n.cluster indicates the
  number of parallel nodes to use.  This is reapply only useful with
  method="likelihood"; with method="countxo", it can be slower than
  just using one CPU.

- Added a function pull.markers, which is the opposite of
  drop.markers.

- Added a function drop.dupmarkers, for dropping markers with
  duplicate names.

- Added an argument 'bandcol' in plot.scanone, to specify a color for
  alternating bands to indicate chromosomes.  The default
  (bandcol=NULL) is to not plot such bands.  A good choice might be
  bandcol="gray70".

- Added an argument 'chr' in est.map, to estimate maps for just a
  subset of chromosomes.

- Revised replace.map (and replacemap.cross) so that the map can have
  just a subset of the chromosomes in the cross, in which case only
  the maps for that subset are replaced.

- Changed the default in plot.info from method="both" to
  method="entropy".  Also, added an argument "fourwaycross", so that
  one can look at the missing information just for the alleles of the
  first parent (A vs B) or the second parent (C vs D).  Finally, added
  an argument "include.genofreq"; if TRUE, the results will include
  estimated genotype frequencies at each position.

- Added argument 'simple' to summary.fitqtl, addint, and addcovarint;
  if TRUE, output includes neither p-values nor sums of squares.

- Added a function allchrsplits(), for testing whether to split a
  linkage group/chromosome into two, by calculating a LOD score for
  each interval comparing the full linkage group to the split into two
  groups at that interval.

- Added a function nqrank() for transforming a numeric vector into the
  corresponding normal quantiles.

- Fixed a problem in MQM code regarding negative or really large
  marker positions.

- Fixed a bug in mqmpermutation(), where it was using the wrong
  phenotype name in the output, if pheno.col is something other than
  1.

- Revised est.map for 4-way crosses slightly.  We'd previously
  randomized the map before starting EM, which seemed a bad idea.

- Fixed a bug in scantwo() for the case that multiple phenotype
  columns are considered but they have different patterns of missing
  data.  An error occurred due to a number of stupid small mistakes.
  [Thanks to Ricardo Verdugo for reporting the problem.]

- Revised summary.cross so that if markers are at the same location,
  the warning message indicates which chromosomes are involved.

- Revised summary.scanone so that if there is one LOD score column in
  the output, but the permutation results have more than one, then
  the first column in the permutation results is used and the others
  are ignored (and a warning, rather than an error, is issued).

- Revised calc.plod() to allow penalties to be infinite.

- Revised read.cross with format="csv" to give an error if the 2nd row
  has all blanks.

- Revised plot.scanoneboot, plot.scanoneperm, and plot.scantwoperm so
  that they can use the ... argument more flexibly.  [I use a scheme
  suggested by Brian Yandell.]

- Revised plot.geno so that the ... argument can include xlim and
  ylim.

- Fixed bug in convert2sa for case of chromosome with just 1 or 2
  markers.

- Fixed bug in refineqtl that resulted in sometimes the rownames in
  the lod profile attribute being messed up.

- Revised scanone and scantwo to trap error if perm.strata is not of
  the correct length.

- Added argument 'ind.noqtl' to scanone, which indicates individuals
  that are assigned no QTL effect.  This is for rare (largely
  internal) use for the case that one is combining multiple crosses.

- Fix bug in internal function create.map, for case of a sex-specific
  map with very small length.


## Version 1.16, 5/23/2010:

### Major changes:

- None.

### Minor changes:

- Revised fitqtl() so that estimated QTL effects in RIL or doubled
  haploid are as they should be (half the difference between the two
  homoyzogotes).

- Added a function rescalemap for rescaling a genetic map (as for the
  case that a cross object has marker positions in basepairs and one
  wishes to convert them to Mbp or some approximation of cM
  locations).

- Added a warning message to summary.cross for the case that there are
  chromosomes > 1000 cM in length (which might indicate that they're
  really in basepairs).

- Revised pull.map to have argument as.table; if as.table=TRUE the map
  is returned as a simple table with chromosome assignments and
  positions.

- Revised fill.geno to have a third option, method="no_dbl_XO", which
  fills in missing genotypes between markers with exactly the same
  genotype.

- Revised est.map so that the output with verbose=TRUE is less verbose
  and more informative.

- Changes in C++ code, to fix problems that prevented the package from
  being compiled on Solaris.


## Version 1.15, 5/2/2010:

### Major changes:

- In collaboration with Danny Arends, Pjotr Prins, and Ritsert Jansen,
  we have incorporated Ritsert Jansen's MQM mapping software within
  R/qtl. (Previously, it was available only through the commercial
  software package, mapqtl).  See the tutorial at
  https://rqtl.org/tutorials/MQM-tour.pdf

### Minor changes:

- Added a function transformPheno for transforming one or more
  phenotypes in a cross object.

- Added a function convert.map for converting a genetic map from one
  map function to another.

- Added functions convert2riself and convert2risib, for converting a
  cross object to be treated as RIL by selfing or sib mating,
  respectively.

- Added function simulateMissingData, for omitting genotypes at random
  from a cross object.

- Revised write.cross so that it can handle non-numeric phenotypes.
  Also changed the default for the "digits" argument to NULL, so that
  phenotypes and map positions are not rounded.

- Revised read.cross to have arguments error.prob and map.function, to
  be used if est.map is called.

- Revised locateXO: it no longer shifts the first marker to position
  0, and it has a new argument, full.info; if this is TRUE,
  the output includes not just the estimated crossover locations but
  also the endpoints of the intervals to which they are known to reside.

- Revised summary.cross so that if there are > 30 phenotypes, we don't
  show the percent missing phenotypes for all traits but just overall.

- Revised stepwiseqtl so that if additive.only=TRUE, you only need to
  give one penalty.

- Revised read.cross with format="csv" so that initial fields in 2nd
  row need not be completely empty, but can have white space.  (This
  was a common problem for users importing csv files.)

- In plot.qtl, added argument "justdots", so that one can plot just
  dots at the QTL rather than arrows and QTL names, and "col", the
  color used to indicate the QTL.

- Fixed a bug in xaxisloc.scanone for the case of multiple chr/pos in
  the input.

- Fixed an apparent bug in read.cross for format="qtlcart"; need to
  treat negative genotype codes as missing.

- Fixed a bug in sim.cross for type="4way"; previously it was just
  using the female map for both female and male meioses.

- Fixed a bug in find.pseudomarkerpos in the case of sex-specific maps
  (as in a 4-way cross).

- lodint and bayesint now stop with an error if chr or qtl.index have
  length > 1.  Before, they gave weird results and a meaningless
  error.

- In read.cross with format="csv", give a better error message if
  there are odd values in the marker positions.

- Fixed a bug in scanone/scantwo for RIL on the X chromosome; the X
  chromosome needs to be treated like to autosomes.  We can't really
  deal with the sexes properly here.

- Revised scanone and scantwo so that, if using multiple CPUs via snow
  and calculations are stopped early, the cluster nodes are stopped on
  exit.

- Revised summary.cross so that, for cross type "riself", there's no
  warning about a chromosome named "X" having class "A".

- Removed some odd erroneous code from the plot.qtl function.

- Fixed a bug in summary.scanone that shows up in the case of multiple
  LOD columns with permutations and format="allpeaks".

- Fixed a bug in scanone permutations in the case of multiple
  phenotypes with missing data.

- Fixed a problem in effectplot that arose from an apparent bug in
  weighted.mean in R version 2.10.1


## Version 1.14, 9/30/2009:

### Major changes:

- None.

### Minor changes:

- Slight changes to read.cross for formats "mm" and "qtlcart",
  regarding use of the function grep().


## Version 1.13, 9/10/2009:

### Major changes:

- Fixed a bug in fitqtl regarding the case of a locus on the
  chromosome...the cited LOD scores for X chromosome loci were not
  right.  Had to make slight modifications to scanqtl, addqtl,
  addpair, addint, addcovarint, refineqtl.

### Minor changes:

- Revised makeqtl, addtoqtl, dropfromqtl, replaceqtl, and reorderqtl
  so that qtl objects now include a "chrtype" component, indicating
  whether QTL are autosomal or for the X chromosome.

- Revised scanone and scantwo for multiple phenotypes with missing
  data so that, with method="hk" or method="imp", the phenotypes are
  grouped into batches with matching patterns of missing data, rather
  than just doing each one at a time.

- Added a function locateXO (formerly the internal function locate.xo)
  to estimate the locations of crossovers on a given chromosome.

- Added a function, c.scantwo (aka cbind.scantwo), for concatenating
  multiple scantwo results.

- Revised summary.map to deal with the case of a single marker on a
  chromosome.

- Revised scantwo so that phenotype names are in dimnames of the lod
  component of the output.

- Fixed bugs in read.cross with format="mm" to deal with changes to
  grep.  [replaced grep("^*", ...) with grep("^\\*", ...)]

- Fixed a bug in reorderqtl; the n.gen component was not getting
  fixed.

- Fixed a bug in switch.order; results of est.rf were getting messed
  up.

- Fixed some minor issues regarding hyperlinks in help files.


## Version 1.12, 6/30/2009:

### Major changes:

- Added the ability to simulate RIL and multiple-strain RIL in the
  sim.cross function.  For the case of multiple-strain RIL, one needs
  genotype data on the founder strains, which may be simulated with
  the new function simFounderSnps.  The encoding of genotypes in
  multiple-strain RIL is quite complicated.  See the help file for
  sim.cross.

- Added the ability to deal with 4- and 8-way RIL in calc.genoprob,
  sim.geno, argmax.geno, est.map, ripple, est.rf, tryallpositions,
  calc.pairprob, and calc.errorlod.

- Added a function, readMWril, for reading data on 4- or 8-way RIL.

- Fixed a minor problem in read.cross, with format="csv", in the case
  of many phenotypes that resulted in *really* slow data import.  (To
  read a file with 200 individuals and 1500 phenotypes, it would take
  about 60 seconds and now takes about 2 seconds.)

- Added the ability to have scanone and scantwo permutations run in
  parallel, if the snow package is installed.  The added argument
  n.cluster indicates the number of parallel nodes to use.

### Minor changes:

- Added a CITATION file; type citation("qtl") within R to get
  information on the citation to use in articles that make use of
  R/qtl.

- You can now subset crosses with brackets, [ ], as with a matrix with
  rows=chromosomes and columns=individuals.  See the examples in the
  help file for subset.cross.

- Added a utility function, findDupMarkers, for identifying groups of
  markers with identical genotype data.  (This is useful for reducing
  the genotype data in the case of a very high marker density.)

- Added a utility function, xaxisloc.scanone, for finding x-axis
  locations for given genomic positions in a plot of scanone results
  (useful for adding annotations, such as text or arrows).

- Added functions subset.map and `[.map` for pulling out selected
  chromosomes from a map object.

- The output of tryallpositions, for testing possible positions for a
  genetic marker, now has class "scanone", so that one may use
  plot.scanone, summary.scanone, etc.

- Added an argument mark.diagonal to plot.rf(), to include black lines
  segments around the pixels on the diagonal.  This helps to separate
  the upper left triangle from the lower right triangle.
  (The default is FALSE.)

- In geno.crosstab, the first argument (mname1) can now be a vector
  with the two marker names; in this case, mname2 should be missing.

- Revised clean.scantwo so that, by default, positions must have at
  least one marker *in between* them.  Added arguments n.mar
  (no. markers that must separate two positions) and distance (cM
  distance between two positions).  These arguments were also added to
  scantwo (as clean.nmar and clean.distance).

- Revised summary.scanone and summary.scantwo so that the perms
  argument can contain a single column of permutation results, in
  which case they are assumed to apply to any  LOD score columns.

- Revised summary.scanone so that the perms argument can be
  scantwo permutations results; added an internal utility function,
  scantwoperm2scanoneperm, for pulling the scanone permutations out of
  the scantwo permutations.

- In plots of output from addpair with a special formula, plot.scantwo
  now gives just one set of numbers on the color scale.

- Fixed a bug in plot.scanone that shows up if the "chr" column is not
  a factor.  Now we convert the column to a factor in advance.

- Fixed a bug in stepwiseqtl in the case that the inferred model
  contains no QTL; deparseQTLformula needs to deal with the NULL case.

- Fixed a bug in c.cross regarding "map" attributes in $prob or $draws.

- Fixed a bug in cim() [reported by Sandy Taylor] that occurred in the
  case of multiple marker covariates within a window (and >3 marker
  covariates on that chromosome.

- Fixed a bug in summary.scantwo in case that scanoneX component is
  numeric(0) and not NULL, which resulted in a major crash.

- Revised +.scantwo and -.scantwo so that if the scanoneX component in
  the input is NULL, the output has scanoneX that is NULL.

- Fixed a bug in addpair() in the case that the user gives a formula
  that includes one but not both of the new QTL.

- Revised ripple() to give a warning message if the chr argument is
  not provided.

- In compareorder(), switch.order(), and ripple(), changed the default
  value for the tol argument to 1e-6 (as in est.map).

- Slight change in compareorder() so that the order argument can be of
  length n.mar+2, but with only the first n.mar items considered.

- Slight change in checks of chr argument in ripple and switch.order;
  added a function testchr for checking that a chromosome argument is
  okay.

- In the output of c.cross, there is a numeric phenotype "cross" that
  indicates which individuals are from which cross, as a single
  column.

- Revised fixXgeno.f2 so that warnings are given the appropriate
  allele labels (if they are provided to read.cross).

- Fixed a few problems in c.cross regarding the attributes to
  QTL genotype probabilities and imputated genotypes

- Added a function chrnames, for pulling out the chromosome names from
  a cross.

- Revised the software license to the GNU General Public License,
  version 3.

- Minor change in plot.cross so that if one types plot(mycross, mymap)
  it is shipped to plot.map rather than giving an error message.

- Revised the R/qtl tutorial to refer specifically to the GPL v3.


## Version 1.11, 3/29/2009:

### Major changes:

- Revised effectplot so that one may refer to "pseudomarkers" by their
  chromosome and position with a construction like "5@32.8", rather
  than having to first call find.pseudomarker to get the name.  Also,
  in the actual plot, we use the form "5@32.8" rather than, for
  example, "c5.loc47".

### Minor changes:

- Revised the software license statements throughout the source code
  (and above), for clarity and consistency.

- Revised write.cross so that it may write doubled haploid (type "dh")
  data.

- Fixed a problem with movemarker regarding the treatment of
  chromsome names.

- Slight change to makeqtl so that QTL names of the form "5@30.0" have
  ending 0's left in, if appropriate (so that if one QTL is referred
  to as "1@12.23", then another like "5@30" will be given as
  "5@30.00", so that all have equal precision.

- Added a function find.pseudomarkerpos for finding the position
  corresponding to a "pseudomarker" name (similar to find.markerpos).

- Added an internal function charround() for rounding numbers, turning
  them into character strings with ending 0's preserved.

- Fixed a slight bug in lodint() regarding the case of multiple
  positions sharing the maximum LOD score.

- In dropfromqtl, addtoqtl, replaceqtl, and reorderqtl, attributes
  "formula" and "pLOD" are now stripped.

- In replaceqtl, changed the argument "indextodrop" to just "index",
  as in dropfromqtl.

- Added a more clear error message in plot.map in the case that a
  sex-specific map and a sex-averaged map are input.

- In plot.info and plot.geno, we now allow one to use main as an
  argument for producing a self-defined title.

- Slight change in plot.qtl, regarding placement of text and size of
  arrows.

- In summary.fitqtl, print.addint, and print.addcovarint, eliminated
  an extraneous blank line after the model formula.

- Added code from Pjotr Prins enabling R/qtl to be linked against Perl
  and Ruby, as part of biolib.


## Version 1.10, 1/11/2009:

### Major changes:

- Revised the way that the 'chr' argument is treated in functions such
  as scanone, scantwo, etc., to give greater consitency.  Numbers are
  interpreted as character strings to be matched to the chromosome
  names.  Negative numbers and character strings that start with "-"
  are interpreted as omitting the corresponding chromosomes, matched
  by name.  One may also use a logical vector (TRUE/FALSE), of the
  same length as there are chromosomes, indicating which chromosomes
  are to be considered.  So if an object has three chromosomes named
  "1", "3", "4", using chr=2 will result in an error, while chr=3 will
  give the second chromosome (named "3").

- Also revised the way that the 'ind' argument is treated in
  subset.cross and plot.geno, in the case that the input cross
  contains individual identifiers in the phenotype data.  The 'ind'
  argument can still be a logical vector, but otherwise we first seek
  to match the values against individual identifiers.  For identifiers
  that are character strings, one may use "-" at the beginning of each
  to indicate all individuals except those given.

- Added an argument 'batchsize' to scanone and scantwo, so that in the
  case that multiple phenotypes (or permutations) are to be run as a
  batch (with method "hk" or "imp"), they can be run in smaller
  batches (indicated by batchsize).  This can speed things up quite a
  bit in the case of a very large number of phenotypes (or
  permutations).

- Added two functions for the de novo construction of a genetic map.
  formLinkageGroups uses pairwise marker linkage information
  (calculated with est.rf) to partition markers into linkage groups.
  orderMarkers uses a quick but not very good algorithm for ordering
  the markers on a chromosome (minimizing the number of obligate
  crossovers).

- Added a function addcovarint, which is similar to addint, but adds
  one QTL x covariate interaction at a time.

### Minor changes:

- Created functions replacemap.scanone and replacemap.scantwo, which
  enable one to plot scanone or scantwo results relative to another
  map (with positions interpolated based on marker locations).  These
  can be used, for example, so that one may plot scanone or scantwo
  results relative to a physical map.

- Added a function replacemap.cross, which is the same as the long
  extant function replace.map.  This was so that I could the functions
  replacemap.scanone and replacemap.scantwo (see above).  All can be
  used with the 'generic' function replacemap.

- Revised the functions nchr, nmar, totmar, so that they work for map
  objects as well as cross objects.  Separated the help files for
  nind, nmar, totmar, nphe, nchr from the summary.cross help file, so
  that the use of the revised functions can be explained.

- Revised the way in which LOD score columns are renamed in c.scanone
  and cbind.scanoneperm. If labels are given, these are appended to
  the end of the names, but if labels are not given and there are no
  repeats in the column names, the column names are left as they were.

- Added functions subset.scanoneperm and subset.scantwoperm for
  pulling out selected LOD columns in the case that permutation tests
  were run with multiple phenotypes.

- Revised calc.penalties so that it can deal with the case of scantwo
  permutations for multiple phenotypes: Added an argument "lodcolumn"
  for selecting the phenotype; if missing, penalties for all phenotype
  are calculated.

- Added a function markernames for pulling the marker names out of a
  cross object (as one long vector).

- summary.cross now checks for duplicate chromosome names and
  chromosome names that start with '-', either of which would cause
  major problems.

- Revised summary.cross so that, if the phenotype data are missing
  column names, an error is given.

- Revised est.map so that the chromosomes are given classes "A" or "X"
  according to the chromosome types in the input cross object.

- Revised geno.crosstab to deal with partially informative genotypes
  in an intercross or 4-way cross.  Added an argument so that (by
  default) columns and rows with no data will not be printed.

- Revised comparegeno to have the option what="both", through which
  the result has the proportion of matches in the lower triangle and
  the number of matches in the upper triangle.  Also, now if
  what="proportion" the diagonal has all missing values (rather than
  all 1's); otherwise the diagonal contains the number of typed
  markers for each individual.

- Revised movemarker so that one can move a marker onto a totally new
  chromosome.

- Revised sim.cross so that if the input map object has no chromosome
  names, the chromosomes in the output cross object still have names.

- Revised summary.cross to give a warning if the chromosomes are not
  named.

- Added a function convert2sa for converting a sex-specific map object
  to a sex-averaged map object by pulling out the female marker
  locations (and issuing a warning if the female and male locations
  are very different).  This is useful for plotting a simpler version
  of a map estimated for a 4-way cross via est.map with sex.sp=FALSE.

- Slight change to countqtlterms(), used by stepwiseqtl(), to skip the
  parsing of interactions in the case that there are 0 or 1
  interactions; this might speed things up slightly.

- Fixed a slight with plotModel; the lines indicating interactions
  were a bit askew.

- Slight changes in ripple and bayesint, changing use of rev(order(.))
  to order(., decreasing=TRUE)

- Added a more clear error message in the case that the number of
  individuals in the cross doesn't match the number of individuals in
  the QTL object in fitqtl, stepwiseqtl, addqtl, addint, addpair,
  refineqtl.

- Fixed a couple of bugs in refineqtl, one concerning convergence and
  the other concerning dropping individuals with missing covariates or
  phenotypes in the case that there is just one covariate.

- Fixed a slight problem in c.cross: if maps are not precisely the
  same, we don't try to combine the genotype probabilities.

- Fixed a bug in summary.scanone with format="allpeaks" for the case
  that there are no peaks meeting the threshold/alpha.

- Fixed a bug in switch.order related to the change in references to
  chromosomes.


## Version 1.09, 7/18/2008:

### Major changes:

- Added a function stepwiseqtl() for performing forward/backward
  selection to identify a multiple QTL model, with model choice made
  via a penalized LOD score, with separate penalties on main effects
  and interactions.

- The documents "Brief tour of R/qtl" and "New functions for exploring
  multiple-QTL models" were revised to discuss the stepwiseqtl
  function.

- calc.penalties() uses permutation results for a 2-dimensional, 2-QTL
  scan to derive penalties for the penalized LOD scores used by
  stepwiseqtl().

- plotModel() is a new function for creating a simple graphical
  representation of a QTL model.

- Added an additional cross type, "dh", for doubled haploids.  This is
  treated like a backcross, though genotypes will be indicated as
  homozygotes.  Changed a whole bunch of functions very slightly to
  accommodate this.

- The argument pheno.col in many functions can now be a vector of
  numeric phenotypes.  This could be useful for studying the results
  with various transformations of a phenotype.  The vector has to be
  numeric, has to have the length equal to the number of individuals
  in the cross, and has to contain either non-integers or values
  outside the range 1,2,3,...n.phe.  (The revised functions are
  scanone, scantwo, addqtl, addint, addpair, cim, effectplot, fitqtl,
  plot.pxg, plot.pheno, refineqtl, scanoneboot, scanqtl, stepwiseqtl.)

- lodint and bayesint were revised to accept qtl objects output by
  refineqtl (with keeplodprofile=TRUE).  An additional argument,
  qtl.index, was added to indicate for which QTL (within such qtl
  objects) the approximate confidence intervals should be derived.
  For scanone output, the functions were modified so that, if the
  results concern just a single chromosome, the chromosome argument is
  not needed.

- In refineqtl, the default for the argument keeplodprofile is now
  TRUE.  The LOD profiles contained in the output of refineqtl are now
  of class "scanone".

### Minor changes:

- Covariates (argument covar) in fitqtl, scanqtl, addqtl, addpair,
  addint, and refineqtl can now be a numeric matrix (with column
  names), and not just a data.frame.

- Permutation results obtained via scantwo now include the maximum LOD
  score from a single-QTL scan.  This is not used in summary.scantwo,
  but is included for completeness.

- Added an argument 'pvalues' to summary.fitqtl and addint; if FALSE,
  the pvalues are not displayed.

- Added a function ntyped() which is just like nmissing() except it
  gives the opposite thing (no. genotypes per individual or marker).

- est.rf, for a backcross, had been replacing estimated rec frac > 0.5
  with 0.5.  This is no longer done.

- The print and summary function for QTL objects now will print the
  formula and penalized LOD ("pLOD") if they exist as attributes.
  They also take into account the case of a null QTL model.

- Revised makeqtl so that QTL names are of the form "1@15.0" rather
  than "Chr1@15.0"; the "Chr" seemed gratuitous.

- Changed some of the examples in the help files to use pull.pheno in
  place of references to cross$pheno.

- Revised summary.cross so that it pays attention to options("width")
  and prints things more nicely if there are loads and loads of
  phenotypes or whatever.

- We now allow formulas in fitqtl, refineqtl, addqtl, etc., to be
  character string representations of formulas.  (They are then
  converted.)

- Added a function cbind.scanone (which is identical to c.scanone).

- The output of fitqtl now includes an element "lod" containing the
  LOD score from the fit of the full model.

- In reorderqtl, if the argument neworder is not provided, the QTL are
  ordered by chromosome and then by position within a chromosome.

- We define a new class, "compactqtl", for defining the trace through
  model space in stepwiseqtl() (if called with the argument
  keeptrace=TRUE).  This is similar to the class "qtl" (of QTL objects
  created by the makeqtl() function), but containing just chromosome
  IDs and positions of QTL.

- Fixed a bug in effectplot in the case that not all possible
  genotypes are observed in the imputations.

- Slight change in plot.scanone so that one may use xlab as an
  argument.  Similarly changed plot.map so that one may use xlab and
  ylab to change the x- and y-axis labels.  Similarly changed
  plot.scantwo so that one may use xlab and/or ylab to change the x-
  and y-axis labels. (In plot.scantwo, if just one of xlab or ylab is
  given, the other is assumed to be the same.)

- Slight change to summary.fitqtl and summary.addint so that very long
  formulas are split across multiple lines.

- Slight change to scanqtl to avoid going just outside the defined
  intervals (specifically, so that refineqtl and plotLodProfile do not
  go just past the flanking QTL).

- Slight change in plot.scantwo to avoid affecting par("mfrow") or
  par("mar") in the case zscale=FALSE.  Also fixed a slight bug
  regarding "any(contours)>0" vs "any(contours>0)".

- Slight change in the way refineqtl determines convergence, to try to
  avoid unnecessary iterations.

- Slight change to plot.scanone so that the chr argument can be logical.

- Simplied code for 'cat' statements in many places; we'd used
  cat(paste(...)) and the call to paste wasn't needed.

- Fixed an error in scanqtl that arose if there are multiple markers
  at the same position.  Now give a warning message.

- The title of the plot produced by plot.rf can now be modified by
  including the argument 'main' in the call.

- In plot.geno, if chr is missing, we plot the genotypes for the first
  chromosome. Changed the label "Position (cM)" to "Location (cM)"
  (just for consistency across functions).

- In geno.crosstab, changed the column and row labels for missing data
  from "NA" to "-".

- Made a slight change to plot.pxg regarding the locations of the
  genotype labels on the x-axis.

- Revised print.map so that it doesn't print the log likelihood
  attribute.

- Fixed a bug in addpair if the formula is missing.

- Fixed a bug in est.rf for 4-way crosses.

- Fixed a bug in c.scanoneperm, c.scanone, and c.scantwoperm regarding
  the "df" (degrees of freedom) attributes.

- Fixed a bug in read.cross with format "csv" (or "csvs" or "csvr" or
  "csvsr") so that dec="," can be used as an argument.

- Fixed a potential bug in read.cross with format "mm" for pulling out
  the cross type.

- Fixed a bug in add.threshold for the case of multiple phenotypes.

- Fixed a bug in scanone permutations in the case of a single
  covariate with some individuals with missing phenotypes and/or
  covariates.

- A more meaningful error message is given in makeqtl in the case that
  multiple markers are at identical positions, so that qtl locations
  cannot be determined.

- Fixed a bug in read.cross with format="csvs" or "csvsr" for the case
  that there are individuals with phenotypes but no genotypes.

- Fixed a bug in refineqtl for the case of linked QTL.  (Incorrect
  limits for the search intervals were used.)

- Fixed a bug in plotLodProfile for the case of linked QTL.  (LOD
  profiles were not placed correctly.)

- Slight change to the internal function reviseqtlnuminformula(),
  so that the input formula can be a character string.

- Fixed slight bug in scanone and scantwo [any(weights)<=0 changed to
  any(weights<=0)].  Fixed a similar bug in plot.geno [any(errors)
  changed to any(errors != 0)].


## Version 1.08-56, 4/8/2008:

- Fixed a bug in scanqtl.


## Version 1.08-55, 4/3/2008:

### Major changes:

- For all functions taking a "pheno.col" argument (including scanone
  and scantwo), this argument can now be a character string indicating
  the name of a phenotype. (Previously, it had to be a numeric index
  indicating the phenotype).

- fitqtl now takes a cross object and pheno.col (as with scanone,
  scantwo and scanqtl), rather than a column of phenotypes.

- Implemented Haley-Knott regression for the fit of multiple-QTL
  models in fitqtl and scanqtl.

- Added functions addint, addqtl, and addpair, for exploration of
  multiple QTL models.  addint tries adding all possible pairwise
  interactions, one at a time, to a multiple QTL model.  addqtl scans
  for an additional QTL to be added to a multiple QTL model.  addpair
  scans for an additional pair of QTL to be added to a multiple QTL
  model.

- Added functions addtoqtl, dropfromqtl, and replaceqtl, for
  manipulating a qtl object created by makeqtl().

- Added a function refineqtl(), for getting the maximum likelihood
  estimates of QTL positons (as best we can) in the context of a
  multiple-QTL model.  Added a function plotLodProfile which can
  create a figure with 1-dimensional LOD profiles for each QTL, in the
  context of a multiple QTL model, as is commonly created for multiple
  interval mapping.

- Added a function, tryallpositions(), for testing all possible
  positions for a given marker, keeping the order of all other markers
  fixed.

- Added a function, compareorder(), for comparing a given order of
  markers on a single chromosome to the current one contained within a
  cross object.

- Added some additional marker genotype codes for the phase-known
  4-way cross, for a dominant marker with both parents being
  heterozygous: 11 = not AC, 12 = not BC, 13 = not AD, 14 = not BD.

- Revised est.rf, for estimating recombination fractions between all
  pairs of markers, so that it can give results for many of the
  incompletely informative markers in a 4-way cross.

### Minor changes:

- Added an additional argument (expandtomarkers) to lodint, bayesint,
  and summary.scanoneboot.  If TRUE, the intervals provided are
  expanded to the nearest flanking markers.

- Added a function geno.crosstab for creating a cross-tabulation of
  the genotypes at two markers.

- Added a function pull.pheno for pulling out the data for a phenotype
  or phenotypes.

- Added a function countXO for counting the number of obligate
  crossovers for each individual across the genome or on individual
  chromosomes.

- Added a function plot.qtl, for plotting the locations of QTL in a
  qtl object against the genetic map.

- Added a function checkformula for checking the formula in
  fitqtl/scanqtl, to ensure that it satisfies the hierarchical
  structure we assume: if a term is involved in an interaction, its
  main effect should also be included.

- Added an additional argument to fitqtl, run.checks.  If TRUE, we
  check the input formula and look for individuals with missing
  phenotype or covariates.  This is included so that the checks are
  not repeated multiple times when scanqtl calls fitqtl.

- Added the ability to calculate joint QTL probabilities assuming
  conditional independence of QTL genotypes given markers genotypes.
  (An approximation, but it speeds up scantwo slightly, for a chr
  versus itself.) Added an argument assumeCondIndep to scantwo().

- Added an argument "zmax" to the plot.rf function, for controlling
  the color scale of LOD scores.  Values at zmax are red; values above
  zmax are thresholded at zmax.

- Added functions plot.scanoneperm and plot.scantwoperm for plotting
  histograms of the permutation results from scanone and scantwo.

- Added a function plot.scanoneboot, for plotting a histogram of the
  results of scanoneboot.

- scanoneboot now stops with an error if the argument pheno.col
  indicates multiple phenotypes.

- Revised find.marker to have an argument "index" which may be used in
  place of the "pos", to find marker names by their numeric order
  within a chromosome rather than by map position.

- In read.cross, with formats "csvs" or "csvsr", we now allow that
  some individuals have phenotypes but no genotypes and vice versa,
  and the individuals in the genotype and phenotype files are not
  required to be in the same order.

- In the getsex() function, for pulling out sex and pgm for all
  individuals, we now attempt to infer the status of individuals with
  missing information.  Warnings are given.

- Fixed a bug in read.cross.qtx, concerning the case that genotypes
  are like H:B or A:B, and need to be converted to A:H.

- Fixed a bug in effectscan for the X chromosome in the case of an
  intercross with both directions but just one sex.

- Slight change to plot.geno() to make individual IDs shown rather
  than just numbers, if they are available.

- Slight change to effectscan, to pass the "..." argument to the plot
  function, so that, for example, you can use the 'main' argument to
  put a title on the plot.

- Slight changes to top.errorlod() and getid() to deal with a bug for
  the case that there is an "ID" phenotype column with names like
  "1_F1".

- Fixed some code in scanqtl regarding dropping individuals with
  missing phenotypes and/or covariates that really slowed things
  down.  Revised the analogous code in fitqtl.

- Added an argument to est.map: omit.uninformative.  If TRUE (which
  is the default, and which was previously the only option),
  individuals with fewer than two typed markers are omitted.  This was
  added for use by the new function tryallpositions().

- Added a function markerloglik, for calculating the log likelihood
  for a fixed marker.  This was added for the use of the new function
  tryallpositions().

- read.cross now prints the cross type at the end.

- read.cross (with format="mm" or one of the "csv" formats) gives a
  more explicit warning message if phenotypes are to be treated as
  missing.

- Revised summary.qtl to also indicate the number of imputations, in
  the case that the qtl object contains them.

- cim() had previously used a single column for each covariate in an
  intercross (that is, it assumed additivity of alleles); this is
  fixed: it now uses two columns.  A revised forward selection
  algorithm for intercrosses was written, to select these pairs of
  columns together.

- Fixed a slight bug in lodint and bayesint that changed the marker
  names if the LOD peak was at one end of the interval.

- subset.scantwo will now drop X chromosome related stuff from the df
  attribute if the X chromosome has been omitted.

- Changed the default tolerance in est.rf and est.map to 10^-6 (rather
  than 10^-4).  Increased the default maxit (maximum no. iterations)
  to 10000.

- Added an extra column in the results of summary.map, giving the
  maximum distance between markers on each chromosome and overall.

- QTL objects produced by makeqtl now include a component "altname",
  which will be like "Q1", "Q2", ...   Changed fitqtl to look at this
  rather than at the column names of qtl$geno.

- Made a slight change to plot.scantwo regarding z-limits when zlim
  is missing and allow.neg=TRUE.

- The objects of class "scanoneperm" and "scantwoperm" now have a
  secondary class (either "matrix" or "list")

- Fixed a bug in scantwo perms for the use of the argument clean.scantwo.

- Fixed a bug in summary.cross in the case of invalid genotypes in a
  cross.

- Fixed a slight bug in subset.cross regarding the recombination
  fractions from est.rf().

- Fixed slight bugs in effectplot and reviseXdata regarding the X
  chromosome in an intercross with both sexes and one cross direction.

- Fixed a bug in CIM for the case of multiple marker covariates on the
  same chromosome.

- Fixed a bug in revisecovar() regarding dropping of covariates for the X
  chromosome.

- Fixed a bug in add.threshold.

- Fixed a bug in movemarker for the case of sex-specific maps and
  with the marker being moved to the middle of a chromosome with
  exactly two markers.

- Fixed a bug in geno.table; missing genotypes weren't shown for the X
  chromosome.

- Fixed a bug in fit.stahl.

- Subtle modification to a few of the help pages to conform to a
  change in R.

- Changed the default for plot.geno back to a horizontal plot
  (horizontal=TRUE); changing it to vertical was a bad idea.  Fixed a
  slight bug regarding F2-type markers in 4-way cross.

- In write.cross with the "qtlcart" format, if there is a previous
  file that would otherwise be overwritten, it is now moved to a file
  with extension ".bak" rather than ".mov".  Also made some slight
  revisions to get things in the map file to line up.

- Made a minor change to plot.map, so that if the "..." contain xlim
  or ylim, they are used in place of the defaults.


## Version 1.07, 9/20/2007:

### Major changes:

- Completely rewrote the effectscan function, so that it now uses
  multiple imputation results and deals with the X chromosome
  appropriately.

- Fixed an important bug in fitqtl, in which incorrect results could
  be obtained if covariates were placed before QTL terms in the
  formula.

### Minor changes:

- Added an argument "alternate.chrid" to plot.scanone, plot.scantwo,
  plot.info, plot.missing, geno.image, plot.map, plot.cross,
  effectscan, plot.errorlod, and plot.rf.  If TRUE, the placement of
  chromosome ID axis labels is alternated, so that they may be more
  easily distinguished.  For plot.cross, alternate.chrid=TRUE has been
  made the default; for the other functions, FALSE is the default.

- In the output from makeqtl, "pos" is now the precise position of the
  pseudomarkers (rather than just the input values), and the QTL names
  reflect that (though they are rounded).  replaceqtl and addqtl were
  similarly revised.

- In makeqtl, added an argument what=c("draws","prob"); we now pull
  out either the results of sim.geno or the results of calc.genoprob,
  and not both.  (Only the former is needed at this point; the latter
  will be used once we have implemented EM/HK/eHK in fitqtl.

- plot.geno now works for a 4-way cross; we changed the default to be
  a vertical plot (ie, horizontal=FALSE).

- Added functions print.qtl, summary.qtl, print.summary.qtl, for
  getting simple information about a QTL object.

- Fixed slight bugs in pull.map and replace.map.

- Fixed a slight bug in top.errorlod, for the case that there are IDs
  (e.g., in cross$pheno$id) that are not numeric.

- Fixed a slight bug in scanqtl.

- Changed checks regarding the class of the input to various functions
  to be a bit more permissive.

- Fixed a potential problem (which shouldn't be realized) in reviseXdata.


## Version 1.06, 8/7/2007:

### Major changes:

- Revised the method for calculating genotyping error LOD scores.  For
  each individual and each marker, the error LOD score is calculated
  assuming that all other genotypes for that individual on that
  chromosome are correct.  The new procedure requires much more
  computation time (especially in the case of dense markers), but
  identifies many additional potential errors.  A new argument,
  version, allows one to specify use of the "new" or "old" version of
  the error LOD score calculations.

- Added a function geno.image for plotting an image of the genotype
  data.  This is much like plot.missing, but gives the genotypes in
  color, rather than just black/white indicating missing/not.

- Revised geno.table so that it gives reasonable p-values for the X
  chromosome and for the case of dominant markers in an intercross.
  Added a chr argument to obtain results for only selected
  chromosomes.

- Added a function cim() for performing composite interval mapping by
  one of the schemes used in QTL Cartographer: forward selection at
  the markers, to a fixed number of markers, followed by interval
  mapping using those marker as covariates, and dropping any markers
  within some fixed window around the position under test.  The
  results may be plotted or summarized using the functions for output
  from scanone().  Also added a function add.cim.covar, for adding
  dots, to a plot from plot.scanone(), to indicate the selected marker
  covariates.

- Extended the code for fitting the Stahl model for crossover
  interference to the case of intercross data.  (Modified the
  functions est.map and fitstahl, and the underlying C code.)

- Added functions scanoneboot and summary.scanoneboot, for deriving
  bootstrap confidence intervals for the location of a QTL, but we
  recommend using lodint or bayesint, instead.

### Minor changes:

- Added a function find.markerpos(), for finding the chromosome and
  position of a marker (or vector of markers).

- Added arguments 'xlab', 'ylab', and 'col' to effectplot(), so that
  you can override the defaults.

- Added a function add.threshold() for adding a significance threshold
  (estimated via permutation results) to a plot created by
  plot.scanone().

- Revised plot.pxg() so that unobserved genotypes will not be
  displayed.  This was needed for the plot of two-locus genotypes on
  the X chromosome.

- Changed the names in two of the columns in the output from the 2d
  permutation test, to be "fv1" and "av1" rather than "2v1.int" and
  "2v1.add", to correspond more closely to the names in the
  summary.scantwo output.  Also, we now allow the argument "alphas" to
  be a single number, in which case it is assumed that the same
  significance level is to be used for all five LOD scores.

- Made a slight change in summary.scanone(), so that when p-values are
  provided, the rownames don't get lost.

- Revised the hyper data set slightly; changed the "alleles"
  attribute, to be c("B","A") rather than c("A","B"), as this was a
  backcross to the B strain.

- In scanqtl, if no a fixed model (with no scanning) is fitted, the
  output is now just the LOD score for that fitted model.

- clean.cross was dropping any "alleles" attribute; this is now fixed.

- Added a "lodcolumn" argument to lodint() and bayesint().

- Fixed a bug in scanqtl; in two-dimensional scans, the first row of
  the results was wrong.

- Fixed a bug in c.cross concerning combining backcross and
  intercrosses.

- Fixed a bug in effectplot() regarding pseudomarkers with names like
  "c3.loc42.5".

- Fixed a bug in discan, for the case of method="mr".

- Fixed a bug in write.qtlcart regarding RIL.

- Fixed a bug in max.scanone for the case that there is more than one
  locus with the maximum LOD score.  (In that case, we print a random
  locus, among those having the maximum LOD.)

- Revised summary.scanone so that the rule is to pick out LOD scores >
  (rather than >=) the threshold.

- Revised -.scanone, -.scanoneperm so that very small differences get
  set to 0.

- Added ability to halt calculations via Ctrl-c in many of the C
  routines. (Previously, you'd have to wait for an exit from the C
  code.)

- There was a bug in sim.cross() regarding the QTL effects for a
  backcross.

- Fixed a bug in drop.markers, for a 4-way cross.

- Fixed a bug in est.map, for a 4-way cross and the case of
  sex-averaged maps.

- Made a slight change regarding estimating rec fracs in 4-way cross
  (possibly immaterial).

- A slight change in the C code for imf_stahl.

- Modified locate.xo() slightly...when there is no crossover, it
  gives numeric(0) rather than NULL.

- Fixed a slight problem with the format of the degrees of freedom in
  scanone permutation results.

- Slight change in the color scheme in effectplot()

- Fixed a slight bug in fitstahl() that made it crash if none of
  m, p, and error.prob were specified.  Also revised the function so
  that we look only for error.prob <= 0.5.

- Fixed a slight bug regarding sex-specific maps in the create.map()
  function.

- Slight revision in plot.geno, so that if the "ind" argument contains
  duplicates of individuals, only the unique individuals are
  plotted.

- Slight revisions to scanone, scantwo, and discan, so that, in the
  midst of permutations, just one instance of various warnings is
  printed.

- Changed the names of some of the sample data files.

- Revised the fitstahl function to use of the estimated map for one
  value of m as the starting point for the next value.  This can
  really cut down on the required EM iterations and so speeds things
  up.

- In the "map" component of the output from scantwo, the name of the
  second column is now "pos" rather than "map", to correspond more
  closely to the scanone output (and because it is more appropriate).

- Fixed a bug in write.cross; in an intercross with all males and all
  pgm==1, all X chromosome genotypes got converted to AA.

- Added a few more verbose error messages in read.cross for the "csvs"
  format (contributed by Steffen Moller, University of Lubeck).

- Fixed a slight bug in scanone with model="2part" or model="binary",
  that showed up if one first used jittermap().

- Added arguments maxit, tol and sex.sp to switch.order()

- Fixed a bug in movemarker() for the case of a 4-way cross.

- In write.cross.gary(), changed a couple of uses of 'T' and 'F' to
  'TRUE' and 'FALSE', respectively.

- Fixed find.markerpos so that it works with a 4-way cross.

- Revised plot.scantwo so that the upper and lower arguments can take
  values "fv1" and "av1" as aliases for "cond-int" and "cond-add",
  respectively.


## Version 1.05, 11/8/2006:

### Major changes:

- Fixed a problem with permutations in scanone and scantwo: covariates
  weren't being permuted to match the phenotypes.

### Minor changes:

- Fixed a slight bug in scanone for model="binary".

- Revised write.cross for the "qtlcart" format; there was a problem in
  the case that there were many markers on a chromosome.

- Revised the C code for imf_stahl, to better deal with void pointers.


## Version 1.04, 10/28/2006:

### Major changes:

- R/qtl has a new web site: rqtl.org

- Revised the format for the output from scantwo.  Added a function
  convert.scantwo for converting from the previous format to the new
  format.  For scantwo results calculated with R/qtl version 1.03 and
  earlier, you'll need to use convert.scantwo to convert them to the
  new format in order to use the summary.scantwo and plot.scantwo
  functions.  (In the previous format, joint and epistasis LOD scores
  were stored; now we store the joint LOD and the LOD from the
  additive QTL model.  This is so that, if there is a problem with the
  joint model, it won't corrupt the results for the additive model.)

- Eliminated the 'run.scanone' argument from the scantwo() function.
  scanone is always run.  The summaries and permutation tests require
  these results.

- plot.scantwo now has an 'upper' as well as a 'lower' argument, for
  complete control over what gets plotted.

- Completely revised the summary.scanone and summary.scantwo
  functions.  I have written documents to explain the use of the new
  functions.  These are distributed with the code and are also
  available at the R/qtl website (https://rqtl.org), under
  "Tutorials".

- summary.scanone:  There is now a format argument, useful for the
  case that the scanone result contains multiple LOD score columns
  (for example, for multiple phenotypes).  We may focus on a single
  LOD column (format="onepheno"), as was done before; include
  different rows for the peaks in each LOD column (format="allpheno");
  or have one row per chrosome, containing the the position and LOD
  score for each the peak from each LOD column (format="allpeaks").
  The function now also can take permutation results in order to
  automatically calculate LOD thresholds or to calculate
  genome-scan-adjusted p-values.

- summary.scantwo: This was quite radically changed.  For each pair of
  chromosomes (including a chromosome with itself), we calculate five
  LOD scores: the maximum LOD for the full model (2 QTL +
  interaction), the maximum LOD for the additive model, the difference
  between these (which concerns a test of whether the two loci
  interact), and two LOD scores concerning 2  vs 1 QTL: the difference
  between the full LOD and the best single-QTL LOD for the pair of
  chromosomes, and the difference between the additive LOD and the
  best single-QTL LOD for the pair of chromosomes.  This is the
  recommended output, indicated via the argument what="best".  One may
  also set the 'what' argument to "full", "add", or "int".  (See the
  help file for summary.scantwo.)  The 'thresholds' argument now
  requires five values, or one may provide permutation results plus a
  set of five 'alphas' (significance levels).  There is also an
  argument 'allpairs'; the default is TRUE, in which case all pairs of
  chromosomes are considered.  If allpairs=FALSE, only the self-self
  chromosomes are considered, so that one may look more easily for
  cases of possible linked QTL.

- summary.scanone and max.scanone now can just just one object, rather
  than multiple such, as before.  However, we have added functions
  c.scanone and cbind.scanoneperm for combining the columns in
  multiple runs of scanone (generally either multiple phenotypes or
  multiple methods).

- Completely revised the summary.scantwo function.  The old version is
  saved as the function summary.scantwo.old().

- Added an argument perm.strata to scanone and scantwo, to allow
  stratified permutation tests.  If provided, it should be a vector
  of length the number of individuals in the cross; unique values in
  perm.strata will specify the strata in which the permutations should
  be performed.  (For example, this could be an indicator of the sexes
  of the individuals, in which case the individuals will be shuffled
  separately within males and within females.)

- Permutations in scanone can now be done to give separate thresholds
  for the autosomes and the X chromosome.  The argument perm.Xsp is
  used to indicate that this should be done, in which case many more
  permutations will be run for the X chromosome than for the
  autosomes, to ensure similar accuracy.  The output of scanone when
  n.perm>0 is now given class "scanoneperm", and we've written a
  function summary.scanoneperm for getting LOD thresholds.  (This is
  necessary, since the calculation autosome- and X-chromosome-specific
  thresholds is a bit complicated.)  We've also added a function
  c.scanoneperm for combining the results of multiple permutation
  runs.  (This because their combination is not so simple as before.)

- Added functions -.scanoneperm, +.scanoneperm, -.scantwoperm, and
  +.scantwoperm, for taking sums or differences of permutation results
  from scanone or scantwo.  This is particularly useful for getting
  LOD thresholds for QTL x covariate interactions, though one must be
  careful to ensure that the permutations are perfectly linked, which
  can be achieved with set.seed.

- The permutation results from scantwo are completely changed.  Rather
  than keep track of the maximum LOD score for the full model (two
  QTLs + interaction) and the interaction LOD score, we keep track the
  genome-wide maxima of the 5 LOD scores calculated in
  summary.scantwo. (See above.)  The output is now given a class
  scantwoperm, and there is a summary.scantwoperm function for
  calculating LOD thresholds.  There is also a c.scantwoperm function
  for combining results from multiple runs, largely for the case that
  multiple sets of permutations were run in parallel.

- Removed the ability, in scanone and scantwo, to use the snow package
  to do parallel analysis on a linux cluster.  With the new changes in
  scanone, I found the code to be too cumbersome.

- Added the extended Haley-Knott method (see Feenstra et al., Genetics
  173:2269-2282, 2006) to the scanone function.  This is faster and
  more robust than standard interval mapping, and is a better
  approximation (but slower) than regular Haley-Knott regression.

- Revised the sim.cross() function so that, with a backcross, the
  effect in the "model" argument is to be specified as the difference
  between the average phenotypes for the heterozygotes and the
  homozygotes.  (Previously, it was 1/2 this, which is
  different from the typical parameterization.)

- Added the ability, for a backcross, to estimate a genetic map under
  the Stahl model for crossover interference (of which the chi-square
  model is a special case).  Also added a function, fitstahl, for
  getting the maximum likelihood estimates in the Stahl model (or the
  chi-square model); the genotyping error probability may be treated
  as known or may also be estimated.

- Added an argument, "use", to scanone and scantwo, for indicating, in
  the case that multiple phenotypes are to be run, whether only
  individuals with complete data on all phenotypes
  (use="complete.obs") or all individuals (use="all.obs") are to be
  used.

- The degrees of freedom are added as attributes in the output from
  scanone and scantwo, including the case of permutations.

## Minor changes:

- In read.cross, the symbol "#" is no longer treated as a comment
  character by default.  The default is to use comment.char=""; that
  is, no symbol is treated as an indicator of comments.  For the
  comma-delimited file formats, one may have a character interpreted
  as indicating comments using the comment.char argument.

- Fixed a bug for the new "batch mode" permutations in scanone and
  scantwo with method="hk" or ="imp".  The trick only works if there
  is no X chromosome or all individuals have the same sex and cross
  direction or permutations are done stratified within sex and
  direction.

- Added an argument "show.marker.names" to plot.map and plot.scanone,
  so that marker names can be added to these plots.

- write.cross can now write data in the "csvr", "csvs" and "csvsr"
  formats.

- Added a function plot.pheno for plotting a histogram or barplot of
  a phenotype distribution.

- Added a function condense.scantwo for producing condensed versions
  of scantwo output, containing just the maximum LOD scores on each
  pair of chromosomes.  One can get summaries from these but not
  plots.

- In plot.info, added step, off.end, error.prob and map.function
  arguments.  The function now always calls calc.genoprob rather than
  relying on the values in the data.

- Changed the default cutoff for genotyping error LOD scores in
  top.errorlod and plot.geno to 4.

- Changed the name of the function clean() to clean.cross() and made
  a new function clean() that will dispactch a cross object to
  clean.cross.  Added a function clean.scantwo() for cleaning up the
  output of scantwo: any values that are missing or are < 0 are
  replaced by 0 and any LOD scores for pairs of loci that do not have
  a marker between them are set to 0.

- The output of scantwo() now contains the original genetic map as an
  attribute; this is needed for clean.scantwo().

- Slightly revised effectplot() so that if sim.geno hasn't been run,
  it is run (with a single imputation) before the plot is created.
  I also changed the names in the output, so that it is "Means" and
  "SEs".  Changed the examples in the help file for effectplot.

- Added an argument "lodcolumn" to max.scanone, so that it behaves
  like summary.scanone.

- Added a function subset.scanone() for pulling out particular
  chromosomes or LOD score columns from scanone output.

- Added a function find.pseudomarker() for identifying the name of a
  pseudomarker that is closest to a specified position.  This is
  useful for the effectplot() function.

- In scanone and scantwo, there is now a warning printed when
  individuals with missing phenotype are dropped.  Also, a slightly
  better error message is printed if addcovar or intcovar are not
  numeric.  Also added a warning message if addcovar or intcovar
  appear to be over-specified (having columns that need to be
  dropped).

- Fixed a slight problem in the column names from scanone() in the
  case of multiple phenotypes; changed the convention for this.  Now
  the columns will just have the phenotype names.

- Fixed the help file for read.cross concerning the individual
  identifiers and about which files are used for the csvs and csvsr
  formats.

- Revised nmissing() so that if what="ind" and individual IDs are
  included as a cross phenotype (named "id" or "ID"), these are used
  as names in the output.

- Revised summary.cross() to include a check of whether the individual
  IDs (in a phenotype named "id" or "ID") are unique.  If they are
  not, a warning is issued.

- Changed the "pheno" argument in plot.cross to "pheno.col", to be
  more consistent with other functions.

- Changed the name of the argument "which" in plot.rf and nmissing to
  "what".

- Added an argument "verbose" to the ripple function; if
  verbose=FALSE, the function doesn't print anything.  Also modified
  print.summary.ripple so that it prints no more than 6 rows.

- Revised plot.missing so that if reorder=TRUE, the reordering of
  individuals is done by the average of only the numeric phenotypes
  (rather than all phenotypes, which gave an error).

- Fixed a slight bug in summary.cross() regarding duplicate marker names.

- Fixed a bug in scanone() that messed up the X chromosome label in
  the case of method="imp" with multiple phenotypes.

- Fixed a bug in scanone() regarding method="mr-imp" with
  permutations.

- Revised the listeria data set slightly; added an "alleles"
  attribute, with alleles "C" and "B", as this was a cross between
  BALB/cByJ and C57BL/6ByJ, and those strains were coded C and B in
  the original paper.  This takes advantage of a feature added in
  version 1.03.  Also added a phenotype "sex" that indicates all
  individuals are female.

- Added a phenotype "sex" in the hyper data, indicating that all
  individuals are male.

- Modified checkAlleles() so that it doesn't give an error if one
  inputs data for only the X chromosome.

- Revised calc.errorlod() so that it won't give a warning if it has to
  run calc.genoprob() because such probabilities aren't available.  It
  will still give a warning if it has to re-run calc.genoprob() with a
  new error.prob value.

- Slight revision of plot.scanone() to include "Chromosome" as an
  x-axis label if multiple chromosomes are plotted.

- Fixed a slight bug in plot.info for the case that results are only
  at the marker positions.

- Changed the name of the "cols" argument in plot.pxg to "col".  (This
  argument was added in version 1.03.)

- Revised summary.cross so that if the "jittermap" warning is printed,
  it's only printed once.

- Made a slight change to summary.map and print.summary.map, so that
  "sexsp" is an attribute.

- A couple of changes were made to the write.cross.qtlcart():
  round the map locations, make sure backcross code is correct, and
  make sure RIL data is written as 0/2 rather than 0/1.

- Revised plot.info() so that if method="entropy" or
  method="variance", only the column requested is returned.
  (Previously, the other was also given, but with all 0's.)

- Fixed a bug in sim.cross for type="4way"; it would stop with an
  error if QTLs were to be simulated.

- Fixed a bug in fill.geno for the case that a chromosome has just one
  marker.

- Removed the function convert.cross(), which converted cross data
  from the format used in R/qtl version < 0.65 to the current format.
  This shouldn't be needed anymore.

- Fixed a bug in plot.pxg; the plot was messed up one requested an
  autosomal marker and an X chromosome marker, but with the X
  chromosome marker listed first.


## Version 1.03, 7/20/2006:

### Major changes:

- Fixed subset.cross() so that the attribitutes in the results of
  calc.genoprob, calc.errorlod, sim.geno, and argmax.geno, don't get
  lost on subsetting.  This was important to ensure that the package
  will conform to a change that will occur in the next release of R.

- Added a function checkAlleles() for identifying loci that might have
  their alleles switched (in an intercross, if AA and BB are switched,
  or in a backcross if AA and AB are switched).  The X chromosome is
  ignored.  An internal function, checkrf() was removed; this was
  previously called by est.rf() but wasn't well written.

- Added an argument, alleles, for read.cross(), which takes two
  single-character allele labels that are included as an attribute in
  the cross and will be used as labels throughout the program (for
  example, in geno.table, effectplot and plot.pxg).  This required
  numerous small changes throughout the package.

- The maps used for the results from argmax.geno(), calc.genoprob()
  and sim.geno() are now saved as an attribute on the data they
  create (within the cross object) so that create.map() doesn't have
  to be called repeatedly.

- Moved the code for simulating genotype data in sim.cross() into C,
  to increase speed, and modified it so that one may simulate under
  Frank Stahl's interference model (which includes the chi-square
  model as a special case).

- If one includes a phenotype named "id" or "ID", this will be used in
  top.errorlod(), plot.errorlod(), and plot.geno() as identifiers for
  the individual.

### Minor changes:

- Added a warning in summary.cross() if there are multiple X
  chromosomes; the summary now includes the names of the autosomes
  and X chromosome (for diagnostic purposes).

- Revised plot.rf() so that, if the results of est.rf aren't
  available, that function is run.

- Made a slight modification to the bayesint() function for getting
  Bayes credible intervals from scanone() results.

- Added a "chr" argument to pull.geno() and pull.map(), so that you
  can pull out the genotype data or map for a selected set of
  chromosomes.

- Fixed a slight bug in locatemarker(), used by makeqtl().

- Added a function print.map() so that when you print a map, all of
  the class stuff doesn't get in the way.

- Included an additional argument [stepwidth = c("fixed", "variable")]
  in the internal create.map function, for Brian Yandell and the
  R/bmqtl package.  This argument was also added to calc.genoprob(),
  sim.geno() and argmax.geno().  We also added a "stepwidth"
  attribute to the bits that these functions produce.

- Added an argument, "lodcolumn", to the function max.scantwo, for
  picking out a single phenotype in the case that the results concern
  multiple phenotypes.  As with summary.scantwo() and plot.scantwo(),
  this function will only give the results for a single phenotype.

- Made a slight change regarding the sizes of the labels in
  plot.pxg() and added an argument concerning colors.

- Fixed a bug in scanqtl() regarding the X chromosome.

- Revised plot.map() slightly so that one can use "main" as an
  argument, to create a custom title (such as "").

- Fixed a bug in sim.map() regarding sex.sp=TRUE

## Version 1.02, 6/2/2006:

### Major changes:

- Modified scanone() and scantwo() so that they may analyze multiple
  phenotypes simultaneously.  This can greatly speed up the analyses
  with Haley-Knott regression (method="hk") and imputation
  (method="imp").  We can use this trick to speed up permutation
  tests: create multiple permuted phenotypes and then analyze them all
  at once.

- scanone() results no longer include parameter estimates (such as the
  phenotypic averages for each genotype group and the residual SD).
  The bookkeeping for this became too painful as we moved to allow
  multiple phenotypes to be analyzed simultaneously.  To get such
  estimates, use fitqtl().  fitqtl() currently works only via the
  imputation method, but we expect to soon implement multiple interval
  mapping (MIM) and also Haley-Knott regression.

- scantwo() with method="em" now works for the X chromosome, including
  with model="binary".

- Modified the "map10" dataset: a genetic map modeled after the mouse,
  with markers having an approximately 10 cM spacing.  We've revised
  the chromosome lengths to match those in the Mouse Genome Database.

- Dropped support for intercrosses with sex-specific maps (class
  "f2ss").

- Changed the meaning of the argument "lodcolumn" in scanone(). It now
  should be an index starting at 1 rather than 3.  I think this will
  be more clear for the case of LOD scores from multiple phenotypes,
  though perhaps it will be less clear.

- summary.scanone() now takes an argument "lodcolumn"; if a single
  scanone output is given as input, it picks off peaks using those LOD
  scores.  (As with plot.scanone, this is indexed starting at 1.)

- Similarly, added a "lodcolumn" argument to plot.scantwo() and
  summary.scantwo(), for the case that the scantwo results are for
  multiple phenotypes; only results for a single phenotype are used by
  these functions, and "lodcolumn" indicates which one.

## Minor changes:

- Modified summary.cross() to give a warning if there are markers at
  precisely the same position.  Added a function jittermap() to assist
  in fixing this.  Numerous functions run into problems if there are
  markers on top of each other.  Revised the hyper dataset so that it
  doesn't have this problem.

- Switched the order of the arguments "model" and "method" in scantwo()
  to match that for scanone().

- Fixed a slight bug in max.scanone() that led to an unnecessary
  warning message.

- Fixed a slight bug in locate.xo(), used by plot.geno() to identify
  the locations of crossovers.

- Made a slight revision to read.map.qtlcart(), so that chromosome
  names are not required in the map file, and so that more informative
  messages are displayed if marker or chromosome names are not found.

- Modified qtlversion(), which prints the installed version of the
  package, to use library(help=qtl) rather than installed.packages().
  This is a lot faster.

- Added, to summary.cross(), a check on whether multiple phenotypes
  have the same name.

- Modified plot.cross() so that it will work if multiple phenotypes
  have the same name (though that shouldn't happen).

- Fixed a bug in fitqtl() regarding the names of coefficients for
  interactions between QTLs and covariates when get.ests=TRUE.

- Fixed a slight problem in read.cross.csv() and read.cross.csvs() so
  they will read in 4-way cross data.  You need to use genotypes=NULL.

- Modified c.cross() so that it can combined crosses typed at
  different numbers of markers.  The number of chromosomes must be the
  same, and the genetic maps must be consistent.

- Modified plot.scanone() so that you can give it a "ylab" argument to
  override the y-axis labels.  Eliminated the "main" argument, as that
  can be passed via "...".

- Modified the threshold for est.rf() to print warnings about possibly
  switched genotype data.

- Fixed a slight typo in the help file for plot.scantwo().

- Fixed a bug in calc.pairprob() (used by scantwo) for RILs.

- Fixed a slight bug in scanone() for permutations; it dropped
  covariates in permutation tests with model="binary".

- Fixed a slight bug in read.cross with format "csvs", for the
  case that the argument genotypes=NULL.

- Fixed a slight bug in checkcovar() used by scanone() that arose when
  there were missing values in the phenotype data.

- Added a utility function chrlen() for pulling out the lengths of all
  of the chromosomes.


## Version 1.01, 10/25/2005:

### Major changes:

- Revised read.cross to include three additional data formats:

  ```
  "csvr"   The format "csv", but with rows and columns
           interchanged. I call that a "rotated" version of the CSV
           format, but it's really a transposed version.

  "csvs"   The format "csv", but with separate files for the
           genotype and phenotype data.  Note that the first column
           in the phenotype data should specify the individuals'
           IDs, and that there should be a column in the phenotype
           data with precisely the same name, and the individuals
           should be in precisely the same order.

   "csvsr"  The format "csvs", but with the genotype and phenotype
            files rotated (really transposed).
   ```

-  Added example files (of the listeria data) in these formats in the
   "sampledata" directory.

### Minor changes:

- plot.scantwo can now plot the additive LOD scores in the lower
  triangle (with the argument lower="add")

- Fixed a bug regarding the X chromosome in scanone() with the use of
  covariates.

- Modifying plot.geno() to put X's at inferred crossover locations.

- Modified plot.rf() so that the lines between chromosomes are white.

- Added a "chr" argument to plot.map, so that a selected set of
  chromosomes may be plotted.

- Changed the default value for the error.prob argument from 0 to
  0.0001 in est.map(), calc.genoprob(), sim.geno(), and other
  functions.


## Version 1.00, 9/10/2005:

### Major changes:

- Revised fitqtl() so that it can provide estimated QTL effects (using
  the imputation method), though this is working completely only for
  autosomes in backcrosses and intercrosses at this point.

- Revised fitqtl() so that it treats the X chromosome appropriately.

### Minor changes:

- Fixed a bug in read.cross for the "qtlcart" format.  There was a
  problem in the reading of map files if the inter-marker distances
  were at all out of alignment (which occurred when two markers were
  100 cM apart).

- est.rf() was revised to treat the X chromosome properly in an
  intercross.

- Fixed a bug in scanqtl() in the case that covar=NULL was specified;
  this should be treated just like if it were missing.

- Modified plot.scantwo() so that the default is lower="joint"; I've
  become suspicious of lower="cond-int" and "cond-add".

- Added another argument to plot.scantwo(): point.at.max, for plotting
  an X at the maximum LOD.

- Added an argument, "verbose", to scanqtl(), to give feedback about
  progress.

- Revised scanqtl() so that makeqtl() doesn't get called repeatedly,
  but rather the imputations get copied over.  This sped the thing up
  immensely.

- Replaced calls to print.matrix() in print.summary.ripple(), as the
  function will be dropped from R ver 2.2.

- Fixed a bug in plot.pxg(): it halted with an error in the case of
  multiple markers on the same chromosome.

- Fixed a bug in fitqtl() that made it die with only one QTL.

- plot.scantwo() now gives cM locations on the axes in the case of
  just one chromosome.

- Fixed a bug in find.marker.

- Fixed a slight bug in plot.scanone()

- Added a function qtlversion() to print the version number of the
  currently installed version of the package.

- Fixed a bug in plot.scantwo() regarding the X chromosome in the case
  that markers were included in scan but are not to be plotted.

- Fixed a bug in scanone() and scantwo() regarding the use of sex
  and/or pgm as covariates on the X chromosome; this affected only
  results with method="imp" and only for the X chromosome.

- Revised some of the help files so that they conform to the rules for
  the latest version of R.

- Revised scantwo() by imputation so that the interaction LOD score is
  obtained by combining across imputations for each of the full and
  additive model and then subtracting, rather than subtracting and
  then combining.

- Fixed some typos in the help files.

## Version 0.99, 4/26/2005:

### Major changes:

- scanone() now allows additive and interactive covariates in the
  case model="binary" (that is, for a binary trait).  This uses
  a logit link.

- scantwo() now allows analysis of binary traits (model="binary"),
  for method="em" only.

- Added a few new utility functions: bayesint() for calculating
  Bayesian probability intervals [cf lodint()], comparegeno for
  comparing individuals' genotypes, and strip.partials()
  for removing partially informative genotypes.

- Added +.scanone() and -.scanone() for adding and subtracting the
  output of scanone()

- Revised summary.scanone() and print.summary.scanone() so that it can
  summarize the result of multiple scanone() results together.

- Added code (the function sim.cc) for simulating the "Collaborative
  Cross" (8-way RILs) and for calculating QTL genotype probabilities
  and identifying the most likely genotypes (using Viterbi) with SNP
  data on such lines.

### Minor changes:

- Revised plot.scanone() so that when one chromosome is plotted, it
  starts whereever the first marker on the chromosome was placed, and
  not necessarily at 0.

- Revised plot.map() so that shifting chromosomes so that the first
  marker is at 0 cM is optional.  (Use shift=FALSE to not do such a
  shift.)

- Fixed the issue regarding use of sex and/or pgm as covariates for
  the X chromosome; such covariates are dropped just for the X
  chromosome.

- Revised plot.pxg() so that it no longer returns Rsq, fit and so
  forth.  I've gone back to returning just information about the data
  that are plotted.  I added an argument "main" in the case one wishes
  to use a plot title different from the default.

- Changed the default for plot.geno from horizontal=FALSE to
  horizontal=TRUE.

- switch.order() no longer resets the location of the first marker on
  the chromosome to 0 cM, but retains the offset from the original
  map.

- Changed "trace" as an argument to "verbose", to avoid clashes with
  the built-in R function trace()

- Revised summary.cross() to check the class of chromosomes (the
  components in cross$geno); they should each have class "A" or "X".

- Revised subset.cross() so that it won't produce duplicate
  chromosomes, and to treat missing values in the ind argument.

- Fixed a slight bug in fitqtl() regarding the drop one analysis
  when there is just one QTL.

- There was a slight bug in plot.pxg() regarding the X chromosome.

- Fixed a bug in scanone() regarding the X chromosome for RILs.

- Fixed a slight, silly bug in est.map() for sex-specific maps, in
  which the location of the initial marker was randomized.

- Revised plot.scanone(), plot.scantwo(), and other plotting functions
  so that they don't show their NULL return values.

- Fixed a bug in makeqtl() regarding the X chromosome if there are
  genotype probabilities.

- Revised the internal checkcovar() function to check that the chosen
  phenotype for scanone(), scantwo(), etc., is numeric.  Previously,
  an rather uninformative error message was given.

- Modified the convergence criterion for scantwo() by the EM
  algorithm; we now just look at the log likelihood, and not at the
  parameters.

- Fixed a bug in read.map.qtlcart().

- Fixed bugs in scanone and scantwo for method="imp", for the case
  that one has exactly one imputation.  Fixed a memory-overwrite
  problem in scanone() with method="imp".

- Fixed a bug in fixXgeno.f2; a slight error that shows up in the case
  the pgm values are switched.

- Added code to simulate Collaborative Cross data and to do
  calc.genoprob() and argmax.geno() on such data, though none of this
  is documented yet.

- movemarker() now updates the results of est.rf, calc.genoprob,
  calc.errorlod, sim.geno, argmax.geno, if they are available.
  (est.rf is simply re-run; the others are re-run for just the
  relevant chromosomes).

- calc.errorlod had neglected to include the map function as an
  attribute.  This is now fixed.

- Revised scantwo with method="em" so that it prints the verbose
  warning messages only if the verbose argument is > 1.  With
  verbose=TRUE, only the chromosomes are printed.


## Version 0.98, 9/11/2004:

### Major changes:

- There is no longer an argument "sep" for the function read.cross.
  Instead, we use "...", which is passed to the function read.table
  (all this just for the "csv" format).  (sep="," is still assumed for
  that format).  This change allows one to use sep=";" and dec=","
  which many people prefer.

- The function read.cross now automatically converts X chromosome
  genotype data into the standard internal format (with all
  individuals coded like an autosome in a backcross).

- The functions scanone() and scantwo() were revised so that they
  treat the X chromosome appropriately.  The argument "x.treatment"
  has been deleted.  What had been x.treatment="full" (namely, that
  hemizygous genotypes are considered different from homozygous
  genotypes) is now forced.  The big change concerns the "null
  hypothesis" for the X chromosome, which includes sex and/or "pgm" as
  covariates, in order to avoid spurious linkage on the X due to sex
  differences in the phenotype.  See the help files for these
  functions for details.  Note that the output of the scantwo function
  has changed somewhat; it would be best to re-run scantwo with this
  new version of the software.

- Added a function movemarker() for moving a marker from one
  chromosome to another.

- Revised scanone(), scantwo(), discan(), vbscan() and plot.info() so
  that inter-marker positions are cited as "c*.loc*" rather than
  "loc*.c*" or "loc*".  Added a function convert.scanone() for
  converting scanone output to the new format.  (The new plot.scantwo
  will interpret, for the old version of scanone output, every
  inter-marker position as a marker.)

### Minor changes:

- Added a function find.pheno() [from Brian Yandell] for finding the
  phenotype column with a particular name.

- Added a function find.flanking() [from Brian Yandell], which is
  similar to find.marker(), but gives not just the closest marker
  but also the left- and right-flanking markers.

- Added a function print.cross() which prints a short message and then
  calls summary.cross().  This was added in order to avoid the
  essentially always unintentional printing of an entire (generally
  quite large) cross object.

- Modified summary.scantwo to allow a type option to get summaries
  based on peak "joint" or "inter" and to allow negative thresholds
  relative to max joint, inter, or individual LODs [from Brian
  Yandell; I'm not really sure what this means].

- Modified plot.scantwo to modify the contour option, making 1.5 from
  the peak the default on each half-image, and to take a numeric
  vector of drops from the peaks.  Added arguments col.scheme and
  gamma for different color schemes [from Brian Yandell].

- Hugely sped up plot.scantwo for the case of lower="cond-int" and
  lower="cond-add"!  (By hugely, I mean by a factor of 50 or more.)

- Added a function summary.map() for giving summary information about
  a genetic map.

- In read.cross, if chromosome names start with "chr" or "chromosome"
  (ignoring case) these initial strings are removed.

- est.rf now returns a warning if a marker appears to have its
  genotypes miscoded (so that it shows a rec. frac. > 0.5 with
  LOD > 3).

- Deleted the function pull.chr(), which was deprecated in version
  0.89 in Nov, 2001.  You can use subset.cross() instead.

- Fixed a bug in calc.errorlod regarding the X chromosome in an
  intercross.

- read.cross with format="qtlcart" printed on screen the number of
  individuals but called it the number of phenotypes.  Brian Yandell
  provided another revision to this, to fix bugs regarding phenotypes
  that are factors.

- Fixed a slight bug in summary.cross() in the case that a cross
  contains no autosomal data.

- Fixed a bug in read.cross for the mapmaker format; if a marker or
  phenotype name were listed without any data, an error resulted.

- Fixed a bug in summary.scanone which resulted in multiple rows being
  returned for a single chromosome if multiple positions shared the
  same maximal LOD score.

- Fixed bugs in plot.scanone and plot.scantwo for the case that the
  "chr" argument was used with chromosomes not in the usual order.

- read.cross with format="csv" now prints a warning if unusual entries
  are seen in the genotype data.

- Revised the geno.table() function so that its first column is the
  chromosome number

- Modified read.cross.qtx(); "X" or "x" in a phenotype is a missing
  value.

- Changed a call to print.coefmat() to a call to printCoefmat(), as
  the former is being discarded in favor of the latter.

- Changed the name of the utility function fixXdata() to
  reviseXdata().  This function deals with the X chromosome genotypes
  in scanone() and scantwo().

- Modified the source code in hmm_bc.c, hmm_f2.c, hmm_main.h so that
  we don't repeatedly calculate log(0.5), log(2.0) and log(0.25), but
  rather rely on #define statements

- Got rid of the WIN32 stuff in addlog and subtractlog in util.c,
  which were used in version 0.97 as I'd not been able to get log1p to
  work.

- plot.missing now gives a more meaningful error if the argument
  "reorder" is greater than the number of phenotypes.

- I split up the read.cross.R file into several smaller files, for
  more easy revisions.

- Revised the function comparecrosses(), adding an argument "tol" so
  that the genetic maps and phenotypes do not need to be *exactly* the
  same, but can be identical to within the specified tolerance, "tol".
  Also, a warning (rather than an error) is produced if the
  inter-marker distances are the same but the position of the initial
  marker is different.

- Fixed a slight bug in write.cross for the case of phenotypes that
  are factors.

- Fixed a problem in read.cross for format "gary"; now I pass the
  "na.strings" argument for reading the phenotype data.

- Fixed a problem in read.cross with format "gary" or "csv" for
  ensuring that phenotypes that appear to be numeric are read as
  numeric and not as factors.

- Modified the function getsex(), which finds and interprets the sex
  and pgm columns in the phenotype data, for the case where sex is
  read as a factor with just one level (either "F", "f", "M", or "m").

- Added more understandable warning/error messages to plot.scanone if
  a chromosome ID in the "chr" argument doesn't match those in the
  scanone output.

- Fixed a problem in the example data "badorder"; chromosomes 2 and 3
  were switched.

- Fixed a bug in read.cross for the QTL Cartographer format for the
  case that there is just one individual.

- Modified some of the C calls to FORTRAN subroutines, using the
  macro F77_CALL(), as this is recommended in the R documentation.
  (Fortran routines currently used: dqrls, dpoco, dposl)

- Modified est.map so that it removes individuals with fewer than two
  typed markers.  Their presence is of no value in estimating the map,
  and can really slow things down.

- Revised a number of the help files to make the automated tests of
  the integrity of the software much faster.

- Revised read.cross.mm() so that (a) the cross type is taken as the
  4th word (rather than the last word) on the first line, in case
  someone includes additional characters, and (b) if the sample map
  format is used but chromosome assignments are not provided, a
  meaningful error message is displayed.

- Made a few very minor revisions to the "rqtltour.pdf" tutorial.

- Revised getgenonames() and reviseXdata() to remove the x.treatment
  argument, which is now assumed to be "full".  Revised effectplot()
  and plot.pxg() accordingly.

- Modified the utility function to do map expansion in RIs for the X
  chromosome.  Need to modify things further to take account of the
  lack of balance on the X chromosome...it's quite different from a
  backcross.

- Revised plot.pxg() so that it can take a vector of marker names
  [from Brian Yandell].

- Revised plot.scantwo to allow different color schemes and to give
  contours at 1.5 (or other specified values) below the maximum LOD
  [from Brian Yandell].

- Revised summary.scanone and print.summary.scanone so that
  summary.scanone will never print anything, but will leave it to
  print.summary.scanone to do so.  (Previously, summary.scanone
  printed a message if there were no peaks above the LOD threshold.)
  summary.scantwo and print.summary.scantwo were revised similarly.

- Modified plot.scanone to allow NAs in the LOD scores.

- Revised scanone and scantwo so that n.perm=0 is treated the same as
  if it were not provided.

- Fixed a bug in plot.rf() in the case that chromosomes are provided
  out of order, in which case the thing plotted garbage.  The fix
  involved revising subset.cross() so that the chr argument is sorted.

- Replaced the example data sets with compressed versions.

- Had to modify part of the C code for scanone with model="2part".
  I'd hard-coded some of the array limits!

- Changed the default line types and colors in plot.scanone().

- Added a warning for scanone() and scantwo() about the number of
  individuals that are omitted due to missing phenotype or covariate
  information.

- Revised write.cross so that X chromosome data is converted back
  from our internal format to the standard format.


## Version 0.97, 6/19/2003:

### Major changes:

- Added an argument "weights" to the functions scanone and scantwo, to
  allow differential weighting of individuals in a genome scan; these
  are only used with model="normal".

- Modified the function "plot.scantwo" for plotting the results from a
  two-dimensional, two-QTL genome scan.  There are two new arguments:
  lower controls what LOD scores are plot in the lower triangle
  (lower="joint" corresponds to the previous version of this
  function), while nodiag controls whether to plot the scanone results
  on the diagonal.  Using lower="cond-int" (the default) gets rid of
  the "coattail" effect often seen when plotting the joint LOD
  scores.  Ted Lystig for suggested this modification.

- Added a new function effectscan() for plotting the estimated allelic
  affect at all markers on selected chromosomes.

### Minor changes:

- Modified write.cross so that it will output data in "gary" format.

- Added a function lodint, for calculating LOD support intervals based
  on results from scanone.

- Added a function nmissing(), which calculates the total number of
  missing genotypes for each individual in a cross, or for each
  marker.

- Added a function pull.geno() for pulling out the set of genotype
  data for a cross as a single big matrix.

- Added a function comparecrosses() for verifying that two objects of
  class "cross" are identical.

- The results of geno.table() now includes P-values from chi-square
  tests for Mendelian segregation.

- Modified the function c.cross for combining crosses.  You can now
  combine backcrosses and intercrosses, provided that they have
  exactly the same genetic maps.  Further, we no longer discard the
  results of sim.geno and calc.genoprob, provided that the same step,
  off.end, and error.prob arguments were used.

- Added an additional argument, cex, to the function plot.geno, for
  control of the size of the points in the plot.  Also changed the
  orientation of the plot when horiz=FALSE, so that the centromere is
  at the top of the figure rather than the bottom.

- Fixed calc.pairprob so that it will work for RI lines ("risib" or
  "riself").  Thus, scantwo should work with RI lines now.

- Updated read.cross for format="gary" so that the marker positions
  file ("mapfile") and phenotype names file ("pnamesfile") are not
  necessary.  Set these arguments to NULL (e.g., mapfile=NULL) if the
  corresponding files are not available.

- Added a "chr" argument to max.scanone, so you can get the maximum
  LOD score for a particular chromosome.

- Revised the function switch.order() so that, if estimated
  recombination fractions are present (i.e., est.rf() was used), these
  are revised appropriately; previously they had been removed.
  Also added err and map.function arguments, to be passed to est.map()
  when the map is re-estimated.

- Revised scanone() and scantwo() slightly; the statement for
  producing a warning regarding the use of method=="im" (vs "imp" or
  "em") was slightly wrong.

- Fixed a slight bug in scanqtl() for the case that a fixed position
  is provided rather than a range (commented out two lines).

- summary.cross() now prints a warning if $data objects are data
  frames. (They should be simple matrices.)

- summary.cross() now prints a warning message if the genetic maps are
  not matrices with 2 rows for "f2ss" and "4way" crosses, or are
  matrices for other crosses.

- drop.markers() now prints a warning if some markers were not found.

- Added arguments ylim and add.legend to the function effectplot().

- Added arguments xlim and mtick to function plot.scanone().  (mtick
  allows marker locations to be indicated by triangles rather than
  line segments.)

- Fixed a bug in read.cross for the case that phenotypes have values
  like "1x2".

- Fixed a slight bug in write.cross for the qtlcart format.

- Fixed a bug in read.cross for the qtlcart format regarding the
  determination of whether a chromosome is autosomal or the X.
  (Previously, looked for an "X" or "x" in the marker names; now look
  at whether the chromosome names contains an "X" or "x".)

- Fixed a bug in makeqtl() for the case of a four-way cross.  (Hadn't
  dealt properly with sex-specific maps.

- fitqtl() now stops with a more meaningful error message if imputed
  genotypes are not available in the input "qtl" object.

- Revised the marker names for the X chromosome in the map10 dataset
  that is included.

- Revised est.map() for the case of a sex-specific f2 (cross type
  "f2ss"); the starting map for the EM algorithm is randomized a bit.

- Revised a bunch of the R code files so that paste() is not included
  within stop() or warn().

- In a couple of utility functions for the hidden Markov model engine,
  I need access to the log1p() function, but I'm having trouble with
  that in Windows.  Thus, in Windows only, I use log(1+x) in place of
  the preferred log1p() function.

- Added tests of input/output that are run when doing a check of the
  package.

## Version 0.96, 9/13/2002:

### Major changes:

- None.

### Minor changes:

- Added Listeria data in QTL Cartographer format to the sampledata
  directory.

- Revised read.cross and write.cross for QTL Cartographer format, so
  that the cross types are converted between those of R/qtl ("f2",
  "bc", "riself", "risib") and those of QTL Cartographer ("RI0",
  "RI1", "RI2", "B1", "B2", "SF2", "RF2").

- Revised read.cross for the Mapmaker format. The map file can now be
  in a second format, ".maps", which is created by Mapmaker/exp.  The
  function determines whether it has been presented with the .maps
  format or the 2- or 3-column tabular format that has been
  available.  Brian Yandell wrote the function to read ".maps" files.

- Fixed a small bug in write.cross for the Mapmaker format, and
  modified the ".prep" file that is created, so that marker distances
  are no longer included, and including lines "framework chr*".

- Revised read.cross with format="csv", so that it gives more clear
  error messages in some cases.

- Updated the R/qtl tutorial, rqtltour.pdf.  This is now in a
  directory "docs" in the R/qtl distribution.

## Version 0.95, 8/1/2002:

### Major changes:

- Modified the functions scanone and scantwo in order to treat the X
  chromosome appropriately.  Each has a new argument, x.treatment,
  which indicates how to treat the X chromosome (in particular,
  whether hemizygous males should be treated the same as homozygous
  females).  For analysis to proceed properly, there should be
  columns "sex" and "pgm" in the phenotype data, indicating the sex
  of each individual, and the direction of the cross.  See the X
  chromosome section of the help file for read.cross for more
  information.

- Added another argument to plot.scanone, "lodcolumn", an integer (or
  a vector of 3 integers) indicating which columns of the scanone
  output should be plotted (generally column 3).

- Added two functions for plotting phenotypes against marker
  genotypes.  plot.pxg() plots the phenotypes against the genotypes
  at a single marker.  effectplot() plots the average phenotypes
  against genotypes at one or two markers (or covariates).  Also
  added a function find.marker() which returns the name of the marker
  closest to a specified position.

- Added facilities for analyzing recombinant inbred lines.  We now
  allow two additional cross types, "riself" (RI lines from selfing)
  and "risib" (RI lines from sibling matings).  Added an internal
  function expand.rf.ri and made important modifications to
  calc.genoprob, sim.geno, argmax.geno, and est.map.  Also modified
  summary.cross, print.summary.cross, geno.table, replace.map,
  discan, ripple, scanone, scantwo, makeqtl, calc.errorlod, est.rf,
  write.cross.mm, write.cross.csv.

- Replaced the example fake.f2 with some new data, which includes
  both males and females and both directions of the intercross, in
  order to illustrate the proper analysis of the X chromosome.

- Modified read.cross and write.cross (and added code from Brian
  Yandell) to read and write data in QTL Cartographer format.

### Minor changes:

- Fixed a bug in read.cross; the "genotypes" argument needs to have
  "C" and "D" reversed. "C" = "not BB" = 5 internally; "D" = "not AA"
  = 4 internally.  Thanks to Martin Grandona for identifying the
  problem.

- Fixed a bug in read.cross for format "csv": an error occurred if
  marker positions were not given and the first individual was
  missing a phenotype.  Thanks to Justin Borevitz and Norman
  Warthmann for identifying the problem.

- Fixed a bug in fitqtl.  Also added type III sums of squares table
  and also nominal P-values.

- Revised the read.cross functions so that if the X chromosome data
  is coded as A:B, it gets re-coded appropriately.

- Revised read.cross.mm and read.cross.csv so that if marker
  positions are included, marker order is taken according to those
  positions.  (Previously, read.cross stopped with an error.)

- Added some additional example data, bristle3 and bristleX, from
  Long et al. (1995) Genetics 139:1273-1291.

- Added an example genetic map, map10, containing 19 autosomes and an
  X chromosome, with chromosome lengths approximately as in the
  mouse and markers at approximately 10 cM spacing.

- Changed web references "biosun01.biostat.jhsph.edu" to
  "www.biostat.jhsph.edu".

- Revised summary.cross to print an error if the cross type is not
  one of "f2", "f2ss", "bc", "4way", "risib", or "riself".

- Revised plot.map so that, when genetic maps are plotted vertically,
  the 1st marker is at the top (rather than at the bottom).


## Version 0.94, 5/30/2002:

### Major changes:

- Modified scanone and scantwo to include the methods "mr-imp" and
  "mr-argmax", for performing "marker regression" by first filling
  in any missing genotypes by one imputation ("mr-imp") or using the
  Viterbi algorithm ("mr-argmax").

- Edited functions read.cross.*, argmax.geno, sim.cross, and
  sim.draws, so that the data, argmax and draws portion of a cross
  object are stored as integers.  This can save considerable space.

- Added functions to perform a general scan by imputation:
  fitqtl(), makeqtl(), scanqtl().

### Minor changes:

- Edited the read.cross.* functions again; now by default dir = ""
  rather than ".", and I no longer remove any trailing "/" from dir.

- Added checks of genotype values in the function summary.cross.

- Rather than edit the read.cross function to read in X chromosome
  data appropriately, I instead edited its help file, to explain that
  X-linked data should be coded as an autosome in a backcross (with
  genotypes A and H).

- Fixed slight errors in the functions scanone, calc.genoprob,
  discan, calc.pairprob, and sim.geno regarding the naming of the
  genotypes for the 4-way cross.

- Changed a couple of apostrophes to double-quotes in the function
  summary.scantwo.

- Added a Morgan map function (mf.m and imf.m), with revisions to
  argmax.geno, calc.genoprob, calc.pairprob, calc.errorlod, est.map,
  ripple, sim.geno, sim.cross, sim.cross.bc, sim.cross.f2,
  sim.cross.4way, fill.geno.

- Fixed a slight bug in ripple regarding the estimated chromosome
  lengths in the case of a 4-way cross; I was picking out the wrong
  element of the map.

- Fixed a bug in plot.scantwo regarding chromosome labels.

- Fixed bugs in max.scantwo and scantwo.perm regarding infinite LOD
  scores (which comes up especially when running scantwo with
  method="mr").

- Fixed a bug in read.cross.qtx regarding the determination of
  whether a cross is a backcross or an intercross.  Also fixed the
  case of a backcross coded as H:B rather than A:H.

- Added a warning in summary.cross regarding duplicate markers.

- Modified the example cross data (such as hyper and listeria) so
  that genotypes are stored as integers.

- Updated the "README.txt" file, to include explanations for
  installation of R and R/qtl on Mac OS.

- Changed the default for the na.strings argument in read.cross and
  read.cross.csv to include "NA".

- Changed a couple of lines in write.cross.csv and write.cross.mm for
  the treatment of NA strings.

- Modified ripple so that orders considered are printed in a way that
  the left-most marker in the original order is always to the left
  of the right-most marker in the original order.

- In various functions, made sure that 0 < error.prob < 1.

- Edited plot.scanone so that when only one chromosome is plotted,
  the chromosome number doesn't appear at the top, and when multiple
  chromosomes are plotted, the chromosome numbers appear at the
  bottom, rather than cumulative cM position.

- Changed the addcov and intcov arguments to scanone and scantwo to
  addcovar and intcovar, respectively.


## Version 0.93, 4/1/2002:

### Major changes:

- Added ability to read data in Mapmanager QTX format. This may be
  done via the read.cross function by using the argument
  format="qtx".  Added a file in this format to the sampledata
  directory distributed with R/qtl.

- Modified function ripple(), which compares marker orders, so that
  it may evaluate counts of obligate crossovers, which will be
  extremely quick relative to performing an exact likelihood
  calculation.  This method has been made the default.

### Minor changes:

- Added functions max.scanone and max.scantwo for getting information
  on the location with the highest LOD or joint and interaction LODs

- Modified summary.scantwo so that if the argument thresholds has
  length 1, the interaction and conditional thresholds are assumed to
  be 0 (so that all chromosome pairs for which the maximum joint LOD
  is greater than the given threshold are printed).

- Revised the C functions emit_bc(), emit_f2(), emit_f2ss() and
  emit_4way() so that unexpected observed genotypes are treated as
  missing.

- Revised read.cross, read.cross.csv and read.cross.mm slightly, so
  that estimate.map is TRUE by default, and so that the genetic maps
  are re-estimated only if both estimate.map is TRUE and the genetic
  map is missing from the input files.  If estimate.map is FALSE and
  the genetic map is missing from the input files, a dummy genetic
  map is inserted.

- Fixed a bug in sim.cross, sorting the "model" matrix in advance of
  performing the simulation, because the results were erroneous if
  QTLs were specified out of order.

- Edited the functions read.cross.* to use the function file.path()
  to create file names.


## Version 0.92, 2/12/2002:

### Minor changes:

- In read.cross.mm and read.cross.csv, when using the function
  read.table, we replaced the use of as.is=TRUE with
  colClasses="character". Apparently as.is=TRUE didn't work in R
  version 1.4.0.

- In read.cross, changed the default of the argument "estimate.map"
  to FALSE.


## Version 0.91, 12/3/2001:

### Minor changes:

- Fixed a problem with chromosome labels in plot.scantwo.

- Fixed a slight bug in summary.ripple.

- Previously forgot to implement the use of the "main" arg for
  plot.scanone.

- Fixed a slight bug in read.cross.gary related to having just one
  marker on a chromosome.

- Fixed a slight bug in plot.cross for the case auto.layout=FALSE.

- Revised read.cross so that, for the csv format, if the argument
  "genotypes" is NULL, the genotypes are assumed to be correct. If
  there are genotypes > 5, it is assumed to be a 4-way cross.

- For some reason, the wrapper for est_map for 4-way crosses got
  deleted.  I've re-written it.  Hopefully it works!

- Fixed a slight bug in plot.map for plotting two sex-specific maps.
  (The function works by pulling apart the sex-specific maps and then
  calling plot.map again twice.  After those calls, it should return.)

- Expanded examples in the help file for fake.4way.

- Fixed a bug in create.map for sex-specific maps.

- Revised calc.genoprob, argmax.geno and sim.geno so that, in the case
  of one marker on a chromosome, off.end is forced to be > 0.

- Revised plot.scanone so that if there is exactly one LOD score for a
  chromosome, a small segment is plotted rather than a dot.

- Fixed a couple of minor bugs in read.cross for the mapmaker format:
  in dealing with the "symbols" information in the mapmaker file, and
  in counting the number of lines in the file.

- Added a utility function checkcovar() to check phenotypes and
  covariates in scanone and scantwo (thus removing some redundancy).


## Version 0.90, 11/24/2001:

### Minor changes:

- Replaced the example data fake.bc with something that will allow the
  illustration of the use of covariates.

- Added print.summary.ripple; I'd forgotten to write it before.

- Added an updated tutorial on R/qtl, distributed as the file
  rqtltour.pdf


## Version 0.89, 11/22/2001:

### Major changes:

- Consolidated scanone, vbscan and discan into the single function
  scanone, with an argument model=c("normal","binary","2part","np").
  The non-parametric "method" is now a "model".

- Buried scanone.perm and scantwo.perm as internal functions.  To do
  permutation tests, one now uses the main functions (such as scantwo)
  and specifies the n.perm argument.

- Similarly, read.cross.* and write.cross.* were buried, so that the
  user is expected to call either read.cross or write.cross rather
  than calling the format-specific functions directly.  This was done
  anticipating an increase in the number of such format-specific
  read.cross functions.

- Got rid of find.errors and plot.errors, as I don't like them.  Use
  calc.errorlod and plot.errorlod instead.

- Wanted to toss pull.chr, but instead just kept an internal version
  which calls subset.cross and prints a warning, in case our one
  official user has code which requires it.

### Minor changes:

- Added an "eq.spacing" argument to sim.map for generating maps with
  equally-spaced markers.  This seems more useful than putting them
  down at random.

- Re-wrote a great deal of the help documentation (especially the
  examples and details).

- Added a new example data set, badorder, with some errors in marker
  order.  (This is to illustrate the functions est.rf, ripple and
  switch.order.)

- Fixed a slight error in summary.scantwo.  We print pairs of loci
  only if their joint LOD exceeds its threshold and either (a) the
  epistasis LOD exceeds its threshold or (b) both conditional LODs
  exceed their thresholds.

- Totally re-wrote print.summary.scantwo.  It was unnecessarily
  complicated before.

- Made a very slight change regarding the zlim in plot.scantwo.

- Fixed scantwo, summary.scantwo and plot.scantwo to deal with cases
  of bad LOD scores (NAs, negative numbers and Infs).  A warning
  message will always be printed.

- Modified scanone_imp.c so that nullRss and altRss don't allocate
  memory each time.  Fixed a very bad bug in dealing with interactive
  covariates.  Fixed a single-character bug in scantwo_mr.c that was
  causing a core dump.

## Version 0.88, 11/20/2001:

### Major changes:

- Added a scantwo function to do two-dimensional genome scans,
  calculating LOD scores for a two-QTL model and to test epistasis
  between each pair, with calculations done by imputation, Haley-Knott
  regression, marker regression or the EM algorithm.  Hao Wu wrote the
  imputation method.

- With Hao Wu, wrote plot.scantwo to plot the results of scantwo,
  summary.scantwo to summarize the results, and scantwo.perm to get
  genome-wide LOD thresholds for a 2-dimensional genome scan by
  permutation tests.  The summary.scantwo function uses a criterion
  due to Gary Churchill and Saunak Sen.

- Added a C function to calculate joint genotype probabilities for
  pairs of putative QTLs on the same chromosome.  Because the
  resulting set of probabilities can take up a lot of memory, we're
  not going to make these accessible to the user.  The function
  calc.pairprob was created, but this is not to be called by the user,
  but rather will be called when needed.

### Minor changes:

- Added a "method" argument to vbscan, even though only method="em" is
  currently available.

- Revised scanone, scantwo, discan, vbscan, and their corresponding
  ".perm" functions so that the output has attribute "method" to
  indicate what method was used and attribute "type" to indicate the
  type of cross that was analyzed.

- Changed method="im" to method="em" in scanone and discan; changed
  method="markreg" again, this time to method="mg".  Changed the order
  of these methods in scanone.

- calc.genoprob now includes an attribute "map.function" with the
  probabilities.

- Changed colors plotted in plot.rf.

- Modified the C function scanone_mr (marker regression) to avoid
  repeatedly running the null model regression in the case of complete
  marker data.

- Changed a good amount of R code like "1:length(x)" to "seq(along=x)"

- Added a function fill.geno for imputing missing marker data by
  simulation (through sim.geno) or by the Viterbi algorithm (through
  argmax.geno), so that one may perform quick-and-dirty (with an
  emphasis on dirty) genome scans by marker regression.

- Fixed a small bug in sim.cross.f2.

- Fixed some problems related to chromosomes with only one marker:
  read.cross.csv, create.map, subset.cross.

- Fixed a bug in the location of chromosome labels in plot.scanone.
  Added an argument "main" for placing a title on the plot.

- Revised lots of little pieces of code using "drop=FALSE" when
  subsetting a matrix or array in order to retain the structure.

- read.cross.csv can now deal with categorical phenotypes, and
  plot.cross was revised to deal with such non-numeric phenotypes.
  Added an argument "auto.layout"; if TRUE, mfrow is set so that the
  many plots produced will all fit in one figure.  par(ask=TRUE) is no
  longer ever set.

- Revised sim.cross so that when keep.qtlgeno=TRUE, the QTL genotypes
  are retained in a component cross$qtlgeno (rather than within the
  data matrices).


## Version 0.87, 11/13/2001:

### Major changes:

- Hao Wu (hao@jax.org) has implemented the imputation method of Sen
  and Churchill (2001) for a genome scan, included as method="imp" in
  the function scanone.

- Added a non-parametric method to the function scanone, using a
  modified version of the Kruskal-Wallis test (cf Kruglyak and Lander
  1995).

- scanone now allows the use of covariates for all methods except the
  non-parametric method.

- Phenotypes in a cross object are now a data.frame.  Modified example
  data files and the following functions to make this work:
  sim.cross.*, read.cross.*, summary.cross, write.cross.csv.

### Minor changes:

- Changed the name of the "anova" method in scanone to "markreg".

- Changed the name of the argument "print.rf" in the est.map function
  to "trace."

- Modified the default cutoff in top.errorlod; allow cuts and colors
  in plot.errorlod to be specified by the user.

- summary.cross() now checks that markers are in increasing order.

- Made the third row (marker positions) in csv file optional in
  read.cross.csv.

- Added a utility function subset.cross() for pulling out specified
  chromosomes or groups of individuals from a cross object.  We should
  not need pull.chr() any longer.

- Added a utility function c.cross() for concatenating multiple cross
  objects.

- Changed stopping rules for discan, discan.perm, vbscan, vbscan.perm,
  est.map, est.rf, ripple, scanone, scanone.perm:
  |x(s+1) - x(s)| < e {|x(s)| + e*100} where by default e = 1e-4

- Fixed the utility function create.map() for the case where the
  genetic map starts at somewhere other than 0.

- Placed help information for discan.perm, scanone.perm and
  vbscan.perm within the files for discan, scanone and vbscan,
  respectively.


## Version 0.86, 11/4/2001:

- Fixed a *real* bug in argmax.geno().

- Added discan() for doing interval mapping with a dichotomous trait.

- Added documentation for the print.summary.* and internal functions.

- Edited documentation files to conform to R guidelines.

- Reduced the minimum value of the error.prob argument in est.map,
  calc.genoprob, argmax.geno and sim.geno from 1e-14 to 1e-50.


## Version 0.85, 10/29/2001:

- Tried to fix up some of the plot.* and summary.* functions so that I
  don't get warning messages in "R CMD check".

- Fixed a few minor problems in the help files.

- Updated the a.starting.point() help file.

- Fixed a couple of problems in marker order in the hyper data.

- Added plot.info() for plotting the proportion of missing information
  in the genotype data.

- Fixed bug in plot.scanone() that led to problems in overlaying LOD
  curves using add=TRUE.  Added an argument, gap, to specify the
  distance between chromosomes.

- Fixed bug in print.summary.scanone() that resulted in an error when
  there was just one chromosome with LOD above the specified
  threshold.


## Version 0.84, 10/10/2001:

- Fixed slight error in sim.cross(); marker genotypes were removed
  rather than qtl genotypes.  We now use the function drop.qtlgeno()
  to do this.

- Changed anova method in scanone() to use observed genotypes.
  Individuals with missing or partially missing genotypes are dropped.

- Added Haley-Knott regression method to scanone().

- Added a function ripple() for comparing marker orders for a single
  chromosome, looking at all permutations of a sliding window of
  markers.  Also added switch.order() to switch the order of markers
  on a specified chromosome.

- Removed null markers from listeria data.

- Fixed bugs in read.cross.mm() and write.cross.mm().

- Added csv and mapmaker format files to sample data directory.

- Allow specification of starting value is scanone and vbscan

- Added a document "rqtltour.pdf" describing the package and giving a
  couple of examples.

## Version 0.83, 09/23/2001:

- Fixed a very slight bug in summary.scanone().

- Changed the argument "which.chr" in plot.scanone() to simply "chr".

- Added a "chr" argument to plot.missing().


## Version 0.82, 09/20/2001:

- Added write.cross.csv(), for writing data in comma-delimited format.
  Changed write.mm() to write.cross.mm() and added write.cross() as a
  wrapper to these two functions.

- All functions that use map functions now allow use of the
  Carter-Falconer map function.

- Changed remove.markers(), remove.nullmarkers(), and remove.qtlgeno()
  to drop.markers(), drop.nullmarkers() and drop.qtlgeno().

- Revised plot.rf() so that missing values appear in gray.

- Added read.cross.gary(), to read data in Gary's format, and
  read.cross.csv(), to read data in comma-delimited format.

- Fixed the bugs in read.cross.mm(); see BUGS.txt.

- Fixed summary.cross() so that it checks marker names in the data and
  the map.

- Added summary.scanone(), giving a summary of the output of
  scanone().

- Added possibility of F2 intercross with sex-specific maps.  Use
  class "f2ss" rather than "f2."  This is in the testing stage.
  The only revised functions, at this point, are est.map() and
  calc.genoprob().

- Added a function convert2ss() to convert a cross object from "f2"
  to "f2ss" format.


## Version 0.81, 09/16/2001:

- plot.scanone can now plot three scanone outputs, and includes an
  "add" argument for adding additional outputs to a current plot.

- Replaced 1e-10 with 1e-14 as tolerance value for error probability
  and minimum map distance.

- Changed the "min.d" argument in plot.geno() to "min.sep", taken to
  be a percent of the chromosome length.

- Added Carter-Falconer map function: mf.cf() and imf.cf().  Note that
  there is no closed-form version of mf.cf(), and so I use the R
  function uniroot().

- Fixed a slight error in replace.map().

- In est.map, calculate the log likelihood at the end; this is saved
  as an attribute, "loglik" for each chromosome's map.  If the "print"
  argument is used, print the loglik, too.

- Made error.prob=0 the default for the functions argmax.geno(),
  calc.genoprob(), est.map(), and sim.geno().

- Fixed the file permissions for many of the files, so that they are
  readable by all users.


## Version 0.80, 08/07/2001:

- Eliminated the map component of the results of calc.genoprob,
  argmax.geno, and sim.geno.  Since we are now including attributes
  "error.prob," "step," and "off.end," we can just use create.map() to
  recreate the map each time, without having to carry it along.

- Changed the name of plot.geno() to plot.missing() and plot.chr() to
  plot.geno().

- Added vbscan() and vbscan.perm() to perform the analysis described
  in V Boyartchuk et al. (2001), for a phenotype where some
  individuals have some quantitative phenotype, while for others it is
  undefined.  (Examples: the size of a lesion, where some individuals
  exhibit no lesion; time-to-death after an infection, where some
  individuals recover from the infection.)


## Version 0.79, 07/27/2001:

- Added map functions (and inversion map functions) for Haldane and
  Kosambi, so that I'm not re-creating them all of the time within
  functions.

- Added a function plot.chr() to plot genotypes for a specific
  chromosome, with likely errors (as determined by calc.errorlod() or
  find.errors()) highlighted.

- Added a warning to the help file for argmax.geno.  The results
  greatly depend on the value of the step argument, and may not be
  terribly trustworthy.  Also, if several sequences (of underlying
  genotypes) are all most likely, our method of randomly choosing
  among them is not right...recombination events are too far to the
  right.


## Version 0.78, 07/24/2001:

- Fixed a small bug in create.map(), which is used by
  calc.genoprob().  An error occurred in the case of a genetic
  map with equally spaced markers, when the argument "step" was set to
  be exactly the inter-marker distance.

- Modified calc.genoprob(), calc.argmax(), sim.geno() and
  calc.errorlod() so that their corresponding components have
  attributes "error.prob", "step" and "off.end" (only "error.prob" for
  calc.errorlod()), specifying the corresponding values used in the
  calculations.  Modified calc.errorlod() to re-run calc.genoprob() if
  the error.prob attribute is different from the corresponding
  argument.


## Version 0.77, 06/22/2001:

- Fixed a small bug in sim.cross(), where dimnames of error component
  was wrong, when simulating genotyping errors with a QTL present.


## Version 0.76, 05/17/2001:

- This is a totally revised version of the package.  Most
  importantly, the data structure for a cross has completely
  changed.  The function convert.cross is included, for converting
  data from the old structure to the new one.  See the help file for
  read.cross for a description of the new data structure.

- The main hidden Markov model engine has been rewritten, to make
  things more flexible and general.  We've now implemented the Viterbi
  algorithm, in the function argmax.geno, to calculate the most likely
  sequence of underlying genotypes, given the observed marker data,
  and we've fixed the calculation of the Lincoln and Lander error LOD
  scores.  The analysis of phase-known four-way crosses is now
  possible.

- The "singlescan" function (to do a genome scan with a single QTL
  model) is now called "scanone" (to save a few keystrokes).  Note
  that this function does not yet allow the use of covariates.  We'll
  add that feature in the near future.

- Saunak Sen and I are now working together on this project, and so
  things will begin to progress more quickly (we hope).
