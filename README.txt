----------------------------------------------------------------------
R/qtl
copyright (c) 2001-2012, Karl W Broman
http://www.rqtl.org

Authors: Karl W Broman <kbroman@biostat.wisc.edu> and Hao Wu, with
         ideas from Gary Churchill and Saunak Sen and contributions
         from Danny Arends, Ritsert Jansen, Pjotr Prins, and Brian
         Yandell
----------------------------------------------------------------------

R/qtl is an extensible, interactive environment for mapping
quantitative trait loci (QTL) in experimental crosses. It is
implemented as an add-on package for the freely available and widely
used statistical language/software R (see http://www.R-project.org). 
The development of this software as an add-on to R allows us to take
advantage of the basic mathematical and statistical functions, and
powerful graphics capabilities, that are provided with R. Further, the
user will benefit by the seamless integration of the QTL mapping
software into a general statistical analysis program. Our goal is to
make complex QTL mapping methods widely accessible and allow users to
focus on modeling rather than computing.

A key component of computational methods for QTL mapping is the hidden
Markov model (HMM) technology for dealing with missing genotype
data. We have implemented the main HMM algorithms, with allowance for
the presence of genotyping errors, for backcrosses, intercrosses, and
phase-known four-way crosses.

The current version of R/qtl includes facilities for estimating
genetic maps, identifying genotyping errors, and performing single-QTL
genome scans and two-QTL, two-dimensional genome scans, by interval
mapping (with the EM algorithm), Haley-Knott regression, and multiple
imputation. All of this may be done in the presence of covariates
(such as sex, age or treatment). One may also fit higher-order QTL
models by multiple imputation and Haley-Knott regression.

----------------------------------------------------------------------
The R/qtl package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
http://www.r-project.org/Licenses/GPL-3
----------------------------------------------------------------------
