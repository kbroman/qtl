\name{scanone}
\alias{scanone}

\title{Genome scan with a single QTL model}

\description{
  Genome scan with a single QTL model, with possible allowance for
  covariates, using any of several possible models for the phenotype and
  any of several possible numerical methods.
}

\usage{
scanone(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
        method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
        addcovar=NULL, intcovar=NULL, weights=NULL,
        use=c("all.obs", "complete.obs"), upper=FALSE,
        ties.random=FALSE, start=NULL, maxit=4000,
        tol=1e-4, n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose,
        batchsize=250, n.cluster=1, ind.noqtl)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes for which LOD
   scores should be calculated.  This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This can be a vector of integers; for methods
   \code{"hk"} and \code{"imp"} this can be considerably faster than doing
   them one at a time.  One may also give a character strings matching
   the phenotype names.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{model}{The phenotype model: the usual normal model, a model for
   binary traits, a two-part model or non-parametric analysis}
 \item{method}{Indicates whether to use the EM algorithm, 
   imputation, Haley-Knott regression, the extended Haley-Knott method,
   or marker regression.  Not all methods are available for all models.
   Marker regression is performed either by dropping individuals with
   missing genotypes (\code{"mr"}), or by first filling in missing data
   using a single imputation (\code{"mr-imp"}) or by the Viterbi
   algorithm (\code{"mr-argmax"}).}
 \item{addcovar}{Additive covariates;
   allowed only for the normal and binary models.}
 \item{intcovar}{Interactive covariates (interact with QTL genotype);
   allowed only for the normal and binary models.}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights.  Used only
   in the case \code{model="normal"}.}
 \item{use}{In the case that multiple phenotypes are selected to be
   scanned, this argument indicates whether to use all individuals, 
   including those missing some phenotypes, or just those individuals
   that have data on all selected phenotypes.}
 \item{upper}{Used only for the two-part model; if true, the
   "undefined" phenotype is the maximum observed phenotype; otherwise,
   it is the smallest observed phenotype.} 
 \item{ties.random}{Used only for the non-parametric "model"; if TRUE,
   ties in the phenotypes are ranked at random.  If FALSE, average ranks
   are used and a corrected LOD score is calculated.}
 \item{start}{Used only for the EM algorithm with the normal model and
   no covariates.  If \code{NULL}, use the usual starting values; if
   length 1, use random initial weights for EM; otherwise, this should
   be a vector of length n+1 (where n is the number of possible
   genotypes for the cross), giving the initial values for EM.}
 \item{maxit}{Maximum number of iterations for methods \code{"em"} and
   \code{"ehk"}.}
 \item{tol}{Tolerance value for determining convergence for methods
   \code{"em"} and \code{"ehk"}.} 
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
 \item{perm.Xsp}{If \code{n.perm} > 0, so that a permutation test will
   be performed, this indicates whether separate permutations should be
   performed for the autosomes and the X chromosome, in order to get an
   X-chromosome-specific LOD threshold.  In this case, additional
   permutations are performed for the X chromosome.}
 \item{perm.strata}{If \code{n.perm} > 0, this may be used to perform a
   stratified permutation test.  This should be a vector with the same
   number of individuals as in the cross data.  Unique values indicate
   the individual strata, and permutations will be performed within the
   strata.}
 \item{verbose}{In the case \code{n.perm} is specified, display
   information about the progress of the permutation tests.}
 \item{batchsize}{The number of phenotypes (or permutations) to be run
   as a batch; used only for methods \code{"hk"} and \code{"imp"}.}
 \item{n.cluster}{If the package \code{snow} is available and
   \code{n.perm} > 0, permutations are run in parallel using this number
   of nodes.}
 \item{ind.noqtl}{Indicates individuals who should not be allowed a QTL
   effect (used rarely, if at all); this is a logical vector of same
   length as there are individuals in the cross.}
}

\section{Models}{
  \bold{The normal model} is the standard model for QTL mapping (see
  Lander and Botstein 1989).  The
  residual phenotypic variation is assumed to follow a normal
  distribution, and analysis is analogous to analysis of variance.

  \bold{The binary model} is for the case of a binary phenotype, which
  must have values 0 and 1.  The proportions of 1's in the different
  genotype groups are compared.  Currently only methods \code{em}, \code{hk}, and
  \code{mr} are available for this model.  See Xu and Atchley (1996) and
  Broman (2003).

  \bold{The two-part model} is appropriate for the case of a spike in
  the phenotype distribution (for example, metastatic density when many
  individuals show no metastasis, or survival time following an
  infection when individuals may recover from the infection and fail to
  die).  The two-part model was described by 
  Boyartchuk et al. (2001) and Broman (2003).  Individuals with QTL
  genotype \eqn{g} have probability \eqn{p_g}{p[g]} of having an
  undefined phenotype (the spike), while if their phenotype is defined,
  it comes from a normal distribution with mean \eqn{\mu_g}{mu[g]} and
  common standard deviation \eqn{\sigma}{s}. Three LOD scores are
  calculated: LOD(\eqn{p,\mu}{p,mu}) is for the test of the hypothesis
  that \eqn{p_g = p}{p[g] = p} and \eqn{\mu_g = \mu}{mu[g] = mu}.
  LOD(\eqn{p}) is for the test that \eqn{p_g = p}{p[g] = p} while the
  \eqn{\mu_g}{mu[g]} may vary. LOD(\eqn{\mu}{mu}) is for the test that
  \eqn{\mu_g = \mu}{mu[g] = mu} while the \eqn{p_g}{p[g]} may vary. 

  \bold{With the non-parametric "model"}, an extension of the
  Kruskal-Wallis test is used; this is similar to the method described
  by Kruglyak and Lander (1995).  In the case of incomplete genotype
  information (such as at locations between genetic markers), the
  Kruskal-Wallis statistic is modified so that the rank for each
  individual is weighted by the genotype probabilities, analogous to
  Haley-Knott regression.  For this method, if the argument
  \code{ties.random} is TRUE, ties in the phenotypes are assigned random
  ranks; if it is FALSE, average ranks are used and a corrected LOD
  score is calculate.  Currently the \code{method} argument is ignored
  for this model.
}

\section{Methods}{
  \bold{\code{em}}: maximum likelihood is performed via the
  EM algorithm (Dempster et al. 1977), first used in this context by
  Lander and Botstein (1989).

  \bold{\code{imp}}: multiple imputation is used, as described by Sen
  and Churchill (2001).

  \bold{\code{hk}}: Haley-Knott regression is used (regression of the
  phenotypes on the multipoint QTL genotype probabilities), as described
  by Haley and Knott (1992).

  \bold{\code{ehk}}: the extended Haley-Knott method is used (like H-K,
  but taking account of the variances), as described in Feenstra et
  al. (2006).

  \bold{\code{mr}}: Marker regression is used.  Analysis is performed
  only at the genetic markers, and individuals with missing genotypes
  are discarded.  See Soller et al. (1976).
}

\section{Covariates}{
Covariates are allowed only for the normal and binary models.  The
normal model is \eqn{y = \beta_q + A \gamma + Z \delta_q + \epsilon}{y = b[q] +
  A g + Z d[q] + e} where \emph{q} is the unknown QTL genotype, \emph{A}
is a matrix of additive covariates, and \emph{Z} is a matrix of
covariates that interact with the QTL genotype.  The columns of \emph{Z}
are forced to be contained in the matrix \emph{A}.  The binary model is
the logistic regression analog.

The LOD score is calculated comparing the likelihood of the above
model to that of the null model \eqn{y = \mu + A \gamma + \epsilon}{y =
  m + A g + e}.

Covariates must be numeric matrices.  Individuals with any missing
covariates are discarded.  
}

\section{X chromosome}{
The X chromosome must be treated specially in QTL mapping.  See Broman
et al. (2006).

If both males and females are included, male hemizygotes are allowed
to be different from female homozygotes.  Thus, in a backcross, we will
fit separate means for the genotype classes AA, AB, AY, and BY.  In such
cases, sex differences in the phenotype could cause spurious linkage to
the X chromosome, and so the null hypothesis must be changed to allow
for a sex difference in the phenotype.

Numerous special cases must be considered, as detailed in the following
table. 

\tabular{lllccc}{
\bold{BC}
 \tab \tab \bold{Sexes} \tab \bold{Null} \tab \bold{Alternative} \tab \bold{df} \cr
 \tab \tab both sexes   \tab sex         \tab AA/AB/AY/BY        \tab 2         \cr 
 \tab \tab all female   \tab grand mean  \tab AA/AB              \tab 1         \cr
 \tab \tab all male     \tab grand mean  \tab AY/BY              \tab 1         \cr
 \tab \tab              \tab             \tab                    \tab           \cr
 
\bold{F2}
 \tab \bold{Direction} \tab \bold{Sexes}\tab \bold{Null}\tab \bold{Alternative} \tab \bold{df} \cr
 \tab Both             \tab both sexes  \tab femaleF/femaleR/male \tab AA/ABf/ABr/BB/AY/BY \tab 3         \cr
 \tab                  \tab all female  \tab pgm        \tab AA/ABf/ABr/BB      \tab 2         \cr
 \tab                  \tab all male    \tab grand mean \tab AY/BY              \tab 1         \cr
 \tab Forward          \tab both sexes  \tab sex        \tab AA/AB/AY/BY        \tab 2         \cr
 \tab                  \tab all female  \tab grand mean \tab AA/AB              \tab 1         \cr
 \tab                  \tab all male    \tab grand mean \tab AY/BY              \tab 1         \cr
 \tab Backward         \tab both sexes  \tab sex        \tab AB/BB/AY/BY        \tab 2         \cr
 \tab                  \tab all female  \tab grand mean \tab AB/BB              \tab 1         \cr
 \tab                  \tab all male    \tab grand mean \tab AY/BY              \tab 1         \cr
}

In the case that the number of degrees of freedom for the linkage test
for the X chromosome is different from that for autosomes, a separate
X-chromosome LOD threshold is recommended.  Autosome- and
X-chromosome-specific LOD thresholds may be estimated by permutation
tests with \code{scanone} by setting \code{n.perm}>0 and using
\code{perm.Xsp=TRUE}. 

}

\details{
  Use of the EM algorithm, Haley-Knott regression, and the extended
  Haley-Knott method require that multipoint genotype probabilities are
  first calculated using \code{\link{calc.genoprob}}.  The
  imputation method uses the results of \code{\link{sim.geno}}.

  Individuals with missing phenotypes are dropped.

  In the case that \code{n.perm}>0, so that a permutation
  test is performed, the R function \code{scanone} is called repeatedly.
  If \code{perm.Xsp=TRUE}, separate permutations are performed for the
  autosomes and the X chromosome, so that an X-chromosome-specific
  threshold may be calculated.  In this case, \code{n.perm} specifies
  the number of permutations used for the autosomes; for the X
  chromosome, \code{n.perm} \eqn{\times \, L_A/L_X}{* L_A/L_X} permutations
  will be run, where \eqn{L_A} and \eqn{L_X} are the total genetic
  lengths of the autosomes and X chromosome, respectively.  More
  permutations are needed for the X chromosome in order to obtain
  thresholds of similar accuracy.  
 
  For further details on the models, the methods and the use of
  covariates, see below.
}

\value{
  If \code{n.perm} is missing, the function returns a data.frame whose
  first two columns contain the chromosome IDs and cM positions.
  Subsequent columns contain the LOD scores for each phenotype.
  In the case of the two-part model, there are three LOD score columns
  for each phenotype: LOD(\eqn{p,\mu}{p,mu}), LOD(\eqn{p}) and
  LOD(\eqn{\mu}{mu}).   The result is given class \code{"scanone"} and
  has attributes  \code{"model"}, \code{"method"}, \code{"df"} and
  \code{"type"} (the latter is the type of cross analyzed). 

  If \code{n.perm} is specified, the function returns the results of a
  permutation test and the output has class \code{"scanoneperm"}.  If
  \code{perm.Xsp=FALSE}, the function returns a matrix with
  \code{n.perm} rows, each row containing the genome-wide 
  maximum LOD score for each of the phenotypes.  In the case of the
  two-part model, there are three columns for each phenotype,
  corresponding to the three different LOD scores. If
  \code{perm.Xsp=TRUE}, the result contains separate permutation results
  for the autosomes and the X chromosome respectively, and an attribute
  indicates the lengths of the chromosomes and an indicator of which
  chromosome is X.
}

\references{
  Boyartchuk,  V. L., Broman, K. W., Mosher, R. E., D'Orazio
  S. E. F., Starnbach, M. N. and Dietrich, W. F. (2001) Multigenic
  control of \emph{Listeria monocytogenes} susceptibility in
  mice. \emph{Nature Genetics} \bold{27}, 259--260.
  
  Broman,  K. W. (2003) Mapping quantitative trait loci in the case of a
  spike in the phenotype distribution. \emph{Genetics} \bold{163},
  1169--1175. 

  Broman, K. W., Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}, Owens, S. E., Manichaikul, A.,
  Southard-Smith, E. M. and Churchill G. A. (2006)  The X chromosome in
  quantitative trait locus mapping.  \emph{Genetics}, \bold{174}, 2151--2158.

  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.

  Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum
  likelihood from incomplete data via the EM algorithm.  \emph{J. Roy.
    Statist. Soc.} B, \bold{39}, 1--38.

  Feenstra, B., Skovgaard, I. M. and Broman, K. W. (2006) Mapping
  quantitative trait loci by an extension of the Haley-Knott regression
  method using estimating equations. \emph{Genetics}, \bold{173},
  2111--2119. 

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Kruglyak, L. and Lander, E. S. (1995) A nonparametric approach for
  mapping quantitative trait loci.  \emph{Genetics} \bold{139},
  1421--1428. 

  Lander, E. S. and Botstein, D. (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  \emph{Genetics}
  \bold{121}, 185--199.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.

  Soller, M., Brody, T. and Genizi, A. (1976) On the power of experimental
  designs for the detection of linkage between marker loci and
  quantitative loci in crosses between inbred lines.
  \emph{Theor. Appl. Genet.} \bold{47}, 35--39. 

  Xu, S., and Atchley, W.R. (1996) Mapping quantitative trait loci for
  complex binary diseases using line crosses. \emph{Genetics}
  \bold{143}, 1417--1424.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}; Hao Wu}

\examples{
###################
# Normal Model
###################
data(hyper)
\dontshow{hyper <- subset(hyper,chr=1:4,ind=1:100)}
# Genotype probabilities for EM and H-K
\dontrun{hyper <- calc.genoprob(hyper, step=2.5)
}\dontshow{hyper <- calc.genoprob(hyper)}
out.em <- scanone(hyper, method="em")
out.hk <- scanone(hyper, method="hk")

# Summarize results: peaks above 3
summary(out.em, thr=3)
summary(out.hk, thr=3)

# An alternate method of summarizing:
#     patch them together and then summarize
out <- c(out.em, out.hk)
summary(out, thr=3, format="allpeaks")

# Plot the results
plot(out.hk, out.em)
plot(out.hk, out.em, chr=c(1,4), lty=1, col=c("blue","black"))

# Imputation; first need to run sim.geno
# Do just chromosomes 1 and 4, to save time
\dontrun{hyper.c1n4 <- sim.geno(subset(hyper, chr=c(1,4)),
                       step=2.5, n.draws=8)
}\dontshow{hyper.c1n4 <- sim.geno(subset(hyper, chr=c(1,4)),
                       step=0, n.draws=4)}
out.imp <- scanone(hyper.c1n4, method="imp")
summary(out.imp, thr=3)

# Plot all three results
plot(out.imp, out.hk, out.em, chr=c(1,4), lty=1,
     col=c("red","blue","black"))

# extended Haley-Knott
out.ehk <- scanone(hyper, method="ehk")
plot(out.hk, out.em, out.ehk, chr=c(1,4))

# Permutation tests
\dontrun{permo <- scanone(hyper, method="hk", n.perm=1000)
}\dontshow{permo <- scanone(hyper, method="hk", n.perm=5)
}

# Threshold from the permutation test
summary(permo, alpha=c(0.05, 0.10))

# Results above the 0.05 threshold
summary(out.hk, perms=permo, alpha=0.05)


####################
# scan with square-root of phenotype
#   (Note that pheno.col can be a vector of phenotype values)
####################
out.sqrt <- scanone(hyper, pheno.col=sqrt(pull.pheno(hyper, 1)))
plot(out.em - out.sqrt, ylim=c(-0.1,0.1),
     ylab="Difference in LOD")
abline(h=0, lty=2, col="gray")


####################
# Stratified permutations
####################
extremes <- (nmissing(hyper)/totmar(hyper) < 0.5)

\dontrun{operm.strat <- scanone(hyper, method="hk", n.perm=1000,
                       perm.strata=extremes)
}\dontshow{operm.strat <- scanone(hyper, method="hk", n.perm=5,
                       perm.strata=extremes)
}

summary(operm.strat)



####################
# X-specific permutations
####################
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"))}
\dontrun{fake.f2 <- calc.genoprob(fake.f2, step=2.5)
}\dontshow{fake.f2 <- calc.genoprob(fake.f2, step=0)}

# genome scan
out <- scanone(fake.f2, method="hk")

# X-chr-specific permutations
\dontrun{operm <- scanone(fake.f2, method="hk", n.perm=1000, perm.Xsp=TRUE)
}\dontshow{operm <- scanone(fake.f2, method="hk", n.perm=5, perm.Xsp=TRUE)
}

# thresholds
summary(operm)

# scanone summary with p-values
summary(out, perms=operm, alpha=0.05, pvalues=TRUE)



###################
# Non-parametric
###################
out.np <- scanone(hyper, model="np")
summary(out.np, thr=3)

# Plot with previous results
plot(out.np, chr=c(1,4), lty=1, col="green")
plot(out.imp, out.hk, out.em, chr=c(1,4), lty=1,
     col=c("red","blue","black"), add=TRUE)

###################
# Two-part Model
###################
data(listeria)
\dontshow{listeria <- subset(listeria,chr=c(1,5,13),ind=!is.na(listeria$pheno[,1]))}
\dontrun{listeria <- calc.genoprob(listeria,step=2.5)
}\dontshow{listeria <- calc.genoprob(listeria,step=0)}
out.2p <- scanone(listeria, model="2part", upper=TRUE)
summary(out.2p, thr=c(5,3,3), format="allpeaks")

# Plot all three LOD scores together
plot(out.2p, out.2p, out.2p, lodcolumn=c(2,3,1), lty=1, chr=c(1,5,13),
     col=c("red","blue","black"))

# Permutation test
\dontrun{permo <- scanone(listeria, model="2part", upper=TRUE,
                 n.perm=1000)
}\dontshow{permo <- scanone(listeria, model="2part", upper=TRUE,
                 n.perm=3)
}

# Thresholds
summary(permo)

###################
# Binary model
###################
binphe <- as.numeric(pull.pheno(listeria,1)==264)
out.bin <- scanone(listeria, pheno.col=binphe, model="binary")
summary(out.bin, thr=3)

# Plot LOD for binary model with LOD(p) from 2-part model
plot(out.bin, out.2p, lodcolumn=c(1,2), lty=1, col=c("black", "red"),
     chr=c(1,5,13))

# Permutation test
\dontrun{permo <- scanone(listeria, pheno.col=binphe, model="binary",
                 n.perm=1000)
}\dontshow{permo <- scanone(listeria, pheno.col=binphe, model="binary",
                 n.perm=3)
}

# Thresholds
summary(permo)

###################
# Covariates
###################
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(2,5,10))}
\dontrun{fake.bc <- calc.genoprob(fake.bc, step=2.5)
}\dontshow{fake.bc <- calc.genoprob(fake.bc, step=0)}

# genome scans without covariates
out.nocovar <- scanone(fake.bc)

# genome scans with covariates
ac <- pull.pheno(fake.bc, c("sex","age"))
ic <- pull.pheno(fake.bc, "sex")

out.covar <- scanone(fake.bc, pheno.col=1,
                     addcovar=ac, intcovar=ic)
summary(out.nocovar, thr=3)
summary(out.covar, thr=3)
plot(out.covar, out.nocovar, chr=c(2,5,10))
}

\seealso{ \code{\link{plot.scanone}}, 
  \code{\link{summary.scanone}}, \code{\link{scantwo}},
  \code{\link{calc.genoprob}}, \code{\link{sim.geno}},
  \code{\link{max.scanone}},
  \code{\link{summary.scanoneperm}},
  \code{\link[=-.scanone]{-.scanone}}, \code{\link[=+.scanone]{+.scanone}} }

\keyword{models}
