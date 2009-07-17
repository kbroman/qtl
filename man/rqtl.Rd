\name{rqtl}
\alias{rqtl}

\title{R/qtl}

\description{
}

\section{Introduction}{
}

\references{
  Boyartchuk V. L., Broman, K. W., Mosher, R. E., D'Orazio
  S. E. F., Starnbach, M. N. and Dietrich, W. F. (2001) Multigenic
  control of \emph{Listeria monocytogenes} susceptibility in
  mice. \emph{Nature Genetics} \bold{27}, 259--260.
  
  Broman,  K. W. (2003) Mapping quantitative trait loci in the case of a
  spike in the phenotype distribution. \emph{Genetics} \bold{163},
  1169--1175. 

  Broman, K. W., Sen, \'S, Owens, S. E., Manichaikul, A.,
  Southard-Smith, E. M. and Churchill G. A.  The X chromosome in
  quantitative trait locus mapping.  \emph{Genetics}, to appear.

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

  Sen, \'S. and Churchill, G. A. (2001) A statistical framework for quantitative
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
hyper <- calc.genoprob(hyper, step=2.5)
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
hyper.c1n4 <- sim.geno(subset(hyper, chr=c(1,4)),
                       step=2.5, n.draws=8)
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
}\dontshow{permo <- scanone(hyper, method="hk", n.perm=100)
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
}\dontshow{operm.strat <- scanone(hyper, method="hk", n.perm=100,
                       perm.strata=extremes)
}

summary(operm.strat)



####################
# X-specific permutations
####################
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"))}
fake.f2 <- calc.genoprob(fake.f2, step=2.5)

# genome scan
out <- scanone(fake.f2, method="hk")

# X-chr-specific permutations
\dontrun{operm <- scanone(fake.f2, method="hk", n.perm=1000, perm.Xsp=TRUE)
}\dontshow{operm <- scanone(fake.f2, method="hk", n.perm=100, perm.Xsp=TRUE)
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
listeria <- calc.genoprob(listeria,step=2.5)
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
plot(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)

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

\seealso{ \code{\link[qtl]{plot.scanone}}, 
  \code{\link[qtl]{summary.scanone}}, \code{\link[qtl]{scantwo}},
  \code{\link[qtl]{calc.genoprob}}, \code{\link[qtl]{sim.geno}},
  \code{\link[qtl]{max.scanone}},
  \code{\link[qtl]{summary.scanoneperm}},
  \code{\link[qtl]{-.scanone}}, \code{\link[qtl]{+.scanone}} }

\keyword{models}
