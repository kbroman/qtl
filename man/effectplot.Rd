\name{effectplot}
\alias{effectplot}

\title{Plot phenotype means against genotypes at one or two markers}

\description{
  Plot the phenotype means for each group defined by the genotypes at
  one or two markers (or the values at a discrete covariate).
}

\usage{
effectplot(cross, pheno.col=1, mname1, mark1, geno1, mname2, mark2,
           geno2, main, ylim, xlab, ylab, col, add.legend=TRUE,
           legend.lab, draw=TRUE, var.flag=c("pooled","group"))
}

\arguments{
  \item{cross}{An object of class \code{cross}.}
  \item{pheno.col}{Column number in the phenotype matrix to be
    drawn in the plot.  One may also give a character string matching
   a phenotype name. Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
  \item{mname1}{Name for the first marker or pseudomarker.
  Pseudomarkers (that is, non-marker positions on the imputation grid)
  may be referred to in a form like \code{"5@30.3"}, for position 30.3
  on chromosome 5.}
  \item{mark1}{Genotype data for the first marker.  If unspecified,
    genotypes will be taken from the data in the input cross object,
    using the name specified in \code{mname1}.}
  \item{geno1}{Optional labels for the genotypes (or classes in a covariate).}
  \item{mname2}{Name for the second marker or pseudomarker (optional).}
  \item{mark2}{Like \code{mark1} (optional).}
  \item{geno2}{Optional labels for the genotypes (or classes in a covariate).}
  \item{main}{Optional figure title.}
  \item{ylim}{Optional y-axis limits.}
  \item{xlab}{Optional x-axis label.}
  \item{ylab}{Optional y-axis label.}
  \item{col}{Optional vector of colors for the different line segments.}
  \item{add.legend}{A logical value to indicate whether to add a legend.}
  \item{legend.lab}{Optional title for the legend.}
  \item{draw}{A logical value to indicate generate the plot or not. If
    FALSE, no figure will be plotted and this function can be used to
    calculate the group means and standard errors.}
  \item{var.flag}{The method to calculate the group variance. "pooled"
    means to use the pooled variance and "group" means to calculate from
    individual group.}
}

\value{
  A data.frame containing the phenotype means and standard errors
  for each group.
}

\details{
  In the plot, the y-axis is the phenotype.  In the case of one marker,
  the x-axis is the genotype for that marker. In the case of two
  markers, the x-axis is for different genotypes of the second marker,
  and the genotypes of first marker are represented by lines in
  different colors.  Error bars are plotted at \eqn{\pm}{+/-} 1 SE.

  The results of \code{sim.geno} are used; if they are not available,
  \code{sim.geno} is run with \code{n.draws=16}.  The average phenotype
  for each genotype group takes account of missing genotype data by
  averaging across the imputations.  The SEs take account of both the
  residual phenotype variation and the imputation error.
}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(1, 13, "X"))}

# impute genotype data
\dontrun{fake.f2 <- sim.geno(fake.f2, step=5, n.draws=64)
}\dontshow{fake.f2 <- sim.geno(fake.f2, step=5, n.draws=8)
}

########################################
# one marker plots
########################################
### plot of genotype-specific phenotype means for 1 marker
mname <- find.marker(fake.f2, 1, 37) # marker D1M437
effectplot(fake.f2, pheno.col=1, mname1=mname)

### output of the function contains the means and SEs
output <- effectplot(fake.f2, mname1=mname)
output

### plot a phenotype
# Plot of sex-specific phenotype means,
# note that "sex" must be a phenotype name here
effectplot(fake.f2, mname1="sex", geno1=c("F","M"))
# alternatively:
sex <- pull.pheno(fake.f2, "sex")
effectplot(fake.f2, mname1="Sex", mark1=sex, geno1=c("F","M"))

########################################
# two markers plots
########################################

### plot two markers
# plot of genotype-specific phenotype means for 2 markers
mname1 <- find.marker(fake.f2, 1, 37) # marker D1M437
mname2 <- find.marker(fake.f2, 13, 24) # marker D13M254
effectplot(fake.f2, mname1=mname1, mname2=mname2)

### plot two pseudomarkers
#####  refer to pseudomarkers by their positions
effectplot(fake.f2, mname1="1@35", mname2="13@25")

#####  alternatively, find their names via find.pseudomarker
pmnames <- find.pseudomarker(fake.f2, chr=c(1, 13), c(35, 25))
effectplot(fake.f2, mname1=pmnames[1], mname2=pmnames[2])

### Plot of sex- and genotype-specific phenotype means
mname <- find.marker(fake.f2, 13, 24) # marker D13M254
# sex and a marker
effectplot(fake.f2, mname1=mname, mname2="Sex",
           mark2=sex, geno2=c("F","M"))

# Same as above, switch role of sex and the marker
# sex and marker
effectplot(fake.f2, mname1="Sex", mark1=sex,
           geno1=c("F","M"), mname2=mname)

# X chromosome marker
mname <- find.marker(fake.f2, "X", 14) # marker DXM66
effectplot(fake.f2, mname1=mname)

# Two markers, including one on the X
mnames <- find.marker(fake.f2, c(13, "X"), c(24, 14))
effectplot(fake.f2, mname1=mnames[1], mname2=mnames[2])
}

\author{Hao Wu; Karl W Broman, \email{broman@wisc.edu} }

\seealso{ \code{\link{plotPXG}}, \code{\link{find.marker}},
  \code{\link{effectscan}}, \code{\link{find.pseudomarker}} }

\keyword{hplot}
