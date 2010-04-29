\name{mqmtestnormal}
\alias{mqmtestnormal}

\title{ Shapiro normality test used for MQM }

\description{
  Wraps a shapiro's normality test from the nortest package. 
  This function is used in MQM to test the normality of the
  trait under investigation
}

\usage{
  mqmtestnormal(cross, pheno.col = 1,significance=0.05, verbose=FALSE)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{pheno.col}{ 
Column number in the phenotype matrix which should be used as the phenotype. This can be a vector of integers. % \phenocol
  }
  \item{significance}{ 
    Significance level used in the normality test. Lower significance levels will accept larger deviations from normality.
  }
  \item{verbose}{
    If TRUE, print result as well as return it.
  }
}

\value{
	Boolean indicating normality of the trait in pheno.col.  (FALSE
	indicates that the phenotype appears to not be normally distributed.)
}

\details{
        For augmented data (as from \code{\link{mqmaugment}}), the cross
        is first reduced to distinct individuals. Furthermore the shapiro used to 
        test normality works only for 3 <= nind(cross) <= 5000
}

\author{ 
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{
	\itemize{
    \item \code{\link[stats]{shapiro.test}} - Function wrapped by our mqmtestnormal
% \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item MQM-tutorial - The online MQM tutorial
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Find the most informative cofactors through backward elimination
    \item \code{\link{mqmsetcofactors}} - Set cofactors
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
data(multitrait)

# test normality of 7th phenotype
mqmtestnormal(multitrait, pheno.col=7)

# take log
multitrait <- transformPheno(multitrait, pheno.col=7, transf=log)
  
# test again
mqmtestnormal(multitrait, pheno.col=7)
}

\keyword{ utilities }
