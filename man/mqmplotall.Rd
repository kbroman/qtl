\name{mqmplot_multitrait}
\alias{mqmplot_multitrait}

\title{Plotting routine to display the results from a mqm scan on multiple phenotypes}

\description{
  Plotting routine to display the results from a mqm scan on multiple phenotypes. 
  It supports four diffferent visualizations: a contourmap, heatmap, 3D graph or a multiple 
  QTL plot from \code{\link[qtl]{plot.scanone}}
}

\usage{
mqmplot_multitrait(result, type = "C", theta = 30, phi = 15, \dots)
}

\arguments{
  \item{result}{ Result object from scanall }
  \item{type}{ "C" for contourplot, "P" for stacked R/QTL plot, "I" for a heatmapplot, "D" for a 3D plot  }
  \item{theta}{ Horizontal axis rotation in a 3D plot}
  \item{phi}{ Vertical axis rotation in a 3D plot}
  \item{\dots}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\seealso{ 
  \itemize{
  	\item \code{\link[qtl]{mqmall}} - Function called to analyse multiple traits in the cross object
  	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for MQM
  }
}

\examples{
data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
multitrait$pheno <- multitrait$pheno[,1:3]}
multitrait <- fill.geno(multitrait) # impute missing genotype data
result <- mqmall(multitrait, doLOG=1)
mqmplot_multitrait(result,"P")
mqmplot_multitrait(result,"C")
mqmplot_multitrait(result,"I")
mqmplot_multitrait(result,"D")
}

\references{
\itemize{
	\item Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988), \emph{The New S Language.} Wadsworth & Brooks/Cole.
}
}

\keyword{ hplot }
