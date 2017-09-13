\name{mqmplot.circle}
\alias{mqmplot.circle}
\title{Circular genome plot for MQM}

\description{
  Circular genome plot - shows QTL locations and relations.
}

\usage{
mqmplot.circle(cross,result,highlight=0,spacing=25, interactstrength=2,
               axis.legend=TRUE, col.legend=FALSE, verbose=FALSE, transparency=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross} with optionally phenotype locations. See
    \code{\link{read.cross}} for details on reading in cross objects, and optionally
    \code{\link{addloctocross}} for adding phenotype locations.
  }
  \item{result}{An object of class \code{mqmmulti} or \code{scanone}.
    See \code{\link{mqmscanall}} \code{\link{scanone}} for details.
  }
  \item{highlight}{ With a mqmmulti object, highlight this phenotype (value between
    one and the number of results in the mqmmultiobject)
  }
  \item{interactstrength}{ When highlighting a trait, consider interactions
  significant they have a change of more than interactstrength*SEs. A higher value
  will show less interactions.  However the interactions reported at higher
  interactstrength values will generaty be more reliable.
  }
  \item{spacing}{ User defined spacing between chromosomes in cM }
  \item{axis.legend}{ When set to FALSE, suppresses the legends. (defaults to plotting legends besides the axis.  }
  \item{col.legend}{ With a mqmmulti object, plots a legend for the non-highlighed version }
  \item{transparency}{ Use transparency when drawing the plots (defaults to no transparency) }
  \item{verbose}{ Be verbose }
}

\details{
  Depending on the input of the result being either \code{\link{scanone}} or \code{mqmmulti} a different plot is drawn.
  If model information is present from \code{\link{mqmscan}} (by setting cofactors) This will be highlighted in
  red (see example).
  If phenotypes have genetic locations (e.g. eQTL) they will be plotted on the genome otherwise
  phenotypes will be plotted in the middle of the circle (with a small offset)
  Locations can be added by using the \code{\link{addloctocross}} function.
}

\value{
    Plotting routine, no return
}

\author{
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{
    \itemize{
    % \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{https://rqtl.org/tutorials/MQM-tour.pdf}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/doc/Sources/MQM/mqm/standard_seealso.txt ^^-----
    }
}

\examples{
data(multitrait)
\dontshow{multitrait$pheno <- multitrait$pheno[,1:3]}
data(locations)
\dontshow{locations <- locations[1:3,]}

multifilled <- fill.geno(multitrait)                    # impute missing genotypes
multicof <- mqmsetcofactors(multitrait,10)              # create cofactors
multiloc <- addloctocross(multifilled,locations)        # add phenotype information to cross
multires <- mqmscanall(multifilled,cofactors=multicof)  # run mqmscan for all phenotypes

#Basic mqmmulti, color = trait, round circle = significant
mqmplot.circle(multifilled,multires)

#mqmmulti with locations of traits in multiloc
mqmplot.circle(multiloc,multires)

#mqmmulti with highlighting
mqmplot.circle(multitrait,multires,highlight=3)

#mqmmulti with locations of traits in multiloc and highlighting
mqmplot.circle(multiloc,multires,highlight=3)
}

\keyword{ utilities }
