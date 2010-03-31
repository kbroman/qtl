\name{mqmplot_circle}
\alias{mqmplot_circle}
\alias{mqmplotcircle}

\title{Circular genome plotting - routine to show QTL information}

\description{
  Circular genome plotting - routine to show QTL information in a circular fashion.
}

\usage{
	mqmplot_circle(cross,result,highlight=0,spacing=25, interactstrength=2, legend=FALSE,verbose=FALSE,transparency=FALSE)
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
  \item{legend}{ With a mqmmulti object, plots a legend for the non-highlighed version }
  \item{transparency}{ Use transparency when drawing the plots (defaults to no transparency) }
  \item{verbose}{ Be verbose }
}

\details{
  Depending on the input of the result being either \code{\link{scanone}} or \code{mqmmulti} a different plot is drawn.
  If model information is present from \code{\link{mqm}} (by setting cofactors) This will be highlighted in 
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
    % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmsetcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmpermutation}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
	}
}

\examples{
  \dontrun{
    library(qtl)
    data(hyper)
    data(multitrait)
    data(locations)
    
    hyperfilled <- fill.geno(hyper)
    hypercof <- mqmsetcofactors(hyper,4)
    hyperres <- mqmscan(hyperfilled,hypercof)
    hyperscanone <- scanone(hyperfilled)

    multifilled <- fill.geno(multitrait)
    multicof <- mqmsetcofactors(multitrait,10)
    multiloc <- addloctocross(multifilled,locations)
    multires <- mqmscanall(multifilled,cofactors=multicof)
    
    #Basic scanone, thickness of markers = LODscore
    mqmplot_circle(hyperfilled,hyperscanone)
    Sys.sleep(3)
    #Basic mqm, thickness of markers = LODscore, red = significant
    mqmplot_circle(hyperfilled,hyperres)
    Sys.sleep(3)
    #Basic mqmmulti, color = trait, round circle = significant
    mqmplot_circle(multifilled,multires)
    Sys.sleep(3)
    #mqmmulti with locations of traits in multiloc
    mqmplot_circle(multiloc,multires)
    Sys.sleep(3)
    #mqmmulti with highlighting
    mqmplot_circle(multitrait,multires,highlight=16)
    Sys.sleep(3)
    #mqmmulti with locations of traits in multiloc and highlighting
    mqmplot_circle(multiloc,multires,highlight=16)
  }
}

\keyword{ utilities }
