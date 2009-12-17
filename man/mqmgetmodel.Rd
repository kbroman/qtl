\name{mqmgetmodel}
\alias{mqmgetmodel}
\alias{mqmgetmodel}

\title{ Retrieves selected cofactors from the results of an mqm scan }

\description{
  Retrieves selected cofactors from the results of an mqm scan
}

\usage{
	mqmgetmodel(scanresult)
}

\arguments{
  \item{scanresult}{An object of class \code{scanone} created bij a mqm scan See \code{\link[qtl]{mqm}} for details.} 
}

%\details{
%
%}

\value{
	The model created (cofactors selected in modeling, and used during scanning) during the mqmscan
}  

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl}}

\note{
If you find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}

\seealso{
	\itemize{
		\item \code{\link[qtl]{mqm}} - Genome scan with a multiple QTL model (MQM)
    \item \code{\link[qtl]{mqmsetcofactors}} - Genome scan with a multiple QTL model (MQM)
	}
}

\keyword{ utilities }
