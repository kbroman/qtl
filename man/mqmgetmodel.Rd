\name{mqmgetmodel}
\alias{mqmgetmodel}
\alias{mqmgetmodel}

\title{ Retrieve the QTL model used in mapping from the results of an mqm scan }

\description{
  Retrieves the QTL model used for scanning from the output of an mqm scan. The model 
  only contains the selected cofactors significant at the specified cofactor.significance 
  from the results of an mqm scan
}

\usage{
	mqmgetmodel(scanresult)
}

\arguments{
  \item{scanresult}{An object of class \code{scanone} 
    created bij a mqm scan See \code{\link[qtl]{mqm}} for details.
  } 
}

\details{
  returns NULL when no cofactors reach the significance threshold during QTL modeling.
}

\value{
	The function returns the QTL model created (cofactors selected in modeling found to 
  be significant at cofactor.significance) and used during the mqmscan
  See for the format of the model: \code{\link[qtl]{makeqtl}}
}  

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl}}

\note{
If you find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}

\seealso{
	\itemize{
		\item \code{\link[qtl]{mqm}} - Genome scanning using a multiple QTL model (MQM)
    \item \code{\link[qtl]{mqmsetcofactors}} - Setting multiple cofactors for backward elimination
    \item \code{\link[qtl]{makeqtl}} - Make a qtl object
	}
}

\keyword{ utilities }
