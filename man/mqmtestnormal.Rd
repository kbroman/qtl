\name{mqmtestnormal}
\alias{mqmtestnormal}

\title{ Pearsons normality test used for MQM }

\description{
  Wraps a pearsons normality test from the nortest package. 
  This function is used in MQM to test the normality of the
  trait under investigation
}

\usage{
mqmtestnormal(cross, pheno.col = 1)
}

\arguments{
  \item{cross}{
    % \crossobject
  \item{pheno.col}{ 
    % \phenocol
  }
}

\value{
	Boolean indicating normality of the trait in pheno.col
}

\author{ 
  Danny Arends \email{Danny.Arends@gmail.com} % \dannyauthor
}

\seealso{
	\itemize{
% \input{"inst/docs/mqm/standard_seealso.txt"}
  }
}

\keyword{ utilities }
