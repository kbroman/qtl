\name{write.cross}
\alias{write.cross}

\title{Write data for a QTL experiment to a file}

\description{
  Data for a QTL experiment is written to a file (or files).
}

\usage{
write.cross(cross, format=c("csv", "csvr", "csvs", "csvsr",
                            "mm", "qtlcart", "gary", "qtab"),
            filestem="data", chr, digits=NULL, descr)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
    \code{\link{read.cross}} for details.}
  \item{format}{Specifies whether to write the data in comma-delimited,
    rotated comma-delimited, Mapmaker, QTL Cartographer, Gary
    Churchill's, or QTAB format.} 
  \item{filestem}{A character string giving the first part of the output
    file names (the bit before the dot).  In Windows, use forward
    slashes (\code{"/"}) or double backslashes (\code{"\\\\"}) to
    specify directory trees.} 
  \item{chr}{A vector specifying for which chromosomes genotype data
    should be written.  This should be a vector of character strings
  referring to chromosomes by name; numeric values are converted to
  strings.  Refer to chromosomes with a preceding \code{-} to have all
  chromosomes but those considered.  A logical (TRUE/FALSE) vector may
  also be used.} 
  \item{digits}{Number of digits to which phenotype values and genetic
     map positions should be rounded. If NULL (the default), they are
     not rounded.}
  \item{descr}{Character string description; used only with \code{format="qtab"}.}
}

%\value{
%  Data is written to a single file (comma-delimited format) or two files
%  (Mapmaker format).
%}

\details{
  Comma-delimited formats: a single csv file is created in the formats
  \code{"csv"} or \code{"csvr"}. Two files are created (one for the
  genotype data and one for the phenotype data) for the formats
  \code{"csvs"} and \code{"csvsr"}; if \code{filestem="file"}, the two
  files will be names \code{"file_gen.csv"} and \code{"file_phe.csv"}.
  See the help file for \code{\link{read.cross}} for details on these formats.

  Mapmaker format: Data is written to two files.
  Suppose \code{filestem="file"}.  Then \code{"file.raw"} will contain
  the genotype and phenotype data, and \code{"file.prep"} will contain
  the necessary code for defining the chromosome assignments, marker
  order, and inter-marker distances.

  QTL Cartographer format: Data is written to two files.    Suppose
  \code{filestem="file"}.  Then \code{"file.cro"} will contain 
  the genotype and phenotype data, and \code{"file.map"} will contain
  the genetic map information.  Note that cross types are converted to
  QTL Cartographer cross types as follows: riself to RF1, risib to RF2,
  bc to B1 and f2 to RF2. 

  Gary's format: Data is written to six files. They are: \cr
  \code{"geno.data"} - genotype data; \cr
  \code{"pheno.data"} - phenotype data; \cr
  \code{"chrid.dat"} - the chromosome identifier for each marker; \cr
  \code{"mnames.txt"} - the marker names; \cr
  \code{"markerpos.txt"} - the marker positions; \cr
  \code{"pnames.txt"} - the phenotype names

  QTAB format: See
  \href{https://github.com/qtlHD/qtlHD/blob/master/doc/input/qtab.md}{documentation}.

}


\examples{
\dontrun{data(fake.bc)

# comma-delimited format
write.cross(fake.bc, "csv", "Data/fakebc", c(1,5,13))

# rotated comma-delimited format
write.cross(fake.bc, "csvr", "Data/fakebc", c(1,5,13))

# split comma-delimited format
write.cross(fake.bc, "csvs", "Data/fakebc", c(1,5,13))

# split and rotated comma-delimited format
write.cross(fake.bc, "csvsr", "Data/fakebc", c(1,5,13))

# Mapmaker format
write.cross(fake.bc, "mm", "Data/fakebc", c(1,5,13))

# QTL Cartographer format
write.cross(fake.bc, "qtlcart", "Data/fakebc", c(1,5,13))

# Gary's format
write.cross(fake.bc, "gary", c(1,5,13))} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}; Hao Wu; Brian
  S. Yandell; Danny Arends}

\seealso{ \code{\link{read.cross}} }

\keyword{IO}
