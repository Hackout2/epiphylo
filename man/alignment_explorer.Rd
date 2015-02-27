\name{alignment_explorer}
\alias{alignment_explorer}
\title{Extracts and plots information on DNA alignments}
\description{
  This function extracts information from the DNA sequences in an
  \linkS4class{obkData} object and prints or plots them.
}
\usage{
alignment_explorer(x)
}
\arguments{
  \item{x}{an \linkS4class{obkData} object.}
}
\details{
  The function scans the data for the different genes and plots the output successively (the user is asked to type enter).
}
\author{Emmanuel Paradis}
\value{
  nothing, the results are printed on the current graphical device.
}
\examples{
\dontrun{
data(ToyOutbreak)
alignment_explorer(ToyOutbreak)
}}
\keyword{hplot}
