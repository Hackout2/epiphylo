\name{tree_time}
\alias{tree_time}
\alias{tree_time_1}
\alias{tree_time_2}
\title{Plot Phylogenies and Time}
\description{
  These two functions plot a tree built from the sequences availabl in
  the \linkS4class{obkData} object together with the time of sampling of
  each sequence (also taken from the same object).
}
\usage{
tree_time_1(x)
tree_time_2(x)
}
\arguments{
  \item{x}{an \linkS4class{obkData} object.}
}
\details{
  The first function \code{tree_time_1} plots the NJ tree growing
  upwards, without the branch length information, and the time axis is
  represented on top with from left to right. The second function
  \code{tree_time_2} plots the NJ tree in the usual way, from left to
  right, with the branch length information, and the time axis is
  represented on bottom also from left to right.

  Coloured lines show the link between the tips of the tree (i.e.,
  unique sequences) and the time of sampling of the individuals. The
  colour scale goes from blue (oldest dates) to red (youngest dates).

  The functions scan the data for the different genes and plots the
  output successively (the user is asked to type enter).
}
\note{
  The neighbor-joigning (NJ) method reconstructs unrooted trees
  whereas the present functions represent them as rooted.
}
\value{nothing the results are printed on the current graphical device.}
\author{Emmanuel Paradis}
\examples{
\dontrun{
data(ToyOutbreak)
tree_time_1(ToyOutbreak)
tree_time_2(ToyOutbreak)
}
}
\keyword{hplot}
