\name{tree_obkData}
\alias{tree_obkData}
\title{Plot phylogenies with frequencies of unique sequences}
\description{
  This function plots a phylogeny built from the unique sequences
  together with a representation of the abundance (frequencies) of
  these in the data.
}
\usage{
tree_obkData(x, ncatmax = 10, colfun = topo.colors)
}
\arguments{
  \item{x}{an \linkS4class{obkData} object.}
  \item{ncatmax}{an integer giving the number of categories of
    haplotype frequency.}
  \item{colfun}{a function to define the colours.}
}
\details{
  First, the unique sequences are extracted using the function
  \code{haplotype} (in \pkg{pegas}). Then, a neighbor-joigning (NJ) tree
  is built. The tree is plotted on the left-hand side of the graph, and
  the frequencies of each unique sequence is represented in two ways:
  coloured labels at the tips of the tree (with the colour scale drawn
  at the top), and a horizontal barplot on the right-hand side of the
  graph.

  The function scans the data for the different genes and plots the
  output successively (the user is asked to type enter)
}
\value{nothing the results are printed on the current graphical device.}
\author{Emmanuel Paradis}
\examples{
\dontrun{
data(ToyOutbreak)
tree_obkData(ToyOutbreak)
}
}
\keyword{hplot}
