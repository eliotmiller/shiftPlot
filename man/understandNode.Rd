% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/understandNode.R
\name{understandNode}
\alias{understandNode}
\title{Investigate a node}
\usage{
understandNode(tree, node)
}
\arguments{
\item{tree}{A phylogeny in ape format.}

\item{node}{A node number in the tree.}
}
\value{
Nothing. Results are printed to screen.
}
\description{
Print out which taxa do and do not descend from a given node.
}
\details{
This is a very simple function, with no accounting made for non-bifurcating
trees or for if the user inputs a node without an ancestor in the tree.

#' @author Eliot Miller
}
\examples{
data(phy)
understandNode(phy, 102)
}
