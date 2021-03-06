% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCoords.R
\name{getCoords}
\alias{getCoords}
\title{Get coordinates for plotting tree}
\usage{
getCoords(tree)
}
\arguments{
\item{tree}{A phylogeny in ape format. This should be the polytomy tree, i.e. the
result of a call to polytomyBind.}
}
\value{
A vector providing the time of the first branching event for each clade in the
tree that was collapsed.
}
\description{
Should be an unexported utility function but I'm too lazy.
}
\author{
Eliot Miller and Bruce Martin
}
