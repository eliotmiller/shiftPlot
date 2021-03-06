% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiftPlotter.R
\name{shiftPlotter}
\alias{shiftPlotter}
\title{Plot a phylogeny with triangles representing clades}
\usage{
shiftPlotter(
  orig.tree,
  poly.tree,
  coords,
  identifyShifts.obj,
  presence.color,
  absence.color,
  pt.cex
)
}
\arguments{
\item{orig.tree}{A phylogeny in ape format. This should be the ORIGINAL tree, not the
collapsed tree nor the polytomy tree.}

\item{poly.tree}{A phylogeny in ape format. This should be the polytomy tree, i.e. the
result of a call to polytomyBind.}

\item{coords}{The result of a call to getCoords. These are returned when the trianglePlotter
function is used, which allows a user to save those and re-use them without overwriting the
triangle plot by calling getCoords again.}

\item{identifyShifts.obj}{Result of a call to identifyShifts.}

\item{presence.color}{What color clades for which the trait is present should be colored.}

\item{absence.color}{What color clades for which the trait is absent should be colored.}

\item{pt.cex}{The size of the shift points.}
}
\description{
Use the optimal collapses to create a tree where triangles represent clades.
}
\details{
Stuff.
}
\examples{
#start with a corHMM output and build up a states.df
#load data. these are the results of following the corHMM precursor model vignette
data(Precur_res.corHMM)
data(phy)
nodeStates <- data.frame(present=Precur_res.corHMM$states[,3]+Precur_res.corHMM$states[,4])
tipStates <- data.frame(present=Precur_res.corHMM$tip.states[,3]+Precur_res.corHMM$tip.states[,4])

#note that tip states comes first here!
states <- rbind(tipStates, nodeStates)

#binarize this. choosing to call 0.5 chance of having trait present
states$present[states$present >= 0.5] <- 1
states$present[states$present < 0.5] <- 0

#flip node 103 and all nodes towards tips from there to trait = absent
induced <- states
induced[geiger:::.get.descendants.of.node(103, phy),] <- 0
induced[103,] <- 0

#get rid of row names
row.names(induced) <- NULL

#find the optimal collapse configuration
collapsed <- optimalCollapse(phy, induced, FALSE)

shifts <- identifyShifts(phy, induced, FALSE)

#collapse the tree
dropped <- dropManyTips(phy, collapsed)

#use the polytomyBind function
polyTree <- polytomyBind(dropped, collapsed)

branchingResults <- firstBranches(phy, dropped, collapsed)

saveMe <- trianglePlotter(tree=polyTree, dropped.results=dropped, clade.table=collapsed,
branches=branchingResults, identifyShifts.obj=shifts,
presence.color="red", absence.color="black",
label.offset=0.3, text.cex=0.09,
root.state="present")

shiftPlotter(phy, polyTree, saveMe, shifts, "red", "black", 1)
}
\author{
Eliot Miller
}
