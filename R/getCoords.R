#' Get coordinates for plotting tree
#'
#' Should be an unexported utility function but I'm too lazy.
#'
#' @param tree A phylogeny in ape format. This should be the polytomy tree, i.e. the
#' result of a call to polytomyBind.
#' 
#' @author Eliot Miller and Bruce Martin
#'
#' @return A vector providing the time of the first branching event for each clade in the
#' tree that was collapsed.
#'
#' @export
#' 

getCoords <- function(tree)
{
	plot(tree, show.tip.label=FALSE)
	coords <- get("last_plot.phylo", envir=.PlotPhyloEnv)
	coords
}
