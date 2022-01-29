#' Get coordinates for plotting tree
#'
#' Should be an unexported utility function but I'm too lazy.
#'
#' @param tree A phylogeny in ape format.
#' 
#' @author Eliot Miller and Bruce Martin
#'
#' @return A list providing details on the plotted tree structure.
#'
#' @export
#' 

getCoords <- function(tree)
{
	plot(tree, show.tip.label=FALSE, plot=FALSE)
	coords <- get("last_plot.phylo", envir=.PlotPhyloEnv)
	coords
}
