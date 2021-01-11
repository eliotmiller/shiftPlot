#' Investigate a node
#'
#' Print out which taxa do and do not descend from a given node.
#'
#' @param tree A phylogeny in ape format.
#' @param node A node number in the tree.
#' 
#' @details This is a very simple function, with no accounting made for non-bifurcating
#' trees or for if the user inputs a node without an ancestor in the tree.
#' 
#' #' @author Eliot Miller
#'
#' @return Nothing. Results are printed to screen.
#'
#' @export
#' 

understandNode <- function(tree, node)
{
	print("This node is the MRCA of this/these taxon/taxa")
	print(tree$tip.label[geiger:::.get.descendants.of.node(node, tree, tips=TRUE)])
	print("The sister clade contains this/these taxon/taxa")
	#if the user inputs the root node this will fail
	anc <- tree$edge[,1][tree$edge[,2]==node]
	#if a node gives rise to more than two other nodes, this will fail
	other <- tree$edge[,2][tree$edge[,1]==anc]
	other <- other[other != node]
	print(tree$tip.label[geiger:::.get.descendants.of.node(other, tree, tips=TRUE)])
}
