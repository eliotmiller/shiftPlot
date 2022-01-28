#' Set up a table for node-matching 
#'
#' Set up a table for node-matching between two different trees.
#'
#' @param tree1 The first tree.
#' @param tree2 The second tree.
#' 
#' @details Some things.
#'
#' @author Eliot Miller
#' 
#' @export

tableSetup <- function(tree1, tree2)
{
  output <- data.frame(matrix(nrow=min(c(length(tree1$tip.label),length(tree2$tip.label)))-1, ncol=2))
  names(output) <- c("tree1", "tree2")
  output
}
