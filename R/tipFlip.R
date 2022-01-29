#' Flip tips to have the states of their parent nodes
#'
#' Ensure tips have the same states as their parent nodes, to preclude shifts on terminal
#' branches.
#'
#' @param orig.tree Phylogeny in ape format, corresponding to states.df.
#' @param states.df Data frame in the specified shiftPlot format. Should contain one
#' column named "state", and a number (state) for every node in the phylogeny, with the tips
#' above the internal nodes, and no row names. See details and examples. 
#' 
#' @details states.df should have one column titled "state". This column should take the form
#' of a numeric vector, indicating the statee of the trait. states.df should have
#' as many rows as there are nodes in phylogeny, and the tip nodes should come first in
#' the data frame. For example, you might rbind the $tip.states and $states objects from
#' a corHMM output together to create states.df The main reason one would use this
#' function is to preclude shifts on terminal branches. You might want to do that if, for
#' example, your phylogeny was so large that corHMM was having trouble calculating the
#' marginal probabilities of a tip having a hidden trait, and the result would otherwise
#' be that it appeared there were many shifts in the hidden trait towards the tips.
#'
#' @return An updated states.df with the tip values flipped as needed.
#' 
#' @author Eliot Miller
#'
#' @export
#' 
#' @examples
#' #start with a corHMM output and build up a states.df
#' #load data. these are the results of following the corHMM precursor model vignette
#' data(Precur_res.corHMM)
#' data(phy)
#' nodeStates <- data.frame(state=Precur_res.corHMM$states[,3]+Precur_res.corHMM$states[,4])
#' tipStates <- data.frame(state=Precur_res.corHMM$tip.states[,3]+Precur_res.corHMM$tip.states[,4])
#'
#' #note that tip states comes first here!
#' states <- rbind(tipStates, nodeStates)
#' 
#' #binarize this. choosing to call 0.5 chance of having trait present
#' states$state[states$state >= 0.5] <- 1
#' states$state[states$state < 0.5] <- 0
#' 
#' #make the model a little more exciting (all tips were origianlly
#' #inferred to have the hidden precursor trait)
#' tipsInQ <- tips(phy, 103)
#' states[tipsInQ,] <- 0
#'
#' #' #get rid of row names
#' row.names(states) <- NULL
#' 
#' #try running the function
#' result <- tipFlip(orig.tree=phy, states.df=states)

tipFlip <- function(orig.tree, states.df)
{
	#create a node column for simplicity sake
	states.df$node <- row.names(states.df)

	#create a shifts object
	shifts <- data.frame(orig.tree$edge)
	names(shifts) <- c("parent.node","daughter.node")
	
	#merge in the state data with the parent node and rename for ease
	shifts2 <- merge(shifts, states.df, by.x="parent.node", by.y="node")
	names(shifts2)[3] <- "parent.state"
	
	#same but with the daughter node
	shifts3 <- merge(shifts2, states.df, by.x="daughter.node", by.y="node")
	names(shifts3)[4] <- "daughter.state"
	
	#the merge re-arranges the column orders. put back here so it looks prettier
	shifts3 <- shifts3[,c("parent.node","daughter.node","parent.state","daughter.state")]

	#set all the daughter.states of the tips to those of their parents
	tipNos <- 1:length(orig.tree$tip.label)
	shifts3$daughter.state[shifts3$daughter.node %in% tipNos] <- shifts3$parent.state[shifts3$daughter.node %in% tipNos]

	#recreate the states.df in the binarized form and return
	results <- shifts3[,c("daughter.node","daughter.state")]

	#don't forget to include the root node and state
	rootNode <- max(tipNos)+1
	toBind <- data.frame(daughter.node=rootNode,
		daughter.state=unique(shifts3[shifts3$parent.node==rootNode,"parent.state"]))
	results <- rbind(results, toBind)
	results <- results[order(results$daughter.node),]
	row.names(results) <- NULL
	names(results) <- c("delete","state")
	results$delete <- NULL
	results
}
