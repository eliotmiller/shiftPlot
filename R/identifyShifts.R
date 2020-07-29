#' Identify shift points along the phylogeny in trait state
#'
#' Identify gains and losses in the trait by the node that subtends the branch a shift
#' is inferred to have occurred on.
#'
#' @param tree Phylogeny in ape format, corresponding to states.df.
#' @param states.df Data frame in the specified shiftPlot format. Should contain one
#' column named "present", and a 0 or a 1 for every node in the phylogeny, with the tips
#' above the internal nodes, and no row names. See details and examples. 
#' @param flip.tips Whether or not to flip tips such that they have the values of the
#' parent node. This will preclude shifts being detected on branches leading to tips. 
#' 
#' @details states.df should have one column titled "present". This column should be
#' coded either as a 0 or a 1, indicating the presence of the trait. states.df should have
#' as many rows as there are nodes in phylogeny, and the tip nodes should come first in
#' the data frame. For example, you might rbind the $tip.states and $states objects from
#' a corHMM output together to create states.df
#'
#' @return A list of lists, with gains and losses of the trait summarized by the node
#' subtending the branch on which a shift was inferred to have occurred. 
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
#' nodeStates <- data.frame(present=Precur_res.corHMM$states[,3]+Precur_res.corHMM$states[,4])
#' tipStates <- data.frame(present=Precur_res.corHMM$tip.states[,3]+Precur_res.corHMM$tip.states[,4])
#'
#' #note that tip states comes first here!
#' states <- rbind(tipStates, nodeStates)
#' 
#' #binarize this. choosing to call 0.5 chance of having trait present
#' states$present[states$present >= 0.5] <- 1
#' states$present[states$present < 0.5] <- 0
#' 
#' #flip node 103 and all nodes towards tips from there to trait = absent
#' induced <- states
#' induced[geiger:::.get.descendants.of.node(103, phy),] <- 0
#' induced[103,] <- 0
#'
#' #get rid of row names
#' row.names(induced) <- NULL
#' 
#' #try running the function
#' result <- identifyShifts(phy, induced, FALSE)
#' result
#' tipsInQ <- tips(phy, 103)
#' states[tipsInQ,] <- 0
#' row.names(states) <- NULL
#' result <- identifyShifts(phy, states, TRUE)
#' result

identifyShifts <- function(tree, states.df, flip.tips)
{
	if(flip.tips==TRUE)
	{
		states.df <- tipFlip(tree, states.df)
	}

	#create a node column for simplicity sake
	states.df$node <- row.names(states.df)

	#create a shifts object
	shifts <- data.frame(tree$edge)
	names(shifts) <- c("parent.node","daughter.node")
	
	#merge in the state data with the parent node and rename for ease
	shifts2 <- merge(shifts, states.df, by.x="parent.node", by.y="node")
	names(shifts2)[3] <- "parent.state"
	
	#same but with the daughter node
	shifts3 <- merge(shifts2, states.df, by.x="daughter.node", by.y="node")
	names(shifts3)[4] <- "daughter.state"
	
	#the merge re-arranges the column orders. put back here so it looks prettier
	shifts3 <- shifts3[,c("parent.node","daughter.node","parent.state","daughter.state")]

	#now if the parent state goes from 0 to daughter state 1, it's a gain, and vice versa
	shifts3$gain <- 0
	shifts3$gain[shifts3$parent.state==0 & shifts3$daughter.state==1] <- 1
	shifts3$loss <- 0
	shifts3$loss[shifts3$parent.state==1 & shifts3$daughter.state==0] <- 1

	#create a list to summarize the descendents from each node subtending each gain and loss
	gains <- list()
	tempNodes <- shifts3$daughter.node[shifts3$gain==1]

	#if there aren't any gains the below loop will throw an error, so deal with that
	if(length(tempNodes) == 0)
	{
		gains <- "No gains were found in this reconstruction"
	}

	else
	{
		for(i in 1:sum(shifts3$gain==1))
		{	
			#you want the daughter node that 
			gains[[i]] <- geiger::tips(tree, tempNodes[i])
		}
	}

	#name the gains with the node number
	names(gains) <- tempNodes

	#create a list to summarize the descendents from each node subtending each gain and loss
	losses <- list()
	tempNodes <- shifts3$daughter.node[shifts3$loss==1]

	#if there aren't any losses the below loop will throw an error, so deal with that
	if(length(tempNodes) == 0)
	{
		losses <- "No losses were found in this reconstruction"
	}

	else
	{
		for(i in 1:sum(shifts3$loss==1))
		{	
			#you want the daughter node that 
			losses[[i]] <- geiger::tips(tree, tempNodes[i])
		}
	}

	#name the gains with the node number
	names(losses) <- tempNodes

	#bundle these and return
	results <- list(gains, losses)
	names(results) <- c("gains","losses")

	results
}
