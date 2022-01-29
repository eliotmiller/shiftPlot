#' Determine which clades can be collapsed
#'
#' Optimal collapse
#'
#' @param orig.tree Phylogeny in ape format, corresponding to states.df.
#' @param states.df Data frame in the specified shiftPlot format. Should contain one
#' column named "state", and number for every node in the phylogeny, with the tips
#' above the internal nodes, and no row names. See details and examples. 
#' @param flip.tips Whether or not to flip tips to the state of their parent nodes.
#' 
#' @details states.df should have one column titled "present". This column should be
#' coded either as a 0 or a 1, indicating the presence of the trait. states.df should have
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
#' @import ape corHMM geiger
#' @importFrom graphics plot points polygon segments text
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
#' #flip node 103 and all nodes towards tips from there to trait = absent
#' induced <- states
#' induced[geiger:::.get.descendants.of.node(103, phy),] <- 0
#' induced[103,] <- 0
#'
#' #get rid of row names
#' row.names(induced) <- NULL
#' 
#' #run the function and don't flip those tips
#' ocResult <- optimalCollapse(orig.tree=phy, states.df=induced, flip.tips=FALSE)

optimalCollapse <- function(orig.tree, states.df, flip.tips)
{
  #add some checks to make sure input looks about right
  if(dim(states.df)[1] != length(orig.tree$tip.label) & dim(states.df)[2] != 1)
  {
    stop("Your input doesn't look right")
  }
  
	#find the root node
	rootNode <- length(orig.tree$tip.label)+1

	if(flip.tips==TRUE)
	{
		states.df <- tipFlip(orig.tree, states.df)
	}
	
	#now create a node column for simplicity sake
	states.df$node <- row.names(states.df)

	#and a species col for simplicity sake (leave it blank for internal nodes)
	states.df$species <- c(orig.tree$tip.label, rep("", dim(states.df)[1]-length(orig.tree$tip.label)))

	#also set up a clade column in the states table
	states.df$clade <- ""

	#create a table like the identifyShifts fxn, but use differently
	shifts <- data.frame(orig.tree$edge)
	names(shifts) <- c("parent.node","daughter.node")

	#start at the first tip, bump down to the parent node, figure out all the taxa that
	#descend from that node, and figure out their states. if they are all the same, bump
	#one node further towards the root. if the states aren't all the same, go back up and
	#log that parent node as one you can collapse.
	for(i in 1:length(orig.tree$tip.label))
	{
		#this way as you log what taxa you can collapse, you can skip ahead to the next
		#clade you'll need to deal with (and not have to solve for all the other taxa in
		#a clade you just discovered you can collapse)
		if(states.df$clade[i] != "")
		{
			next()
		}

		#set aside this vector of states (it'll be only one state at first) for use in the
		#while loop
		allStates <- states.df$state[i]

		#log which node you are considering (it's a tip at first but it'll update below)
		node <- i

		#set up a placeholder
		placeHolder <- rootNode + 1

		#this while loop will run while all the tips under consideration have the same state
		#and while the node number is > the root node, so that you don't error out. have to use
		#this placeholder instead of the node itself, because on the first run of the while loop 
		#node will be equal to a tip, which is greater than the root node.
		while(length(unique(allStates))==1 & placeHolder > rootNode)
		{
			#log the parent node of the node in question
			parentNode <- shifts$parent.node[shifts$daughter.node==node]

			#find all the tips that descend from the parent node
			allTips <- geiger::tips(orig.tree, parentNode)

			#pull the state data for all the tips that descend from parentNode
			allStates <- states.df$state[states.df$species %in% allTips]

			#if allStates are the same, and if parent node is greater than the root node,
			#then set node and the placeholder to parent node, because we'll be starting
			#this while loop over again and walking down the tree in this way
			if(length(unique(allStates))==1 & parentNode > rootNode)
			{
				node <- parentNode
				placeHolder <- parentNode
			}
		}

		#if the while loop worked, then all tips descending from node should be able to be collapsed.
		allTips <- geiger::tips(orig.tree, node)
		states.df$clade[states.df$species %in% allTips] <- paste("node", node, sep="_")
	}

	#prep the object for export
	temp <- unique(states.df$clade[states.df$clade != ""])
	tempNodes <- as.numeric(unlist(lapply(strsplit(temp, "_"), "[", 2)))
	prepped <- data.frame(node=tempNodes, state=0, collapse=1, clade=temp, stringsAsFactors=FALSE)

	#pull the state of each of these clades
	for(i in 1:length(temp))
	{
	  prepped$state[i] <- states.df$state[states.df$clade == temp[i]][1]
	}

	prepped
}
