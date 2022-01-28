#' Identify shift points along the phylogeny in trait state
#'
#' Identify shifts in the state of a trait by the node that subtends the branch a shift
#' is inferred to have occurred on.
#'
#' @param orig.tree Phylogeny in ape format, corresponding to states.df.
#' @param states.df Data frame in the specified shiftPlot format. Should contain one
#' column named "state", and number for every node in the phylogeny, with the tips
#' above the internal nodes, and no row names. See details and examples. 
#' @param flip.tips Whether or not to flip tips such that they have the values of the
#' parent node. This will preclude shifts being detected on branches leading to tips. 
#' 
#' @details states.df should have one column titled "state". This column should be
#' take the form of a number, indicating the state of the trait. states.df should have
#' as many rows as there are nodes in phylogeny, and the tip nodes should come first in
#' the data frame. For example, you might rbind the $tip.states and $states objects from
#' a corHMM output together to create states.df. In my experience, there are limited cases
#' where the states of internal nodes do not match intuition after running marginal tip
#' reconstruction in corHMM. The optimalCollapse function seem impermeable to these issues,
#' but this function will detect these counterintuitive shifts in state. The trianglePlotter
#' will warn if it detects any issues where shifts do not match between this function and
#' a collapsed clade from optimalCollapse. If that is the case, users may wish to change
#' the state of the internal nodes in question and re-run this function before proceeding.
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
#' shiftObj <- identifyShifts(orig.tree=phy, states.df=induced, flip.tips=FALSE)

identifyShifts <- function(orig.tree, states.df, flip.tips)
{
  #add some checks to make sure input looks about right
  if(dim(states.df)[1] != length(orig.tree$tip.label) & dim(states.df)[2] != 1)
  {
    stop("Your input doesn't look right")
  }
  
  if(flip.tips==TRUE)
	{
		states.df <- tipFlip(orig.tree, states.df)
	}

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

	#figure out how many unique shift types there are
	uniShifts <- unique(shifts3[,c("parent.state","daughter.state")])
	uniShifts <- uniShifts[uniShifts$parent.state!=uniShifts$daughter.state,]
	
	#go through the whole thing and tabulate which shifts occur where
	for(i in 1:dim(uniShifts)[1])
	{
	  #create a temporary column
	  shifts3$temp <- 0
	  
	  #subset this temp column to instances where daughter and parent states match those
	  #in row i of uniShifts
	  shifts3$temp[shifts3$parent.state==uniShifts$parent.state[i] & shifts3$daughter.state==uniShifts$daughter.state[i]] <- 1
	  
	  #change the name of the temp column to the relevant transition
	  names(shifts3)[names(shifts3)=="temp"] <- paste(uniShifts[i,1], uniShifts[i,2], sep="to")
	}
	
  #set aside the columns that record transitions, these will be needed soon
	relCols <- names(shifts3)[grep("to", names(shifts3))]
	
	#set up a list to save these results into
	results <- vector("list", length(relCols))
	names(results) <- relCols
	
  #go into a for loop where each iteration you summarize the descendents from the node
	#subtending each transition
	for(i in 1:length(relCols))
	{
	  temp <- list()
	  tempNodes <- shifts3$daughter.node[shifts3[,relCols[i]]==1]
	  
	  #this nested part will pull the taxa that descend each node in tempNodes
	  for(j in 1:length(tempNodes))
	  {
	    temp[[j]] <- geiger::tips(orig.tree, tempNodes[j])
	  }
	  
	  #set the names to reflect the node in question
	  names(temp) <- tempNodes
	  
	  #set these results into the relevant slot in the final results
	  results[[i]] <- temp
	}

	results
}
