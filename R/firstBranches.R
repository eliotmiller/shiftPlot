#' Drop all tips but one per named clade
#'
#' Collapse a phylogeny by node and rename according to a clade table.
#'
#' @param orig.tree A phylogeny in ape format. This should be the original tree.
#' @param dropped.result The result of a call to dropManyTips.
#' @param oc.result The result of a call to optimalCollapse.
#' Data frame with four columns: node, state, collapse (0/1),
#' and clade. 
#' 
#' @details Stuff.
#'
#' @return A vector providing the time of the first branching event for each clade in the
#' tree that was collapsed.
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
#' ocResult <- optimalCollapse(orig.tree=phy, states.df=induced, flip.tips=FALSE)
#' 
#' #run the dropManyTips fxn
#' dropped <- dropManyTips(orig.tree=phy, oc.result=ocResult)
#' 
#' #run the polytomyBind function
#' polyTree <- polytomyBind(dropped.result=dropped, oc.result=ocResult)
#' 
#' #run the firstBranches fxn
#' branches <- firstBranches(orig.tree=phy, dropped.result=dropped, oc.result=ocResult)

firstBranches <- function(orig.tree, dropped.result, oc.result)
{
	#have bumped into a bug before which I believe is due to using trees with existing node
  #labels, then having corHMM strip them out. throw an error if the tree has node labels
  if(length(orig.tree$node.label > 0))
  {
    stop("Your original tree has node labels. Remove the node labels and rerun.")
  }
  
  #prep a table for the function below
	groupsDF <- data.frame(species=unlist(dropped.result[[2]]),
		group=rep(dropped.result[[1]]$tip.label, unlist(lapply(dropped.result[[2]], length))),
		stringsAsFactors=FALSE)

	#and merge in the clade states for another function further below
	groupsDF <- merge(groupsDF, oc.result[,c("state","clade")], by.x="group", by.y="clade")

	results <- c()
	uniqueGroups <- unique(groupsDF$group)

	#solve once for the branching times
	times <- ape::branching.times(orig.tree)

	for(i in 1:length(uniqueGroups))
	{
		#subset to the taxa in question for that group
		taxa <- groupsDF$species[groupsDF$group==uniqueGroups[i]]

		#if there's only a single species, return NA for now (probably want to revisit this)
		if(length(taxa)==1)
		{
			results[i] <- NA
			next()
		}

		#find the MRCA of those species
		mrcaNode <- ape::getMRCA(orig.tree, taxa)

		#find the branching time of that node
		results[i] <- times[names(times)==mrcaNode]
	}

	results
}
