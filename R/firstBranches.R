#' Drop all tips but one per named clade
#'
#' Collapse a phylogeny by node and rename according to a clade table.
#'
#' @param tree A phylogeny in ape format. This should be the ORIGINAL tree, not the
#' collapsed tree nor the polytomy tree.
#' @param dropped.results The result of a call to dropManyTips.
#' @param clade.table Data frame with four columns: node, present (0/1), collapse (0/1),
#' and clade. Node provides values indicating which nodes you will collapse (all nodes
#' tipwards from the indicated node will be collapsed). Present indicates whether the
#' tips being collapsed do or do not have the trait in question (I think this column is
#' not actually used and can be skipped). Collapse indicates whether to actually collapse
#' that node or not (allowing a user to manually override the results from optimalCollapse).
#' Clade provides a character string which will be used to rename the collapsed clade.
#' Initially, this is given a generic name based on the node being collapsed, but these can
#' be replaced with a name of the user's choosing. This input (clade.table) is the output
#' of optimalCollapse.
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
#' #find the optimal collapse configuration
#' result <- optimalCollapse(phy, induced, FALSE)
#' 
#' #collapse the tree
#' dropped <- dropManyTips(phy, result)
#' 
#' #find the branching times
#' branches <- firstBranches(phy, dropped, result)

firstBranches <- function(tree, dropped.results, clade.table)
{
	#prep a table for the function below
	groupsDF <- data.frame(species=unlist(dropped.results[[2]]),
		group=rep(dropped.results[[1]]$tip.label, unlist(lapply(dropped.results[[2]], length))),
		stringsAsFactors=FALSE)

	#and merge in the clade states for another function further below
	groupsDF <- merge(groupsDF, clade.table[,c("state","clade")], by.x="group", by.y="clade")

	results <- c()
	uniqueGroups <- unique(groupsDF$group)

	#solve once for the branching times
	times <- ape::branching.times(tree)

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
		mrcaNode <- ape::getMRCA(tree, taxa)

		#find the branching time of that node
		results[i] <- times[names(times)==mrcaNode]
	}

	results
}
