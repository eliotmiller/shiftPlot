#' Drop all tips but one per named clade
#'
#' Collapse a phylogeny by node and rename according to a clade table.
#'
#' @param tree Phylogeny in ape format.
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
#' @return A list of result. The first object is the collapsed tree, the second is a list of
#' which tips/species descend from each collapsed node.
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

dropManyTips <- function(tree, clade.table)
{
	#set a tree aside to modify
	modTree <- tree

	#create a list to tabulate which tips (and how many) you're dropping
	drops <- list()

	for(i in 1:dim(clade.table)[1])
	{
		if(clade.table$collapse[i]==0)
		{
			drops[[i]] <- NA
			next()
		}

		else
		{
			#rename the first tip descending from that node. first define all of them.
			allTips <- geiger::tips(tree, clade.table$node[i])
			drops[[i]] <- allTips

			#now in the set aside tree, change just that tip label to the clade name
			modTree$tip.label[modTree$tip.label == allTips[1]] <- clade.table$clade[i]
			
			#now drop all tips in allTips. one of those is no longer in modTree, because
			#you changed its name, but this doesn't seem to matter (fxn flexible to it)
			modTree <- ape::drop.tip(modTree, allTips)
		}
	}

	results <- list(tree=modTree, dropped=drops)
	results
}
