#' Drop all tips but one per named clade
#'
#' Collapse a phylogeny by node and rename according to a clade table.
#'
#' @param dropped.result The result of a call to dropManyTips.
#' @param oc.result Data frame with four columns: node, present (0/1), collapse (0/1),
#' and clade. Node provides values indicating which nodes you will collapse (all nodes
#' tipwards from the indicated node will be collapsed). Present indicates whether the
#' tips being collapsed do or do not have the trait in question (I think this column is
#' not actually used and can be skipped). Collapse indicates whether to actually collapse
#' that node or not (allowing a user to manually override the results from optimalCollapse).
#' Clade provides a character string which will be used to rename the collapsed clade.
#' Initially, this is given a generic name based on the node being collapsed, but these can
#' be replaced with a name of the user's choosing. This input (oc.result) is the output
#' of optimalCollapse.
#' 
#' @details Stuff.
#'
#' @return A tree with the tips from a collapsed node bound in as polytomies.
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

polytomyBind <- function(dropped.result, oc.result)
{
	#prep a table for the function below
	groupsDF <- data.frame(species=unlist(dropped.result[[2]]),
		group=rep(dropped.result[[1]]$tip.label, unlist(lapply(dropped.result[[2]], length))),
		stringsAsFactors=FALSE)

	#and merge in the clade states for another function further below
	groupsDF <- merge(groupsDF, oc.result[,c("state","clade")], by.x="group", by.y="clade")

	uniqueGroups <- unique(groupsDF$group)

	#set this aside here
	tree <- dropped.result[[1]]

	for(i in 1:length(uniqueGroups))
	{
		#subset to the taxa in question for that group
		taxa <- groupsDF$species[groupsDF$group==uniqueGroups[i]]

		#create a tree on the fly
		treeContents <- paste(paste(taxa, ":0", sep=""), collapse=",")
		cat("(", treeContents, ")", ";", file="polytomyBindTemp.tre", sep="\n")

		#read it back in
		toBindIn <- read.tree("polytomyBindTemp.tre")

		#find the correct node for the tip to bind it to
		node <- tree$edge[ape::which.edge(tree, uniqueGroups[i]),2]

		#bind the tree in
		tree <- ape::bind.tree(x=tree, y=toBindIn, where=node)
	}

	#remove that temporary file
	file.remove("polytomyBindTemp.tre")

	#also ladderize the tree, write it out, then read it back and remove
	tree <- ape::ladderize(tree)
	ape::write.tree(tree, "deleteComplete.tre")
	tree <- ape::read.tree("deleteComplete.tre")
	file.remove("deleteComplete.tre")
	tree
}
