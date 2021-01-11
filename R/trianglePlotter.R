#' Plot a phylogeny with triangles representing clades
#'
#' Use the optimal collapses to create a tree where triangles represent clades.
#'
#' @param tree A phylogeny in ape format. This should be the polytomy tree, i.e. the
#' result of a call to polytomyBind.
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
#' @param branches Result of a call to firstBranches.
#' @param presence.color What color clades for which the trait is present should be colored.
#' @param absence.color What color clades for which the trait is absent should be colored.
#' @param label.offset How far away from the tips the clade name should be positioned.
#' @param text.cex The size of the clade labels. 
#' 
#' @details Should change the graphical parameters so there's more space for long clade names.
#' Should not plot the tips. Should allow customization of edge colors. 
#'
#' @author Eliot Miller and Bruce Martin
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
#' #use the polytomyBind function
#' polyTree <- polytomyBind(dropped, result)
#' 
#' branchingResults <- firstBranches(phy, dropped, result)
#' 
#' saveMe <- trianglePlotter(polyTree, dropped, result, branchingResults, "red", "black", 0.05, 0.2)

trianglePlotter <- function(tree, dropped.results, clade.table, branches, identifyShifts.obj,
                            presence.color, absence.color, label.offset, text.cex, root.state)
{
	coords <- getCoords(tree)

	#coords$xx <- coords$xx[-(1:length(tree$tip.label))]
	#coords$yy <- coords$yy[-(1:length(tree$tip.label))]
	
	#print(coords)
	
	#figure out what the branch colors should be
	if(root.state=="present")
	{
	  cols <- rep(presence.color, dim(tree$edge)[1])
	}
	
	else if(root.state=="absent")
	{
	  cols <- rep(absence.color, dim(tree$edge)[1])
	}
	
	else
	{
	  stop("root.state must be one of either present or absent")
	}
	
	#set branches to the correct colors. i suspect this is too simplistic
	#to just do one then the other, but try this for now. probably, there is some complicated
	#nesting that can happen that needs to be anticipated, but seems to work for simplest
	#cases
	for(i in 1:length(identifyShifts.obj$losses))
	{
	  toChange <- which.edge(tree, unlist(identifyShifts.obj$losses[i]))
	  cols[toChange] <- absence.color
	}
	
	for(i in 1:length(identifyShifts.obj$gains))
	{
	  toChange <- which.edge(tree, unlist(identifyShifts.obj$gains[i]))
	  cols[toChange] <- presence.color
	}
	
	#plot the tree structure here. critically, this is going to plot the structure in black,
	#including the tip segments. this might look weird, and it might be a good idea to suppress
	#the plotting of tip segments, and/or customize the color of these segments
	segments(y0=coords$yy[as.vector(t(tree$edge))],
		y1=coords$yy[rep(tree$edge[,2],each=2)],
		x0=coords$xx[rep(tree$edge[,1],each=2)],
		x1=coords$xx[as.vector(t(tree$edge))], lwd=0.1, col=rep(cols,each=2))

	#prep a table for the function below
	groupsDF <- data.frame(species=unlist(dropped.results[[2]]),
		group=rep(dropped.results[[1]]$tip.label, unlist(lapply(dropped.results[[2]], length))),
		stringsAsFactors=FALSE)

	#and merge in the clade states for another function further below
	groupsDF <- merge(groupsDF, clade.table[,c("present","clade")], by.x="group", by.y="clade")

	uniqueGroups <- unique(groupsDF$group)

	#figure out what the right hand x-points of the triangles will be
	rightPoint <- max(ape::branching.times(tree))
	
	for(i in 1:length(uniqueGroups))
	{
		#subset to the taxa in question for that group
		taxa <- groupsDF$species[groupsDF$group==uniqueGroups[i]]

		#identify the tips that belong to those taxa. create a little vector of tip numbers
		tipNos <- 1:length(tree$tip.label)
		theTips <- tipNos[tree$tip.label %in% taxa]

		#if this is NA, then there's only a single taxon, so don't plot a triangle but do label it
		if(is.na(branches[i]))
		{
			text(groupsDF$group[groupsDF$group==uniqueGroups[i]][1], x=rightPoint+label.offset, y=theTips,
				adj=0, cex=text.cex)
			next()
		}

		#figure out what the high and low y-points of the right side of the triangle will be
		highPoint <- min(theTips)
		lowPoint <- max(theTips)
		midPoint <- (highPoint-lowPoint)/2 + lowPoint

		#find the left point
		leftPoint <- rightPoint-branches[i]

		#plot the triangle with the correct color
		if(groupsDF$present[groupsDF$group==uniqueGroups[i]][1]==1)
		{
			polyColor <- presence.color
		}

		else
		{
			polyColor <- absence.color
		}

		polygon(x=c(leftPoint,rightPoint,rightPoint,leftPoint),
			y=c(midPoint,lowPoint,highPoint,midPoint), col=polyColor, border=NA)

		text(groupsDF$group[groupsDF$group==uniqueGroups[i]][1], x=rightPoint+label.offset, y=midPoint,
			adj=0, cex=text.cex)
	}

	coords
}
