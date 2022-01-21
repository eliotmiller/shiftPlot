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
#' @param identifyShifts.obj Result of a call to identifyShifts.
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
#' collapsed <- optimalCollapse(phy, induced, FALSE)
#' 
#' shifts <- identifyShifts(phy, induced, FALSE)
#' 
#' #collapse the tree
#' dropped <- dropManyTips(phy, collapsed)
#' 
#' #use the polytomyBind function
#' polyTree <- polytomyBind(dropped, collapsed)
#' 
#' branchingResults <- firstBranches(phy, dropped, collapsed)
#' 
#' saveMe <- trianglePlotter(tree=polyTree, dropped.results=dropped, clade.table=collapsed,
#' branches=branchingResults, identifyShifts.obj=shifts,
#' presence.color="red", absence.color="black",
#' label.offset=0.3, text.cex=0.09,
#' root.state="present")
#' 

trianglePlotter <- function(tree, poly.tree, states.df, dropped.results, clade.table, branches,
                            identifyShifts.obj, translation.table, label.offset, text.cex)
{
	#set aside the coordinates
  coords <- getCoords(poly.tree)
  
  #identify all shifts you need to change colors downstream of
  allShifts <- as.numeric(unlist(lapply(identifyShifts.obj, names)))

  #split out internal node shifts and tip shifts and handle each separately
  internalShifts <- allShifts[allShifts > length(tree$tip.label)]
  tipShifts <- allShifts[!(allShifts %in% internalShifts)]
  
  #map out the node matches between the main tree and poly tree
  matches <- idMatches(tree, poly.tree, TRUE)
  
  #use this to convert the shifts to the relevant nodes in poly.tree. note that
  #these only refer to internal nodes right now
  internalShiftsPoly <- matches$tree2[matches$tree1 %in% allShifts]
  
  #create a color/state lookup table
  shiftType <- c()
  for(i in 1:length(identifyShifts.obj))
  {
    shiftType <- c(shiftType, rep(names(identifyShifts.obj)[i], length(identifyShifts.obj[[i]])))
  }

  #split and pull the state it changed to
  toState <- strsplit(shiftType, "to")
  
  #pull second element and unlist
  toState <- unlist(lapply(toState, "[", 2))
  
  #create table
  lookup <- data.frame(node=allShifts, state=toState)
  lookup <- merge(lookup, translation.table[,c("new","color")], by.x="state", by.y="new")
  
  #figure out what the tip numbers (nodes) in the poly tree are of the tips that had shifts
  #as calculated on the real tree
  tempVector <- 1:length(poly.tree$tip.label)
  shiftSpp <- tree$tip.label[tipShifts]
  inPoly <- tempVector[poly.tree$tip.label %in% shiftSpp]
  
  #create a quick data frame of these matches and bind onto previous internal node matches
  matchedTips <- data.frame(tree1=tipShifts, tree2=inPoly)
  matches <- rbind(matches, matchedTips)
  
  #merge in the poly nodes
  lookup <- merge(lookup, matches, by.x="node", by.y="tree1")

  #set aside the poly tip nodes for later
  tipShiftsPoly <- lookup$tree2[!(lookup$tree2 %in% internalShiftsPoly)]
  
  #begin by coloring all the branches the designated color for the root state. 
	#find root state then find color and rep
	rootState <- states.df[(length(tree$tip.label)+1),"state"]
	cols <- rep(translation.table$color[translation.table$new==rootState], dim(poly.tree$edge)[1])

  #go into a for loop the length of internalShifts
	for(i in 1:length(internalShiftsPoly))
	{
	  #find all the edges that connect back from the tips that descend from shift i
    spp <- extract.clade(poly.tree, internalShiftsPoly[i])$tip.label
    allEdges <- which.edge(poly.tree, spp)
    
    #find all nodes that descend from this internal node. even though tips is FALSE
    #here, it also returns tip nodes
    allNodes <- geiger:::.get.descendants.of.node(internalShiftsPoly[i], poly.tree, tips=FALSE)
	    
    #subset to any that are in internalShiftsPoly
    allNodes <- allNodes[allNodes %in% internalShiftsPoly]
    
    #if there are any internal nodes with shifts that descend from this node, deal with them
    if(length(allNodes) > 0)
    {
      #go into a for loop where you find the edges for each descendent node j that has a shift
      #on it and discard those edges from allEdges
      for(j in 1:length(allNodes))
      {
        tempSpp <- extract.clade(poly.tree, allNodes[j])$tip.label
        tempEdges <- which.edge(poly.tree, tempSpp)
        
        #remove the edges you identified
        allEdges <- allEdges[!(allEdges %in% tempEdges)]
      }
    }
    
    #figure out which color is relevant for this shift, and set the relevant edges to it
    cols[allEdges] <- lookup$color[lookup$tree2==internalShiftsPoly[i]]
	}

	#go into a for loop for the tip shifts
	for(i in 1:length(tipShiftsPoly))
	{
	  #find the edge for just that tip
	  allEdges <- which.edge(poly.tree, poly.tree$tip.label[tipShiftsPoly[i]])
	  cols[allEdges] <- lookup$color[lookup$tree2==tipShiftsPoly[i]]
	}
	
	#plot the tree structure here. critically, this is going to plot the structure in black,
	#including the tip segments. this might look weird, and it might be a good idea to suppress
	#the plotting of tip segments, and/or customize the color of these segments
	segments(y0=coords$yy[as.vector(t(poly.tree$edge))],
	         y1=coords$yy[rep(poly.tree$edge[,2],each=2)],
	         x0=coords$xx[rep(poly.tree$edge[,1],each=2)],
	         x1=coords$xx[as.vector(t(poly.tree$edge))], lwd=0.1, col=rep(cols,each=2))
	
	#prep a table for the function below
	groupsDF <- data.frame(species=unlist(dropped.results[[2]]),
		group=rep(dropped.results[[1]]$tip.label, unlist(lapply(dropped.results[[2]], length))),
		stringsAsFactors=FALSE)

	#and merge in the clade states for another function further below
	groupsDF <- merge(groupsDF, clade.table[,c("state","clade")], by.x="group", by.y="clade")

	uniqueGroups <- unique(groupsDF$group)

	#figure out what the right hand x-points of the triangles will be
	rightPoint <- max(ape::branching.times(poly.tree))
	
	for(i in 1:length(uniqueGroups))
	{
		#subset to the taxa in question for that group
		taxa <- groupsDF$species[groupsDF$group==uniqueGroups[i]]

		#identify the tips that belong to those taxa. create a little vector of tip numbers
		tipNos <- 1:length(poly.tree$tip.label)
		theTips <- tipNos[poly.tree$tip.label %in% taxa]

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

		#find the color the polygon should be
		tempState <- unique(groupsDF$state[groupsDF$group==uniqueGroups[i]])
		polyColor <- translation.table$color[translation.table$new==tempState]

		polygon(x=c(leftPoint,rightPoint,rightPoint,leftPoint),
			y=c(midPoint,lowPoint,highPoint,midPoint), col=polyColor, border=NA)

		text(groupsDF$group[groupsDF$group==uniqueGroups[i]][1], x=rightPoint+label.offset, y=midPoint,
			adj=0, cex=text.cex)
	}

	coords
}
