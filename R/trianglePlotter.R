#' Plot a phylogeny with triangles representing clades
#'
#' Use the optimal collapses to create a tree where triangles represent clades.
#'
#' @param orig.tree A phylogeny in ape format. This should be the original tree.
#' @param poly.tree A phylogeny in ape format. This should be the polytomy tree, i.e. the
#' result of a call to polytomyBind.
#' @param states.df Data frame in the specified shiftPlot format. Should contain one
#' column named "state", and number for every node in the phylogeny, with the tips
#' above the internal nodes, and no row names. See details and examples. 
#' @param dropped.result The result of a call to dropManyTips.
#' @param oc.result The result of a call to optimalCollapse.
#' Data frame with four columns: node, state, collapse (0/1),
#' and clade. Node provides values indicating which nodes you will collapse (all nodes
#' tipwards from the indicated node will be collapsed). Collapse indicates whether to actually collapse
#' that node or not (allowing a user to manually override the results from optimalCollapse).
#' Clade provides a character string which will be used to rename the collapsed clade.
#' Initially, this is given a generic name based on the node being collapsed, but these can
#' be replaced with a name of the user's choosing. 
#' @param branches Result of a call to firstBranches.
#' @param is.result Result of a call to identifyShifts.
#' @param translation.table A data frame with the following two columns: "state" and "color",
#' the name of the color to use in the plot to indicate that that trait state. Color can be
#' specified in a variety of ways--word, hexadecimal, probably rgb but untested.
#' @param label.offset How far away from the tips the clade name should be positioned.
#' @param text.cex The size of the clade labels.
#' @param pt.cex The size of the dots used to indicate shifts. Set to 0 to suppress points.
#' @param seg.lwd Width of the edges in the plotted phylogeny.
#' 
#' @details Ambitions of making this work for radial trees too.
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
#' 
#' #run the function and don't flip those tips
#' is.result <- identifyShifts(orig.tree=phy, states.df=induced, flip.tips=FALSE)
#' 
#' #create a color, state, and state name translation table
#' translationTable <- data.frame(state=c(0,1), color=c("black","red"))
#' 
#' #make the triangle plot
#' trianglePlotter(orig.tree=phy, poly.tree=polyTree, dropped.result=dropped,
#' oc.result=ocResult, branches=branches, is.result=is.result,
#' label.offset=0.3, text.cex=0.09, translation.table=translationTable,
#' states.df=states, pt.cex=2, seg.lwd=0.1)

trianglePlotter <- function(orig.tree, poly.tree, states.df, dropped.result, oc.result, branches,
                            is.result, translation.table, label.offset, text.cex, pt.cex, seg.lwd)
{
	#set aside the coordinates
  coords <- getCoords(poly.tree)
  
  #identify all shifts you need to change colors downstream of
  allShifts <- as.numeric(unlist(lapply(is.result, names)))

  #split out internal node shifts and tip shifts and handle each separately
  internalShifts <- allShifts[allShifts > length(orig.tree$tip.label)]
  tipShifts <- allShifts[!(allShifts %in% internalShifts)]
  
  #map out the node matches between the main tree and poly tree
  matches <- idMatches(orig.tree, poly.tree, TRUE)
  
  #use this to convert the shifts to the relevant nodes in poly.tree. note that
  #these only refer to internal nodes right now
  internalShiftsPoly <- matches$tree2[matches$tree1 %in% allShifts]
  
  #create a color/state lookup table
  shiftType <- c()
  for(i in 1:length(is.result))
  {
    shiftType <- c(shiftType, rep(names(is.result)[i], length(is.result[[i]])))
  }

  #split and pull the state it changed to
  toState <- strsplit(shiftType, "to")
  
  #pull second element and unlist
  toState <- unlist(lapply(toState, "[", 2))
  
  #create table
  lookup <- data.frame(node=allShifts, state=toState)
  lookup <- merge(lookup, translation.table)
  
  #figure out what the tip numbers (nodes) in the poly tree are of the tips that had shifts
  #as calculated on the real tree
  tempVector <- 1:length(poly.tree$tip.label)
  shiftSpp <- orig.tree$tip.label[tipShifts]
  inPoly <- tempVector[poly.tree$tip.label %in% shiftSpp]
  
  #create a quick data frame of these matches and bind onto previous internal node matches
  matchedTips <- data.frame(tree1=tipShifts, tree2=inPoly)
  matches <- rbind(matches, matchedTips)
  
  #merge in the poly nodes
  lookup <- merge(lookup, matches, by.x="node", by.y="tree1")

  #there can be funny business between identify shift results and the optimal collapse
  #results that I think has to do with marginal tip reconstruction, probably. look for
  #situations where certain shifts found do not match up
  if(length(allShifts) != dim(lookup)[1])
  {
    warning("Consider changing the inferred states of these nodes (or their parents) from the original tree. To do so, change in states.df and re-run identifyShifts")
    print(setdiff(allShifts, lookup$node))
  }
  
  #set aside the poly tip nodes for later
  tipShiftsPoly <- lookup$tree2[!(lookup$tree2 %in% internalShiftsPoly)]
  
  #begin by coloring all the branches the designated color for the root state. 
	#find root state then find color and rep
	rootState <- states.df[(length(orig.tree$tip.label)+1),"state"]
	cols <- rep(translation.table$color[translation.table$state==rootState], dim(poly.tree$edge)[1])

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
	         x1=coords$xx[as.vector(t(poly.tree$edge))], lwd=seg.lwd, col=rep(cols,each=2))
	
	#prep a table for the function below
	groupsDF <- data.frame(species=unlist(dropped.result[[2]]),
		group=rep(dropped.result[[1]]$tip.label, unlist(lapply(dropped.result[[2]], length))),
		stringsAsFactors=FALSE)

	#and merge in the clade states for another function further below
	groupsDF <- merge(groupsDF, oc.result[,c("state","clade")], by.x="group", by.y="clade")

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
		polyColor <- translation.table$color[translation.table$state==tempState]

		#unlist this shift object, but not recursively. see if any shifts correspond to the species
		#defined in taxa. if so, then this is a state shift and you do want to add a point
		unlisted <- unlist(is.result, recursive=FALSE)
		if(any(unlist(lapply(unlisted, function(x) all.equal(taxa,x)))==TRUE))
		{
		  points(x=leftPoint, y=midPoint, pch=20, cex=pt.cex, col=polyColor)
		}
		
		polygon(x=c(leftPoint,rightPoint,rightPoint,leftPoint),
			y=c(midPoint,lowPoint,highPoint,midPoint), col=polyColor, border=NA)

		text(groupsDF$group[groupsDF$group==uniqueGroups[i]][1], x=rightPoint+label.offset, y=midPoint,
			adj=0, cex=text.cex)
	}
  
	#go ahead and add the shift pts. first find all the branching times for the original tree
	times <- ape::branching.times(poly.tree)
	
	#go through the tip shifts first
	for(i in 1:length(tipShiftsPoly))
	{
	  xCoord <- max(times)
	  yCoord <- coords$yy[tipShiftsPoly[i]]
	  #points(x=xCoord, y=yCoord, pch=20, cex=pt.cex, col=lookup$color[lookup$tree2==tipShiftsPoly[i]])
	}
	
	#go through the internal shifts next
	for(i in 1:length(internalShiftsPoly))
	{
	  #find the original branching time of that node
	  theTime <- times[names(times)==internalShiftsPoly[i]]

	  #use the time to find actual x-coordinates
	  xCoord <- max(times)-theTime
	  
	  #pull the y coordinates for that shift
	  yCoord <- coords$yy[internalShiftsPoly[i]]

	  #if xCoord is within some (currently hard-coded) tolerance of the tips, do not
	  #plot the shift, as it looks like a tip shift
    if(all.equal(as.numeric(theTime), 0, tolerance=0.001)==TRUE)
    {
      next()
    }
	  
	  else
	  {
	    points(x=xCoord, y=yCoord, pch=20, cex=pt.cex, col=lookup$color[lookup$tree2==internalShiftsPoly[i]]) 
	  }
	}
}
