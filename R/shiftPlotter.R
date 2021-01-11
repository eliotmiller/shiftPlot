#' Plot a phylogeny with triangles representing clades
#'
#' Use the optimal collapses to create a tree where triangles represent clades.
#'
#' @param orig.tree A phylogeny in ape format. This should be the ORIGINAL tree, not the
#' collapsed tree nor the polytomy tree.
#' @param poly.tree A phylogeny in ape format. This should be the polytomy tree, i.e. the
#' result of a call to polytomyBind.
#' @param coords The result of a call to getCoords. These are returned when the trianglePlotter
#' function is used, which allows a user to save those and re-use them without overwriting the
#' triangle plot by calling getCoords again.
#' @param identifyShifts.obj Result of a call to identifyShifts.
#' @param presence.color What color clades for which the trait is present should be colored.
#' @param absence.color What color clades for which the trait is absent should be colored.
#' @param pt.cex The size of the shift points.
#' 
#' @details Stuff.
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
#' shiftPlotter(phy, polyTree, saveMe, shifts, "red", "black", 1)

#coords should be an object created for poly.tree
shiftPlotter <- function(orig.tree, poly.tree, coords, identifyShifts.obj, presence.color,
	absence.color, pt.cex)
{
	#find all the branching times for the original tree
	times <- ape::branching.times(orig.tree)

	#create a vector to subset for tracking tips
	tipNos <- 1:length(orig.tree$tip.label)

	#save a segment vector to use
	segVector <- 1:length(coords$xx)

	#this is uglier code than it needs to be, but for simplicity, run through all the gains,
	#then all the losses
	for(i in 1:length(identifyShifts.obj$gains))
	{
		#do this in case no gains were found
		if(identifyShifts.obj$gains[[1]][1]=="No gains were found in this reconstruction")
		{
			next()
		}

		#pull taxa for simplicity
		taxa <- identifyShifts.obj$gains[[i]]
		
		#if there's just a single taxon here, then the code below would cause an error.
		#put a dot down on that single taxon and move on to next if that's the case
		if(length(taxa)==1)
		{
		  xCoord <- max(times)
		  tempEdge <- which.edge(poly.tree, as.character(identifyShifts.obj$gains[i]))
		  tempTip <- poly.tree$edge[tempEdge,2]
		  yCoord <- coords$yy[tempTip]
		}

		else
		{
  		#figure out what the tip numbers are that correspond to those taxa in each tree
  		theTipsOrig <- tipNos[orig.tree$tip.label %in% taxa]
  		theTipsPoly <- tipNos[poly.tree$tip.label %in% taxa]
  
  		#find the mrca of those taxa in each tree
  		mrcaNodeOrig <- ape::getMRCA(orig.tree, taxa)
  		mrcaNodePoly <- ape::getMRCA(poly.tree, taxa)
  		
  		#find ALL taxa that descend from poly mrca	
  		allTaxa <- geiger::tips(poly.tree, mrcaNodePoly)
  		allTips <- tipNos[poly.tree$tip.label %in% allTaxa]
  
  		#find the original branching time of that node
  		theTime <- times[names(times)==mrcaNodeOrig]
  
  		#use the time to find actual x-coordinates
  		xCoord <- max(times)-theTime
  
  		#highly experimental way of finding the y-coordinate of the shift.
  		yCoord <- coords$yy[mrcaNodePoly]
		}

		#put the point down
		points(x=xCoord, y=yCoord, cex=pt.cex, col=presence.color, pch=20)
	}	

	for(i in 1:length(identifyShifts.obj$losses))
	{
		#do this in case no losses were found
		if(identifyShifts.obj$losses[[1]][1]=="No losses were found in this reconstruction")
		{
			next()
		}

		#pull taxa for simplicity
		taxa <- identifyShifts.obj$losses[[i]]

		#if there's just a single taxon here, then the code below would cause an error.
		#put a dot down on that single taxon and move on to next if that's the case
		if(length(taxa)==1)
		{
		  xCoord <- max(times)
		  tempEdge <- which.edge(poly.tree, as.character(identifyShifts.obj$losses[i]))
		  tempTip <- poly.tree$edge[tempEdge,2]
		  yCoord <- coords$yy[tempTip]
		}
		
		else
		{
  		#figure out what the tip numbers are that correspond to those taxa in each tree
  		theTipsOrig <- tipNos[orig.tree$tip.label %in% taxa]
  		theTipsPoly <- tipNos[poly.tree$tip.label %in% taxa]
  
  		#find the mrca of those taxa in each tree
  		mrcaNodeOrig <- getMRCA(orig.tree, taxa)
  		mrcaNodePoly <- getMRCA(poly.tree, taxa)
  		
  		#find ALL taxa that descend from poly mrca	
  		allTaxa <- tips(poly.tree, mrcaNodePoly)
  		allTips <- tipNos[poly.tree$tip.label %in% allTaxa]
  
  		#find the original branching time of that node
  		theTime <- times[names(times)==mrcaNodeOrig]
  
  		#use the time to find actual x-coordinates
  		xCoord <- max(times)-theTime
  
  		#highly experimental way of finding the y-coordinate of the shift.
  		yCoord <- coords$yy[mrcaNodePoly]
		}

		#put the point down
		points(x=xCoord, y=yCoord, cex=pt.cex, col=absence.color, pch=20)
	}
}

