#' Drop duplicated matched nodes
#'
#' Identify which of the duplicated nodes you want to retain.
#'
#' @param tree1 The first tree.
#' @param tree2 The second tree.
#' @param id.results Results from a run of the idMatches function.
#' @param least.inclusive Do you want to keep the smallest or largest set of
#' matched taxa. See details. 
#' 
#' @details The point of this function is to drop duplicated nodes that appear more than once
#' in the results for tree 2. If there are multiple nodes in tree 2 which could potentially be mapped to
#' a node in tree 1, do you want to keep the smallest or the largest one? for example:
#' ((A,B),Z) vs (((A,B),Z),C), where C doesn't occur in tree 1. Do you want the matched node
#' from tree 2 to be A,B,Z or A,B,Z,C.
#'
#' @author Eliot Miller

equivCut <- function(tree1, tree2, id.results, least.inclusive)
{
  #identify nodes that appear more than once
  dups <- unique(id.results[,2][duplicated(id.results[,2])])
  
  #if there aren't any, return the original results
  if(length(dups)==0)
  {
    return(id.results)
  }
  
  #loop through, add details on size of clade to data frame
  for(i in 1:length(dups))
  {
    #set aside all the results that pertain to that node
    setAside <- id.results[id.results[,2]==dups[i],]
    
    #go through each node and figure out how many taxa descend from that node
    temp <- data.frame(matrix(nrow=dim(setAside)[1], ncol=2))
    names(temp) <- c("node","no.taxa")
    for(j in 1:dim(setAside)[1])
    {
      temp[j,"node"] <- setAside[j,"tree1"]
      temp[j,"no.taxa"] <- length(extract.clade(tree1, setAside[j,"tree1"])$tip.label)
    }
    
    #now you can decide which of these to keep
    if(least.inclusive)
    {
      keep <- temp$node[temp$no.taxa==min(temp$no.taxa)]
      toDrop <- setdiff(temp$node, keep)
      id.results <- id.results[!(id.results$tree1 %in% toDrop),]
    }
    
    else
    {
      keep <- temp$node[temp$no.taxa==max(temp$no.taxa)]
      toDrop <- setdiff(temp$node, keep)
      id.results <- id.results[!(id.results$tree1 %in% toDrop),]
    }
  }
  
  id.results
}
