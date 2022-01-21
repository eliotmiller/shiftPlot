#' Identify matched nodes between two trees
#'
#' Identify nodes with matching descendant taxa between two phylogenies.
#'
#' @param tree1 The first tree.
#' @param tree2 The second tree.
#' @param least.inclusive Do you want to keep the smallest or largest set of
#' matched taxa. See details. 
#' 
#' @details If there are multiple nodes in tree 2 which could potentially be mapped to
#' a node in tree 1, do you want to keep the smallest or the largest one? for example:
#' ((A,B),Z) vs (((A,B),Z),C), where C doesn't occur in tree 1. Do you want the matched node
#' from tree 2 to be A,B,Z or A,B,Z,C.
#'
#' @author Eliot Miller
#' 
#' @export
#' 
#' @examples
#' first <- sim.bdtree(n=50)
#' second <- sim.bdtree(n=1000)
#' third <- sim.bdtree(n=2000)
#' #always returns equivalent no matter which tree is first if you set it to least inclusive = TRUE
#' system.time(test <- idMatches(first,second))
#' system.time(test2 <- idMatches(second,first))
#' setequal(test[,1], test2[,2])
#' setequal(test[,2], test2[,1])

idMatches <- function(tree1, tree2, least.inclusive=TRUE)
{
  #set the results table up
  results <- tableSetup(tree1, tree2)
  
  #define tree1 nodes here
  nodes1 <- (length(tree1$tip.label)+1):max(tree1$edge)
  
  for(i in 1:length(nodes1))
  {
    #set the correct row in the tree1 column to the node in question
    results[i,1] <- nodes1[i]
    
    #find the descendants of this node
    tips1 <- tips(tree1, nodes1[i])
    
    #drop these to just tips that are in tree 2
    tips1 <- tips1[tips1 %in% tree2$tip.label]
    
    #if there's nothing left, set to no match and move on
    if(length(tips1)==0)
    {
      results[i,2] <- "no.match"
      next()
    }
    
    #find the MRCA of those tips in tree 2
    #if there's only a single taxon left, pull the node it descends from (getMRCA will fail)
    if(length(tips1==1))
    {
      edge2 <- which.edge(tree2, tips1)
      mrca2 <- tree2$edge[edge2,][1]
    }
    else
    {
      mrca2 <- getMRCA(tree2, tips1)
    }
    
    #find the descendants of that node
    tips2 <- tips(tree2, mrca2)
    
    #drop to just taxa that are in tree1
    tips2 <- tips2[tips2 %in% tree1$tip.label]
    
    #if these tips are equivalent, the nodes match
    if(setequal(tips1, tips2))
    {
      results[i,2] <- mrca2
    }
    
    else
    {
      results[i,2] <- "no.match"
    }
  }
  
  #consider what you want to do next. do you want to, e.g., add a row for
  #every not-yet-mentioned tree2 node and add "no.match" to the first column,
  #or do you want to return as is, or do you want to only return matches? for
  #now, drop to only matches and re-class both as integer so can demonstrate
  #equality no matter which tree is first or second
  results <- results[results[,2] != "no.match",]
  results[,2] <- as.integer(results[,2])
  
  #run the equivCut function
  results <- equivCut(tree1, tree2, results, least.inclusive)
  results
}
