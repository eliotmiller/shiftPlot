#' Process phytools simmaps
#'
#' Convert a list of phytools simmaps to a states.df
#'
#' @param simmaps A list of simmaps; the result of a call to phytools::make.simmap
#' where nsim > 1.
#' @param translation.table A data frame with the following two columns: "state" and "color",
#' the name of the color to use in the plot to indicate that that trait state. Color can be
#' specified in a variety of ways--word, hexadecimal, probably rgb but untested.
#' 
#' @details Take a list of simmaps, summarize them, and converts the results into a
#' states.df suitable for downstream use in shiftPlot.
#'
#' @return A new states.df.
#' 
#' @author Eliot Miller
#'
#' @export
#' 

#here is a function to parse simmaps into right format for shiftPlot
processMaps <- function(simmaps, translation.table)
{
  #summarize the simmaps
  summ <- summary(simmaps)
  
  #bind the tips and node internal states like in corHMM
  stateDF <- rbind(summ$tips, summ$ace)
  
  #set this aside and fill accordingly
  temp <- stateDF
  
  #set the max state to 1, all other to 0
  stateDF[] <- 0
  for(i in 1:dim(temp)[1])
  {
    stateDF[i,][temp[i,] == max(temp[i,])] <- 1
  }
  
  #set up a state df in the format you need for shiftPlot
  states <- as.data.frame(matrix(nrow=dim(stateDF)[1], ncol=1, 0))
  names(states) <- "state"
  
  #loop through and fill with the new states
  for(i in 1:dim(translation.table)[1])
  {
    states$state[stateDF[,i]==1] <- i
  }
  
  states
}
