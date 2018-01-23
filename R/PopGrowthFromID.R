#' Calculate population growth rate between years
#'
#' Uses minimal and maximal individual presence to compute changes in population size. 
#' 
#' @param dat Data-frame. Contains all necessary data.
#' @param idname Column name containing individual unique identifiers
#' @param yearname Column name containing years (or other time step). 
#' 
#' #' @return A table of annual population growth rate.
#' 
#' #' @examples
#' data("svdata")
#' PopGrowthFromID(dat=svdata, idname="ID", yearname="Year")
#'
#' @export
PopGrowthFromID <- function(dat, idname="ID", yearname="Year")
{
  if(!is.numeric(dat[,yearname])) # check we don't get nonsense from a factor or something
  {
    stop("The column yearname is not numeric.")
  }
  
  yearpresence <- unlist(tapply(X = dat[,yearname], INDEX = dat[,idname], function(x) {
    seq(min(x), max(x),by = 1)
  }))#expands individual presence to years in between first and last occurence
  
  popsizes <- table(yearpresence, dnn = "From year y to y+1")
  popsizes <- popsizes[order(names(popsizes))] # be sure the years are in the right order
  lambda <- popsizes[-1]/popsizes[-length(popsizes)] # get growth rates
  names(lambda) <- names(popsizes[-length(popsizes)])#rename years to be consistent with opportunity for selection
  
  return(lambda)
}
