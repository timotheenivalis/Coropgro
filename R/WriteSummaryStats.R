#' Write opportunities and population growth to csv file
#'
#' Write year-specific opportunities for selection and population growth rate to a csv file
#'
#' @param Opportunities Output from a call to Fiz().
#' @param PopulationGrowth Output from a call to PopGrowthFromID().
#' @param FileName Name of the file where the output sould be written (including directory if not to be printed in the working directory).
#' @param write Boolean. If TRUE, write the file. If FALSE, return the output instead.
#' @param fitnesstrait Name of the fitness trait used for opportunity for selection.
#' @param family GLM family ("gaussian", "poisson" or "binomial")
#' 
#' @return Write to a file or return a data frame
#'
#' @seealso \code{\link{Fiz}} and \code{\link{PopGrowthFromID}}
#' 
#' @examples
#'
#' data(svdata)
#' Opportunities <- Fiz(dat=svdata, fitnesstrait="Phi",
#'  phenotraits=c("Mass", "BL", "TL"), family="binomial")
#' PopulationGrowth <- PopGrowthFromID(dat=svdata, idname="ID", yearname="Year")
#' 
#' WriteSummaryStats(Opportunities=Opportunities, PopulationGrowth=PopulationGrowth )
#'
#' @export
WriteSummaryStats<-function(Opportunities, PopulationGrowth, FileName="SummaryStatsIpop.csv", write=TRUE, fitnesstrait=NA, family=NA){
  
  outd <- data.frame(years= Opportunities$years, 
             SampleSizes = as.vector(Opportunities$SampleSizes),
             Itot = Opportunities$Itot,
             sq=Opportunities$sq,
             lambda = as.vector(PopulationGrowth),
             fitnesstrait=fitnesstrait,
             family=family)
  
  outd2 <-cbind(outd, t(Opportunities$iz))
  
  if (write)
  {
    utils::write.csv(x = outd2, file = FileName, row.names = FALSE)
  }else{
    return(outd2)
  }
  
}#end function