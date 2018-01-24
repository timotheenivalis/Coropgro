#' Correlations between opportunities for selection and population growth rates
#'
#' Calculate the year-based correlations between the various opportunities for selection (total, and trait-specific), and population growth rates
#' 
#' @param Opportunities A list created by the function Fiz()
#' @param PopulationGrowth A table created by the function PopGrowthFromID()
#' 
#' @return A data frame with correlations and their standard errors for the different types of opportunity for selection
#' 
#' @seealso \code{\link{Fiz}}, \code{\link{PopGrowthFromID}} and \code{\link{cor_with_se}}
#' 
#' @examples
#' data("svdata")
#' 
#' Opportunities <- Fiz(dat = svdata, fitnesstrait = "FitnessTot", 
#' phenotraits = c("Mass", "BL", "TL"), yearname = "Year", 
#' covariates = c("Sex", "Age"), family = "poisson")
#' PopulationGrowth <- PopGrowthFromID(dat=svdata, idname="ID", yearname="Year")  
#' AllCorrelations(Opportunities = Opportunities, PopulationGrowth = PopulationGrowth)
#'
#' @export
AllCorrelations <- function(Opportunities, PopulationGrowth)
{
  if(mean(Opportunities$years == names(PopulationGrowth))<1){stop("Year labels in Opportunities and PopulationGrowth do not match. Please check the two are in pair-wise match.")}
  
  wt <- sqrt((Opportunities$SampleSizes-1)/2) # inverse of the standard error factor for a variance estimator
  
  correlations <- data.frame(Opportunity=NA, CorrelationOL=NA, SEOL=NA, CorrelationOIL=NA, SEOIL=NA)#I for inverse lambda
  
  correlations[1,1] <- "Total"
  correlations[1,2:3] <- c(cor_with_se(x = Opportunities$Itot, y = PopulationGrowth, wt = wt))
  correlations[1,4:5] <- c(cor_with_se(x = Opportunities$Itot, y = 1/PopulationGrowth, wt = wt))
  
  for (i in 1:nrow(Opportunities$iz))
  {
    notNAyears <- which(!is.na(Opportunities$iz[i,]))# there might be missing values for some traits on some years
    correlations[1+i,1] <- rownames(Opportunities$iz)[i]
    correlations[1+i,2:3] <- c(cor_with_se(x = Opportunities$iz[i,notNAyears], 
                                           y = PopulationGrowth[notNAyears],
                                           wt = wt[notNAyears]))
    correlations[1+i,4:5] <- c(cor_with_se(x = Opportunities$iz[i,notNAyears], 
                                           y = 1/PopulationGrowth[notNAyears],
                                           wt = wt[notNAyears]))
  }
  return(correlations)
}