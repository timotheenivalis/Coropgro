#' Reformats data on the opportunity for selection
#'
#' Takes data produced by the Fiz() function and reformats in a long-table to be output to a csv file. Optionally can write direct to file
#' 
#' @param Opportunities A list created by the function Fiz()
#' @param lambda A vector of population growth rates created by the function PopGrowthFromID(), if that is to be included in the file. Defaults to NULL
#' @param outputfile File path for where the reformatted dataframe should be saved to. Default is NULL (no file written)
#' 
#' @return A long form data frame containing information on the opportunity for selection over multiple years.
#' 
#' @seealso \code{\link{Fiz}}
#' 
#' @examples
#' data("svdata")
#' 
#' Opportunities <- Fiz(dat = svdata, fitnesstrait = "FitnessTot", 
#' phenotraits = c("Mass", "BL", "TL"), yearname = "Year", 
#' covariates = c("Sex", "Age"), family = "poisson")
#' PopulationGrowth <- PopGrowthFromID(dat=svdata, idname="ID", yearname="Year")  
#' OpportunitiesCSV(Opportunities, lambda=PopulationGrowth, outputfile="Snowvoleopportunities.csv")
#'
#' @export

OpportunitiesCSV<-function(Opportunities, lambda=NULL, outputfile=NULL){

allI<-rbind(Opportunities$iz, Opportunities$Itot)
rownames(allI)<-c(rownames(Opportunities$iz), "Itot")

opp<-data.frame(Year=rep(Opportunities$years, each=nrow(allI)),
                SampleSizes=rep(Opportunities$SampleSizes, each=nrow(allI)),
                sq=rep(Opportunities$sq, each=nrow(allI)),
                Opportunity=rep(rownames(allI), length(Opportunities$years)),
                I_value=c(allI))

if(!is.null(lambda)){
  opp$lambda<-NA
  opp$lambda<-lambda[match(opp$Year, names(lambda))]
}

if(!is.null(outputfile)){
  utils::write.csv(opp, file=outputfile, row.names = FALSE)
}

return(opp)
}