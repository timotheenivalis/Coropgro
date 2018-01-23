#' Analyses the correlation between the opportunity for selection and population growth from raw  (individual level) data
#'
#' Starting from individual level data, first calculates population growth rates and opportunity for selection in each year (and for each fitness and phenotypic trait specified). The correlation between these measures is then calculated. 
#' 
#'
#' @param dat Data-frame. Contains all necessary data.
#' @param fitnesstrait Character vector (length>=1). Name of the column(s) containing fitness data.
#' @param phenotraits Character vector (length>=1). Name of the column(s) containing phenotypic data.
#' @param idname Character (length=1). Name of the column containing the individual unique identifiers
#' @param yearname Character (length=1). Name of the column containing the time-step identity (typically year). Default="Year".
#' @param MinDataPoints Integer. Minimal number of data points for a trait on one year. If there are fewer data, the trait is not considered for that year. Default=5.
#' @param covariates NULL of character vector (length>=1). Name of the column(s) containing covariates which effect should be removed from calculations. Default=NULL.
#' @param NAsToMean Boolean. Should missing phenotypic values be changed to the trait mean. Default=TRUE.
#' @param glm_family Character vector. Family of the (G)LMs to be fitted. Possible values are "gaussian", "poisson", "binomial". Default="gaussian". If there is more than one fitness measure considered this should be the same length as the number of fitness measures. If only one family is passed, calculations for all fitness measures are assumed to follow this distribution.
#' @param insertmissing boolean. Should missing individual data be added (individuals not present in a year but known to be in the population as they are present in a subsequent year). Defaults to FALSE
#' @param survivalname Character (length=1). Name of the column containing the survival data. Defaults to NULL. Does not need to be changed from default if insertmissing=FALSE
#' @param outputfile Character (length=1). File path that the resulting data should be written to, if so required. Defaults to NULL (no file is saved).
#' @param SummaryStatsfiles Name of a csv file with summary statistics for opportunity for selection, for all fitness traits, and population growth rate. If NULL (default) does not write it.
#' 
#' @return A data frame with correlations and their standard errors for the different types of opportunity for selection for each of the fitness measures passed to the function.
#' 
#' @seealso \code{\link{Fiz}}, \code{\link{PopGrowthFromID}}, \code{\link{cor_with_se}}, \code{\link{InsertMissing}} and \code{\link{AllCorrelations}}
#' 
#' @examples
#' data("svdata")
#' 
#' AnalyseAll(svdata, fitnesstrait=c("FitnessTot","Phi"), phenotraits=c("Mass","BL","TL"), 
#' idname="ID", yearname="Year", glm_family =c("gaussian", "binomial"),
#'  insertmissing=TRUE, survivalname=NULL)
#' 
#' @export
AnalyseAll<-function(dat, fitnesstrait, phenotraits, idname="ID", yearname="Year", MinDataPoints=5, covariates=NULL, NAsToMean=TRUE, glm_family="gaussian", survivalname=NULL, insertmissing=FALSE, outputfile=NULL, SummaryStatsfiles=NULL){
  
  cordata<-data.frame(Fitness = rep(fitnesstrait, each=(length(phenotraits)+1)),
                      Opportunity = rep(c("Total",phenotraits), length(fitnesstrait)),
                      Correlation = rep(NA, length(fitnesstrait)*(length(phenotraits)+1)),
                      SE = rep(NA, length(fitnesstrait)*(length(phenotraits)+1)),
                      glm_family = rep(NA, length(fitnesstrait)*(length(phenotraits)+1)))
  
  if(insertmissing){
  dat <- InsertMissing(data = dat, idname, yearname, survivalname, verbose = TRUE)
  }
  
  if(length(glm_family)==1 & length(fitnesstrait)>1){
    
    glm_family<-rep(glm_family, length(fitnesstrait))
    cordata$glm_family<-rep(glm_family, each=(length(phenotraits)+1))
    
    print("Only one glm distribution family specified - assumming that all fitness components passed are this family and running GLMs using this distribution")
    
  }else if(length(glm_family)==length(fitnesstrait)){
    
    cordata$glm_family<-rep(glm_family, each=(length(phenotraits)+1))
  
  }else if(length(glm_family)!=length(fitnesstrait)){
      stop("The number of glm distribution families does not match the number of fitness traits specified")
    }
  
  PopulationGrowth<-PopGrowthFromID(dat=dat, idname=idname, yearname=yearname)

  if(!is.null(SummaryStatsfiles))
  {
    outsummary<- data.frame(years=NULL, 
                          SampleSizes = NULL,
                          Itot = NULL,
                          sq=NULL,
                          lambda = NULL,
                          fitnesstrait=NULL,
                          family=NULL)
  }
  
  for(i in 1:length(fitnesstrait)){
    
    fit<-fitnesstrait[i]
    fam<-glm_family[i]
    Opportunities<-Fiz(dat = dat, fitnesstrait = fit, phenotraits = phenotraits, yearname="Year",
                       MinDataPoints =  MinDataPoints, covariates = covariates,NAsToMean =  NAsToMean,family =  fam)
    
    if(!is.null(SummaryStatsfiles))
    {
      outsummary <- rbind(outsummary,WriteSummaryStats(Opportunities, PopulationGrowth,write=FALSE, fitnesstrait=fit, family=fam))
    }
    
    AllCor<-AllCorrelations(Opportunities = Opportunities, PopulationGrowth =  PopulationGrowth)
    
    cordata$Correlation[which(cordata$Fitness==fitnesstrait[i] & cordata$glm_family==fam)]<-AllCor$Correlation
    cordata$SE[which(cordata$Fitness==fitnesstrait[i] & cordata$glm_family==fam)]<-AllCor$SE
  }
  
  if(!is.null(outputfile)){
    utils::write.csv(x = cordata, file=outputfile, row.names = FALSE)
  }
  if(!is.null(SummaryStatsfiles)){
    utils::write.csv(x = outsummary, file = SummaryStatsfiles, row.names = FALSE)
  }
  
  return(cordata)
  
}