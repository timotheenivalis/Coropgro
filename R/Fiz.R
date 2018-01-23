#' Snow vole data
#'
#' A dataset containing individual-based phenotypic data from a snow vole population. The variables are as follows:
#'
#' \itemize{
#'   \item ID Individual unique identifier
#'   \item Year Year of measurement
#'   \item Sex Factor, "Female" or "Male"
#'   \item Age Factor, "A" or "J"
#'   \item Mass Body mass in grams
#'   \item BL Body length in cm
#'   \item TL Tail length in cm
#'   \item Phi Boolean, survival to the next year
#'   \item RhoZ Reproductive success on the next year
#'   \item FitnessTot Sum of RhoZ and Phi
#' }
#'
#' @docType data
#' @keywords datasets
#' @name svdata
#' @usage data(svdata)
#' @format A data frame with 1423 rows and 10 variables
NULL

#' Calculates trait-specific opportunity for selection
#'
#' Calculates trait-specific opportunity for selection
#'
#' @param dat Data-frame. Contains all necessary data.
#' @param fitnesstrait Character (length=1). Name of the column containing fitness data.
#' @param phenotraits Character vector (length>=1). Name of the column(s) containing phenotypic data.
#' @param yearname Character (length=1). Name of the column containing the time-step identity (typically year). Default="Year".
#' @param MinDataPoints Integer. Minimal number of data points for a trait on one year. If there are fewer data, the trait is not considered for that year. Default=5.
#' @param covariates NULL of character vector (length>=1). Name of the column(s) containing covariates which effect should be removed from calculations. Default=NULL.
#' @param NAsToMean Boolean. Should missing phenotypic values be changed to the trait mean. Default=TRUE.
#' @param family Character. Family of the (G)LMs to be fitted. Possible values are "gaussian", "poisson", "binomial". Default="gaussian".
#' 
#' @return A list.
#'
#' @examples
#' data("svdata")
#' Fiz(dat=svdata, fitnesstrait="Phi", phenotraits=c("Mass", "BL", "TL"), family="binomial")
#'
#' @export
Fiz <- function(dat, fitnesstrait, phenotraits, yearname="Year", MinDataPoints=5, covariates=NULL, NAsToMean=TRUE, family="gaussian"){
  
  #### Filter out fitness NA ####
  dat <- dat[!is.na(dat[,fitnesstrait]) & !is.nan(dat[,fitnesstrait]),]
  
  #### Checks and warnings ####
  if(family=="poisson" | family=="quasipoisson")
    if(!is.integer(dat[,fitnesstrait]))
    {
      print("Family is (quasi-)Poisson, but fitness data contains non-integer... maybe that's okay, just saying!")
    }
  if(family=="binomial")
    if(mean(unique(dat[,fitnesstrait]) %in% c(0,1))<1)
    {
      print("Family is binomial, but fitness data contains non-binary!")
    }
  
  if(!is.numeric(dat[,yearname])) # check we don't get nonsense from a factor or something
  {
    stop("The column yearname is not numeric.")
  }
  
  #### Create year labels
  years <- unique(dat[,yearname])# extract years
  years <- years[order(years)]# re-order them
  years <- years[-length(years)]# remove the last year (for which we won't have population growth rate)
  
  #### Define containers ####
  # R-square explained by all traits
  sq <- vector(length = length(unique(dat[,yearname]))-1, ) # assumes we discard the last year because we don't have full fitness data yet, and no population growth rate to the next year
  names(sq) <- years
  
  # i_z for traits (partial opportunity for selection)
  iz <- matrix(data = NA, nrow = length(phenotraits), ncol = length(unique(dat[,yearname]))-1, dimnames = list(phenotraits, years))
  
  # total opportunity for selection
  Itot <- vector(length=length(unique(dat[,yearname]))-1)
  names(Itot) <- years
  
  #### Loop calculations on years ####
  loopyear <- min(dat[,yearname]):(max(dat[,yearname])-1)
  for (y in loopyear)
  {
    ThisYear <- dat[dat[,yearname]==y,]
    if(length(phenotraits)==1)
    {
      checktraitdat <- sum(!is.na(ThisYear[,c(phenotraits)]))>MinDataPoints
      names(checktraitdat) <- phenotraits
    }else{
      checktraitdat <- apply(ThisYear[,c(phenotraits)], MARGIN = 2,FUN = function(x){sum(!is.na(x))})>MinDataPoints
    }
    NbGoodTraits <- sum(checktraitdat)
    if( NbGoodTraits>0)
    {
      if(NAsToMean)#replace phenotypic NAs by year mean
      {
        for (trait in 1:NbGoodTraits)
        {
          focaltrait <- names(checktraitdat)[checktraitdat][trait]
          ThisYear[is.na(ThisYear[,focaltrait]),focaltrait] <- mean(ThisYear[,focaltrait], na.rm = TRUE)
        }
      }
      if (!is.null(covariates))#standardize covariates
      {
        for (covar in 1:length(covariates))
        {
          if(length(table( ThisYear[,covariates[covar]] ))>1){
          ThisYear[,covariates[covar]] <- as.numeric(ThisYear[,covariates[covar]])
          ThisYear[is.na(ThisYear[,covariates[covar]]),covariates[covar]]<- mean(ThisYear[,covariates[covar]], na.rm = TRUE)
          ThisYear[,covariates[covar]] <- (ThisYear[,covariates[covar]]-mean(ThisYear[,covariates[covar]]))/stats::sd(ThisYear[,covariates[covar]])
          }
        }
      }
      if(NbGoodTraits==1)
      {
        VZ <- stats::var(ThisYear[,c(phenotraits)], na.rm=TRUE)
      }else{
        VZ <- apply(ThisYear[,c(phenotraits)], MARGIN = 2,FUN = function(x){stats::var(x, na.rm = TRUE)})
      }
      for (trait in 1:NbGoodTraits)#standardization, not sure if it is useful now, but thought could be used to compute partial R2
      {
        focaltrait <- names(checktraitdat)[checktraitdat][trait]
        ThisYear[,focaltrait] <- (ThisYear[,focaltrait] - mean(ThisYear[,focaltrait], na.rm = TRUE))/stats::sd(ThisYear[,focaltrait], na.rm = TRUE)
      }
      
      if(NbGoodTraits==1)
      {
        VZ <- stats::var(ThisYear[,c(phenotraits)], na.rm=TRUE)
      }else{
        VZ <- apply(ThisYear[,c(phenotraits)], MARGIN = 2,FUN = function(x){stats::var(x, na.rm = TRUE)})#Should be 1 or NA
      }
      
      fitness <- ThisYear[,fitnesstrait]
      #the multiple regression
      
      if(is.null(covariates))
      {
        if(NbGoodTraits==1)
        {
          mr <- stats::glm(fitness ~ ThisYear[,c(phenotraits[checktraitdat])],
                    family=family)
        }else{#if more than 1 trait
          mr <- stats::glm(fitness ~ ., data=ThisYear[,c(phenotraits[checktraitdat])],
                    family=family)
        }#end if/else(NbGoodTraits==1)
        
        smr <- summary(mr)
        R2covariates <- 0
        R2pheno <- (smr$null.deviance-smr$deviance)/smr$null.deviance
        Itot[y-min(dat[,yearname])+1] <- stats::var(fitness/mean(fitness))
      }else{#if covariates
        mr <- stats::glm(fitness ~ ., data=ThisYear[,c(phenotraits[checktraitdat], covariates)],
                  family=family)
        smr <- summary(mr)#not used, but safer
        if(length(covariates)==1)
        {
          mCovariates <- stats::glm(fitness ~ ThisYear[,c(covariates)], family=family)
        }else{
          mCovariates <- stats::glm(fitness ~ ., data=ThisYear[,c(covariates)], family=family)
        }
        smCovariates <- summary(mCovariates)
        R2covariates <- (smCovariates$null.deviance-smCovariates$deviance)/smCovariates$null.deviance
        
        Itot[ y - min( dat[,yearname] ) + 1 ] <- stats::var( ( as.matrix( ThisYear [,phenotraits[checktraitdat]] ) %*% as.vector( mr$coefficients[2:(NbGoodTraits+1)])) + mr$residuals )/ (mean(ThisYear[,fitnesstrait], na.rm=TRUE)^2)
        
        R2pheno <- stats::var( as.matrix( ThisYear[,phenotraits[checktraitdat]] ) %*% as.vector( mr$coefficients[2:(NbGoodTraits+1)]) ) / stats::var( ( as.matrix( ThisYear [,phenotraits[checktraitdat]] ) %*% as.vector( mr$coefficients[2:(NbGoodTraits+1)])) + mr$residuals )
        
      }#end if/else(is.null(covariates))
      
      
      betaWZ <- checktraitdat
      betaWZ[1:length(phenotraits)] <- NA
      betaWZ[which(checktraitdat)] <- stats::coefficients(mr)[2:(NbGoodTraits+1)]/mean(ThisYear[,fitnesstrait], na.rm=TRUE)
      
      bWZ <- checktraitdat
      bWZ[1:length(phenotraits)] <- NA
      for (trait in 1:NbGoodTraits)
      {
        focaltrait <- names(checktraitdat)[checktraitdat][trait]
        if(is.null(covariates))
        {
          bWZ[which(checktraitdat)[trait] ] <- stats::coefficients(stats::glm(fitness ~ThisYear[,c(focaltrait)], data=))[2]/mean(ThisYear[,fitnesstrait], na.rm=TRUE)
        }else{
          bWZ[which(checktraitdat)[trait] ] <- stats::coefficients(stats::glm(fitness ~., data=ThisYear[,c(focaltrait, covariates)]))[2]/mean(ThisYear[,fitnesstrait], na.rm=TRUE)
        }
      }#end for (trait in 1:NbGoodTraits)
      
      
      iz[,y-min(dat[,yearname])+1] <- betaWZ*bWZ*VZ
      sq[y-min(dat[,yearname])+1] <- R2pheno
      
    }#end if NbGoodTraits>0
  }#end for years

  SampleSizes <- table(dat[,yearname])
  SampleSizes <- SampleSizes[-length(SampleSizes)]
  
  return(list(sq=sq, iz=iz, Itot= Itot, years = years, SampleSizes = SampleSizes))

}#end function
