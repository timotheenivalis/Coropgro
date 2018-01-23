#' Inserts individuals missing in years they were alive but not measured
#'
#' Find individuals that are present in multiple years, but missing from the phenotypic data when they were known to be alive. Gives those individuals data for survival (1, as they survived) and missing data for fecundity.
#' 
#' @param data - data frame containing individual data on fitness
#' @param idname The name of the column containing individual identification, defaults to "ID"
#' @param yearname The name of the column containing year identity (or other time steps), defaults to "Year"
#' @param survivalname The name of the column containing survival information, defaults to "survival"
#' @param verbose Boolean. If TRUE, print the number of missing entries added
#' 
#' @return A data frame containing all the original data, with missing data added in. If no missing data was found, the original dataframe is returned.
#' 
#' 
#' @examples
#' data("svdata")
#' svdata <- InsertMissing(data = svdata, idname = "ID", yearname = "Year", 
#' survivalname = "Phi", verbose = TRUE)
#' 
#' @export


InsertMissing<-function(data, idname="ID", yearname="Year",
                        survivalname="survival", 
                        verbose=FALSE){
  
  data$IDyear<-paste(data[,idname], data[,yearname], sep=".")
  
    yearmin<-tapply(X = data[,yearname], INDEX = data[,idname], function(x) {
    min(x)})
  yearmax<-tapply(X = data[,yearname], INDEX = data[,idname], function(x) {
    max(x)})
  
  yearframe<-data.frame(ID=unique(factor(data[, idname])),
                        min=NA, 
                        max=NA)
  yearframe$min<-yearmin[match(yearframe$ID, names(yearmin))]
  yearframe$max<-yearmax[match(yearframe$ID, names(yearmax))]
  
  
  allyearid<-vector()
  for(i in 1:length(yearframe$ID)){
  allyearid<-c(allyearid,paste(yearframe$ID[i], seq(yearframe$min[i], yearframe$max[i], 1), sep="."))
  }
  
  missingidyear<-allyearid[which(!allyearid%in%data$IDyear)]
  
  if(length(missingidyear)>0){
  missingdat <- data.frame(matrix(nrow = length(missingidyear), ncol = ncol(data)))
  colnames(missingdat)<-colnames(data)
  
  missingdat[,idname]<-matrix(unlist(strsplit(missingidyear, "[.]")), ncol=2, byrow=T)[,1]
  missingdat[,yearname]<-matrix(unlist(strsplit(missingidyear, "[.]")), ncol=2, byrow=T)[,2]
  missingdat[,survivalname]<-1
  
  
  #missingdat<-subset(missingdat, select= -c(IDyear) )
  missingdat <- missingdat[, ! names(missingdat) %in% ("IDyear")]
  
  #data<-subset(data, select= -c(IDyear))
  data <- data[, ! names(data) %in% ("IDyear")]
  
  if(verbose==TRUE){
    print(paste(
      missingidyear, "were found to be missing and were added"))
  }
  
  newdata<-rbind(data, missingdat)
  return(newdata)

  
   }else{
  print("No missing data found")
  #data<-subset(data, select= -c(IDyear)) #use of subset is OK but 
     #gives a Check message because IDyear is not an object devtools recognizes
  
  data <- data[,! names(data) %in% ("IDyear")]
  
  return(data)  
}

  
}# end function
  
  
  
  
  
  
  
