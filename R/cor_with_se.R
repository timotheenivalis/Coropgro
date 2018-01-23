#' Correlation with standard error
#'
#' Calculates correlation between x and y, possibly with weights, with standard error.
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @param wt Numeric vector or NULL. If not NULL, contains weights.
#' 
#' @return A vector of length 2, containing the estimated correlation between x and y, followed by the standard error of the estimate.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100) + x
#' cor_with_se(x,y)
#'
#' @export
cor_with_se<-function(x,y, wt=NULL){
  if(is.null(wt))
    {
    r <- stats::cor(x,y, use="complete.obs")
  }else{
  r <- stats::cov.wt(cbind(x,y), wt = as.vector(wt), method = "unbiased", cor = TRUE)$cor[1,2]
  }
  
  n<-length(x)
  SEr<-sqrt( (1-r^2) / (n-2))
  df<-c(r, SEr)
  names(df)<-c("correlation","standard_error")
  return(df)
}