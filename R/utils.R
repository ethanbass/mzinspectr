
#' Tidy EI spectrum
#' @description Converts MS-DIAL EI spectrum to data.frame
#' @param x spectrum to be converted
#' @importFrom stringr str_split_fixed
#' @return Returns EI spectrum as a \code{data.frame}.
#' @author Ethan Bass

tidy_eispectrum <- function(x){
  x<-str_split_fixed(strsplit(x,split = " ")[[1]], ":", 2)
  x <- apply(x,2,as.numeric)
  colnames(x) <- c("mz", "intensity")
  as.data.frame(x)
}

