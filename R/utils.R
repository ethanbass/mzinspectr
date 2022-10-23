
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

#' Find peak based on retention time and/or mass
#' @param x An \code{msdial_alignment} object
#' @param rt Retention time
#' @param mz Quant.mass
#' @param rt.tol Tolerance for matching retention time
#' @param mz.tol Tolerance for matching Quant.mass
#' @importFrom stringr str_split_fixed
#' @return Returns EI spectrum as a \code{data.frame}.
#' @author Ethan Bass

find_peak <- function(x, rt, mz, rt.tol = .01, mz.tol =.05){
  pks <- seq_along(x$peak_meta)
  if (!missing(rt)){
    pks1 <- which(abs(rt-as.numeric(x$peak_meta["Average.Rt.min.",])) < rt.tol)
  } else{pks1 <- pks}
  if (!missing(mz)){
    pks2 <- which(abs(mz - as.numeric(x$peak_meta["Quant.mass",])) < mz.tol)
  } else{pks2 <- pks}
  intersect(pks1, pks2)
}
