#' Tidy EI spectrum
#' Converts MS-DIAL EI spectrum to data.frame
#' @param x Spectrum to be converted in MSDIAL character format.
#' @importFrom stringr str_split_fixed
#' @return Returns EI spectrum as a \code{data.frame}.
#' @author Ethan Bass
#' @noRd
tidy_eispectrum <- function(x){
  x<-str_split_fixed(strsplit(x,split = " ")[[1]], ":", 2)
  x <- apply(x,2,as.numeric)
  colnames(x) <- c("mz", "intensity")
  as.data.frame(x)
}

#' Convert EI spectrum back to MSDIAL text format
#' @noRd
condense_eispectrum <- function(x){
  paste(apply(x, 1, function(xx) paste(xx,collapse = ":")), collapse=" ")
}

#' Find peak based on retention time and/or mass
#' @param x An \code{msdial_alignment} object
#' @param rt Retention time
#' @param mz Quant.mass
#' @param rt.tol Tolerance for matching retention time
#' @param mz.tol Tolerance for matching Quant.mass
#' @param plot_it Logical. Whether to plot the spectra.
#' @importFrom stringr str_split_fixed
#' @return Returns EI spectrum as a \code{data.frame}.
#' @author Ethan Bass
#' @export

ms_find_peak <- function(x, rt, mz, rt.tol = .01, mz.tol = .05, plot_it = TRUE){
  pks <- seq_len(ncol(x$tab))
  if (!missing(rt)){
    pks1 <- which(abs(rt - as.numeric(x$peak_meta[,"Average.Rt.min."])) < rt.tol)
  } else{pks1 <- pks}
  if (!missing(mz)){
    pks2 <- which(abs(mz - as.numeric(x$peak_meta[,"Quant.mass"])) < mz.tol)
  } else{pks2 <- pks}
  pks <- intersect(pks1, pks2)
  names(pks) <- colnames(x$tab)[pks]
  if (plot_it)
    for (pk in names(pks)){
      ms_plot_spectrum(x, pk, type = "base")
    }
  pks
}

#' Choose apply function
#' @importFrom parallel mclapply
#' @return Returns \code{\link[pbapply]{pblapply}} if \code{progress_bar == TRUE},
#' otherwise returns \code{\link{lapply}}.
#' @noRd
choose_apply_fnc <- function(progress_bar, parallel = FALSE, cl = NULL){
  if (progress_bar){
    check_for_pkg("pbapply")
    pblapply<-pbapply::pblapply
    fn <- purrr::partial(pblapply, cl = cl)
    # fn <- pbapply::pblapply
  } else if (!progress_bar && parallel){
    fn <- purrr::partial(mclapply, mc.cores = cl)
  } else{
    fn <- lapply
  }
  fn
}


#' Check for suggested package
#' @author Ethan Bass
#' @noRd
check_for_pkg <- function(pkg, return_boolean = FALSE){
  pkg_exists <- requireNamespace(pkg, quietly = TRUE)
  if (return_boolean){
    return(pkg_exists)
  } else if (!pkg_exists) {
    stop(paste(
      "Package", sQuote(pkg), "must be installed to perform this action:
          try", paste0("`install.packages('", pkg, "')`.")),
      call. = FALSE
    )
  }
}

#' Convert retention times to retention indices.
#' @importFrom stats approx
#' @param rts A vector of retention times.
#' @param RIs A matrix or data.frame containing retention times in column one
#' and retention indices in column two.
#' @export
ms_rt_to_ri <- function(rts, RIs){
  if (mean(nchar(RIs[,2])) == 2){
    RIs[,2] <- RIs[,2]*100
  }
  round(approx(x = as.numeric(RIs[,1]), y = as.numeric(RIs[,2]),
               xout=as.numeric(rts))$y)
}

#' Convert retention times to retention indices in alignment object.
#' @param x An \code{msdial_alignment} object.
#' @param Ris A matrix or data.frame containing retention times in column one
#' and retention indices in column two.
#' @export

ms_calculate_RIs <- function(x, Ris){
  Rts <- as.numeric(x$peak_meta$Average.Rt.min.)
  x$peak_meta$Average.RI <- ms_rt_to_ri(Rts, Ris)
  x
}
