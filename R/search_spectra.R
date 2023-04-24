#' Search spectra in MSDIAL alignment against database
#'
#' This function can be used to identify peaks in your peaktable by matching them
#' to a spectral database (\code{db}). It takes several arguments that can
#' be used to customize the matching algorithm, including \code{ri_thresh},
#' \code{spectral weight}, \code{n_results}. The retention index threshold
#' (\code{ri_thresh}) is used to subset the provided database, which greatly
#' improves the search speed. Only database entries with a retention index
#' falling within the specified threshold will be considered. The spectral
#' weight affects the weight given to spectral similarity when calculating the
#' the total similarity score, which is used to rank matches.
#'
#' @param x An \code{msdial_alignment} object.
#' @param db MSP database. The provided object should be a nested list, where the
#' sublists contain the following elements: retention indices in an element named
#' \code{RI} and mass spectra in an element called \code{Spectra}. All other elements
#' are optional.
#' @param cols Index or indices of feature(s) to be identified.
#' @param ... Additional arguments to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}.
#' @param ri_thresh Maximum difference between retention indices for a match.
#' to be considered. Defaults to 100.
#' @param spectral_weight A number between 0 and 1 specifying the weight given.
#' to spectral similarity versus retention index similarity. Defaults to 0.6.
#' @param n_results How many results to return. Defaults to 10.
#' @param parallel Logical. Whether to use parallel processing. (This feature
#' does not work on Windows).
#' @param mc.cores How many cores to use for parallel processing? Defaults to 2.
#' @param print Logical. Whether to print the results after each search. Defaults
#' to FALSE.
#' @param ris Retention indices to use.
#' @param progress_bar Logical. Whether to display  progress bar or not.
#' @note See \href{https://github.com/QizhiSu/mspcompiler}{mspcompiler} for help compiling
#' an msp database.
#' @return Returns a modified \code{msdial_alignment} object with database matches
#' in the \code{matches} slot as a list of data frames. Each \code{data.frame}
#' will contain the database matches as rows and columns corresponding to the
#' elements of the database entry (e.g. "Name", "InChIKey", etc.) as well as
#' match scores for spectral similarity (\code{spectral_match}), retention index
#' similarity (\code{ri_match}) and the total similarity score (\code{total_score}).
#' @author Ethan Bass
#' @export

ms_search_spectra <- function(x, db, cols, ..., ri_thresh = 100, spectral_weight = 0.6,
                           n_results=10, parallel, mc.cores = 2,  print = FALSE,
                           ris, progress_bar=TRUE){
  if (is.null(x$matches)){
    x$matches <- as.list(rep(NA, ncol(x$tab)))
    names(x$matches) <- colnames(x$tab)
  }
  if (missing(ris)){
    ris <- sapply(db, function(x) x$RI)
  }
  if (missing(cols)){
    cols <- seq_len(ncol(x$tab))
  }
  if (missing(parallel)){
    parallel <- .Platform$OS.type != "windows"
  }
  laplee <- laplee <- choose_apply_fnc(progress_bar = progress_bar, cl = mc.cores)
  x$matches[cols] <- laplee(cols, function(col){
    try({
      sp <- ms_get_spectrum(x, col)
      ri_diff <- abs(as.numeric(x$peak_meta[col, "Average.RI"]) - ris)
      idx <- which(ri_diff < ri_thresh)
      sp_score <- search_msp(sp, db[idx], ..., what="scores", parallel = FALSE)
      ri_score <- ri_diff[idx]/ri_thresh
      total_score <- sp_score*spectral_weight + ri_score*(1 - spectral_weight)
      sel <- order(total_score, decreasing = TRUE)[seq_len(n_results)]
      results <- msp_to_dataframe(db[idx][sel])
      results$spectral_match <- sp_score[sel]
      results$ri_match <- ri_score[sel]
      results$total_score <- total_score[sel]
      results
    })
  })
  x
}

#' Get spectrum from MSDIAL alignment object
#' @param x An \code{msdial_alignment} object or matrix with rows as samples and features as columns.
#' @param col Index of the feature (column).
#' @return Returns spectrum as a data.frame with two columns: "mz" and "intensity".
#' @author Ethan Bass
#' @export

ms_get_spectrum <- function(x, col){
  spec <- tidy_eispectrum(x$peak_meta[col, "EI.spectrum"])
  spec
}

#' Search MSP database for spectrum
#' @param x Spectrum, as produced by \code{\link{get_spectrum}}.
#' @param db MSP database as list
#' @param ... Additional arguments to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}
#' @param n_results Number of results to return
#' @param parallel Logical. Whether to use parallel processing. (This feature
#' does not work on Windows).
#' @param mc.cores How many cores to use for parallel processing? Defaults to 2.
#' @param what What kind of object to return. Either \code{msdial_alignment} object,
#'  (\code{msd}), or \code{data.frame} (\code{df}).
#' @param progress_bar Logical. Whether to display progress bar or not.
#' @importFrom pbapply pblapply
#' @importFrom OrgMassSpecR SpectrumSimilarity
#' @author Ethan Bass
#' @noRd

search_msp <- function(x, db, ..., n_results = 10, parallel, mc.cores = 2,
                       what=c("msd", "df", "scores"), progress_bar = FALSE){
  what <- match.arg(what, c("msd", "df", "scores"))
  if (missing(parallel)){
    parallel <- .Platform$OS.type != "windows"
  } else if (parallel & .Platform$OS.type == "windows"){
    parallel <- FALSE
    warning("Parallel processing is not currently available on Windows.")
  }
  laplee <- laplee <- choose_apply_fnc(progress_bar = progress_bar, cl = mc.cores)
  sim <- unlist(laplee(seq_along(db), function(i){
    db[[i]]$Spectra <- as.data.frame(apply(db[[i]]$Spectra, 2, as.numeric))
    try(spectral_similarity(spec.top = x, spec.bottom = db[[i]]$Spectra, ...))
  }))
  sim <- suppressWarnings(as.numeric(sim))
  if (what == "scores"){
    r <- sim
  } else {
    results <- do.call(rbind, db[order(sim, decreasing = TRUE)[seq_len(n_results)]])
    r <- as.data.frame(results[,which(colnames(results)!="Spectra")])
    r$SS <- sim[order(sim, decreasing = TRUE)[seq_len(n_results)]]
  }
  r
}

#' Convert MSP to data.frame
#' @noRd
msp_to_dataframe <- function(db){
  x <- lapply(db, function(xx){
    xx[sapply(xx,length) == 0] <- NA
    xx$Spectra <- condense_eispectrum(xx$Spectra)
    as.data.frame(xx)
    # as.data.frame(xx[-c(which(names(xx) == "Spectra"))])
  })
  x<-as.data.frame(do.call(rbind,x))
  x
}

#' Calculate spectral similarity between two peaks
#' @noRd
spectral_similarity <- function(spec.top, spec.bottom, t = 0.25, b = 10,
                                xlim = c(50, 1200), x.threshold = 0){
  colnames(spec.top) <- c("mz", "intensity")
  spec.top$normalized <- round((spec.top$intensity/max(spec.top$intensity)) * 100)
  spec.top <- spec.top[which(spec.top$mz >= xlim[1] & spec.top$mz <= xlim[2]),]
  top <- spec.top[which(spec.top$normalized >= b),]

  colnames(spec.bottom) <- c("mz", "intensity")
  spec.bottom$normalized <- round((spec.bottom$intensity/max(spec.bottom$intensity)) * 100)
  spec.bottom <- spec.bottom[which(spec.bottom$mz >= xlim[1] & spec.bottom$mz <= xlim[2]),]
  bottom <- spec.bottom[which(spec.bottom$normalized >= b),]

  for (i in 1:nrow(bottom)){
    top[, 1][which(bottom[, 1][i] >= top[,1] - t & bottom[, 1][i] <= top[, 1] + t)] <- bottom[,1][i]
  }

  alignment <- merge(top[,-2], bottom[,-2], by = 1, all = TRUE)
  if (length(unique(alignment[, 1])) != length(alignment[, 1]))
    warning("the m/z tolerance is set too high")
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")

  if (x.threshold < 0)
    stop("x.threshold argument must be zero or a positive number")
  alignment <- alignment[alignment[, 1] >= x.threshold, ]
  u <- alignment[, 2]
  v <- alignment[, 3]
  similarity_score <- sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
  return(similarity_score)
}
