#' Search spectra in MS alignment against database
#'
#' This function can be used to identify peaks in a peak table by matching them
#' to a spectral database (\code{db}). It takes several arguments that can
#' be used to customize the matching algorithm, including \code{ri_thresh},
#' \code{spectral weight}, \code{n_results}. The retention index threshold
#' (\code{ri_thresh}) is used to subset the provided database, which greatly
#' improves the search speed. Only database entries with a retention index
#' falling within the specified threshold will be considered. The spectral
#' weight affects the weight given to spectral similarity (versus retention
#' index similarity) when calculating the the total similarity score, which is
#' used to rank matches.
#'
#' @param x An \code{ms_alignment} object.
#' @param db MSP database. The provided object should be a nested list, where the
#' sublists contain the following elements: retention indices in an element named
#' \code{RI} and mass spectra in an element called \code{Spectra}. All other elements
#' are optional.
#' @param cols Index or indices of feature(s) to be identified.
#' @param ... Additional arguments to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}.
#' @param ri_thresh Maximum difference between retention indices for a match.
#' to be considered. Defaults to 100. Use \code{NULL} to search database without
#' filtering by retention time (this will take a long time for large databases).
#' @param spectral_weight A number between 0 and 1 specifying the weight given.
#' to spectral similarity versus retention index similarity. Defaults to 0.6.
#' @param n_results How many results to return. Defaults to 10.
#' @param parallel Logical. Whether to use parallel processing. (This feature
#' does not work on Windows).
#' @param mc.cores How many cores to use for parallel processing? Defaults to 2.
#' @param print Logical. Whether to print the results after each search. Defaults
#' to FALSE.
#' @param progress_bar Logical. Whether to display progress bar or not.
#' @note See \href{https://github.com/QizhiSu/mspcompiler}{mspcompiler} for help compiling
#' an msp database.
#' @return Returns a modified \code{ms_alignment} object with database matches
#' in the \code{matches} slot as a list of data frames. Each \code{data.frame}
#' will contain the database matches as rows and columns corresponding to the
#' elements of the database entry (e.g. "Name", "InChIKey", etc.) as well as
#' match scores for spectral similarity (\code{spectral_match}), retention index
#' similarity (\code{ri_match}) and the total similarity score (\code{total_score}).
#' @author Ethan Bass
#' @export

ms_search_spectra <- function(x, db, cols, ..., ri_thresh = 100, spectral_weight = 0.6,
                           n_results=10, parallel, mc.cores = 2,  print = FALSE,
                           progress_bar = TRUE){
  if (!is.null(ri_thresh) && all(is.na(x$peak_meta$Average.RI))){
    stop(paste("Retention indices are not present. Please add retention indices using \n \t",
               sQuote("ms_calculate_RIs"), "before proceeding or set", sQuote("ri_thresh = NULL"),"."))
  }
  if (is.null(x$matches)){
    x$matches <- as.list(rep(NA, ncol(x$tab)))
    names(x$matches) <- colnames(x$tab)
  }
  ris <- sapply(db, function(x) x$RI)
  if (missing(cols)){
    cols <- seq_len(ncol(x$tab))
  }
  if (missing(parallel)){
    parallel <- .Platform$OS.type != "windows"
  }
  laplee <- choose_apply_fnc(progress_bar = progress_bar, cl = mc.cores)
  x$matches[cols] <- laplee(cols, function(col){
    try({
      sp <- ms_get_spectrum(x, col)
      target_ri <- as.numeric(x$peak_meta[col, "Average.RI"])
      if (!is.null(ri_thresh)){
        ri_diff <- abs(target_ri - ris)
        idx <- which(ri_diff < ri_thresh)
      } else{
        idx <- seq_along(db)
      }
      # this is still slow
      sp_score <- search_msp(sp, db[idx], ..., what = "scores", parallel = FALSE)
      ri_score <- 1 - ri_diff[idx]/target_ri
      total_score <- sp_score*spectral_weight + ri_score*(1 - spectral_weight)
      sel <- order(total_score, decreasing = TRUE, method = "radix")[seq_len(n_results)]
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
#' @param x An \code{ms_alignment} object or matrix with rows as samples and features as columns.
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
#' @param what What kind of object to return. Either \code{ms_alignment} object,
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
  laplee <- choose_apply_fnc(progress_bar = progress_bar, cl = mc.cores)
  sim <- unlist(lapply(seq_along(db), function(i){
    db[[i]]$Spectra <- as.data.frame(apply(db[[i]]$Spectra, 2, as.numeric))
    try(spectral_similarity(spec.top = x, spec.bottom = db[[i]]$Spectra, ...), silent = TRUE)
  }))
  sim <- suppressWarnings(as.numeric(sim))
  if (what == "scores"){
    r <- sim
  } else {
    results <- do.call(rbind, db[order(sim, decreasing = TRUE)[seq_len(n_results)]])
    r <- as.data.frame(results[,which(colnames(results) != "Spectra")])
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
  })
  x <- as.data.frame(do.call(rbind,x))
  x
}

#' Calculate spectral similarity between two peaks
#'
#' This function is slightly adapted from the \code{SpectrumSimilarity} function
#' in [OrgMassSpecR](https://orgmassspec.github.io/) where it is licensed under
#' BSD-2 (© 2011-2017, Nathan Dodder). The function was re-factored here for
#' increased speed.
#'
#' @param spec.top data frame containing the experimental spectrum's peak list
#' with the m/z values in the first column and corresponding intensities in the
#' second.
#' @param spec.bottom data frame containing the reference spectrum's peak list
#' with the m/z values in the first column and corresponding intensities in the
#' second.
#' @param t numeric value specifying the tolerance used to align the m/z values
#' of the two spectra.
#' @param b numeric value specifying the baseline threshold for peak
#' identification. Expressed as a percent of the maximum intensity.
#' @param xlim numeric vector of length 2, defining the beginning and ending
#' values of the x-axis.
#' @param x.threshold numeric value of length 1 specifying the m/z threshold
#' used for the similarity score calculation. Only peaks with m/z values above
#' the threshold are used in the calculation. This can be used to exclude noise
#' and/or non-specific ions at the low end of the spectrum. By default all ions
#' are used.
#' @author Nathan G. Dodder
#' @author Ethan Bass

spectral_similarity <- function(spec.top, spec.bottom, t = 0.25, b = 10,
                                xlim = c(50, 1200), x.threshold = 0){
  if (x.threshold < 0){
    stop("x.threshold argument must be zero or a positive number")
  }

  top <- normalize_spectrum(spec.top, xlim = xlim, b = b)
  bottom <- normalize_spectrum(spec.bottom, xlim = xlim, b = b)

  for (i in 1:nrow(bottom)){
    top[, 1][which(bottom[, 1][i] >= top[,1] - t & bottom[, 1][i] <= top[, 1] + t)] <- bottom[,1][i]
  }

  mz <- unique(c(top$mz, bottom$mz))
  alignment <- cbind(mz, x = top[match(mz, top[, "mz"]), "normalized"],
                      y = bottom[match(mz, bottom[, "mz"]), "normalized"])

  if (length(unique(alignment[, 1])) != length(alignment[, 1]))
    warning("the m/z tolerance is set too high")
  alignment[, 3][is.na(alignment[, 3])] <- 0
  alignment[, c(2,3)][is.na(alignment[, c(2,3)])] <- 0

  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")

  alignment <- alignment[alignment[, 1] >= x.threshold, ]
  u <- alignment[, 2]
  v <- alignment[, 3]
  similarity_score <- sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
  return(similarity_score)
}

#' Normalize spectrum
#' @noRd
normalize_spectrum <- function(spec, b, xlim){
  colnames(spec) <- c("mz", "intensity")
  spec$normalized <- round((spec$intensity/max(spec$intensity)) * 100)
  spec <- spec[which(spec$mz >= xlim[1] & spec$mz <= xlim[2]),]
  spec[which(spec$normalized >= b),]
}
