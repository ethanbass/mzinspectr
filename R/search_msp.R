#' Get spectrum from MSDIAL alignment object
#' @param x An \code{msdial_alignment} object or matrix with rows as samples and features as columns.
#' @param col Index of the feature (column).
#' @author Ethan Bass
#' @export
get_spectrum <- function(x, col){
  spec <- tidy_eispectrum(x$peak_meta[col, "EI.spectrum"])
  spec
}

#' Search MSP database for spectrum
#' @param x Spectrum, as produced by \code{\link{get_spectrum}}.
#' @param db MSP database as list
#' @param n.results Number of results to return
#' @param parallel Logical. Whether to use parallel processing. (This feature
#' does not work on Windows).
#' @param mc.cores How many cores to use for parallel processing? Defaults to 2.
#' @param ... Additional arguments to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}
#' @importFrom pbapply pblapply
#' @importFrom OrgMassSpecR SpectrumSimilarity
#' @author Ethan Bass
#' @export

search_msp <- function(x, db, n.results = 10, parallel, mc.cores = 2, ...){
  if (missing(parallel)){
    parallel <- .Platform$OS.type != "windows"
  } else if (parallel & .Platform$OS.type == "windows"){
    parallel <- FALSE
    warning("Parallel processing is not currently available on Windows.")
  }
  sim <- as.numeric(unlist(pblapply(seq_along(db), function(i){
    SpectrumSimilarity(spec.top = x, spec.bottom = db[[i]]$Spectra,
                       print.graphic = FALSE, ...)
  }, cl = mc.cores)))
  results <- do.call(rbind, db[order(sim, decreasing = TRUE)[seq_len(n.results)]])
  r <- as.data.frame(results[,which(colnames(results)!="Spectra")])
  r$SS <- sim[order(sim, decreasing = TRUE)[seq_len(n.results)]]
  r
}
