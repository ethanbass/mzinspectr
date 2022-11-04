
#' Search spectra in MSDIAL alignment against database
#' @param x Spectrum, as produced by \code{\link{get_spectrum}}.
#' @param db MSP database as list
#' @param cols Index or indices of feature(s) to select
#' @param ... Additional arguments to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}
#' @param ri_thresh Maximum difference between retention indices for a match
#' to be considered.
#' @param spectral_weight A number between 0 and 1 to specify the weight given
#' to spectral similarity versus retention index similarity.
#' @param n.results How many results to return.
#' @param mc.cores How many cores to use for parallel processing? Defaults to 2.
#' @param ris Retention indices to use
#' @author Ethan Bass
#' @export

search_spectra <- function(x, db, cols, ..., ri_thresh = 100, spectral_weight = 0.6,
                           n.results=10, mc.cores = 2,  ris){
  if (any(is.null(x$matches))){
    x$matches <- as.list(rep(NA, ncol(x$tab)))
    names(x$matches) <- colnames(x$tab)
  }
  if (missing(ris)){
    ris <- sapply(db, function(x) x$RI)
  }
  for (col in cols){
    sp <- get_spectrum(x, col)
    ri_diff <- abs(as.numeric(x$peak_meta[col, "Average.RI"]) - ris)
    idx <- which(ri_diff < ri_thresh)
    sp_score <- search_msp(sp, db[idx], ..., what="scores", mc.cores = mc.cores)
    ri_score <- ri_diff[idx]/ri_thresh
    total_score <- sp_score*spectral_weight + ri_score*(1-spectral_weight)
    sel <- order(total_score, decreasing=TRUE)[1:n.results]
    results <- msp_to_dataframe(db[idx][sel])
    results$spectral_match <- sp_score[sel]
    results$ri_match <- ri_score[sel]
    results$total_score <- total_score[sel]
    x$matches[[col]] <- results
  }
  x
}

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
#' @param ... Additional arguments to \code{\link[OrgMassSpecR]{SpectrumSimilarity}}
#' @param n.results Number of results to return
#' @param parallel Logical. Whether to use parallel processing. (This feature
#' does not work on Windows).
#' @param mc.cores How many cores to use for parallel processing? Defaults to 2.
#' @param what What kind of object to return. Either \code{msdial_alignment} object (\code{msd})
#' or \code{data.frame} (\code{df}).
#' @importFrom pbapply pblapply
#' @importFrom OrgMassSpecR SpectrumSimilarity
#' @author Ethan Bass
#' @export

search_msp <- function(x, db, ..., n.results = 10, parallel, mc.cores = 2,
                       what=c("msd", "df","scores")){
  # result <- match.arg(result, c("msd", "df", "scores"))
  if (missing(parallel)){
    parallel <- .Platform$OS.type != "windows"
  } else if (parallel & .Platform$OS.type == "windows"){
    parallel <- FALSE
    warning("Parallel processing is not currently available on Windows.")
  }
  sim <- as.numeric(unlist(pblapply(seq_along(db), function(i){
    try(SpectrumSimilarity(spec.top = x, spec.bottom = db[[i]]$Spectra,
                       print.graphic = FALSE, ...))
  }, cl = mc.cores)))
  if (what == "scores"){
    r <- sim
  } else{
  results <- do.call(rbind, db[order(sim, decreasing = TRUE)[seq_len(n.results)]])
  r <- as.data.frame(results[,which(colnames(results)!="Spectra")])
  r$SS <- sim[order(sim, decreasing = TRUE)[seq_len(n.results)]]
  }
  r
}

#'@noRd
msp_to_dataframe <- function(db){
  m <- do.call(rbind,db)
  as.data.frame(m[,-which(colnames(m)=="Spectra")])
}
