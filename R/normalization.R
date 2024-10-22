#' @title Probabilistic Quotient Normalization
#'
#' @description Performs Probabilistic Quotient Normalization on peak table.
#'
#' @param x A \code{ms_alignment} object or matrix with rows as samples and
#' features as columns.
#' @param ref Reference for normalization: either \code{median} (default) to use
#' the overall median of variables as the reference, or \code{mean} to use the
#' overall average of variables as the reference.
#' @param QC vector of number(s) to specify samples which average to use as reference
#' (e.g. QC samples)
#'
#' @return A normalized \code{ms_alignment} object or \code{matrix},
#' according to the input.
#'
#' @importFrom stats median
#'
#' @author E. Nevedomskaya
#' @author Rico Derks
#' @author Ethan Bass
#' @note Adapted from the \href{https://github.com/ricoderks/Rcpm}{Rcpm} package
#' by  Rico Derks (licensed under GPL3).
#' @references Dieterle, F., Ross, A., Schlotterbeck, G. & Senn, H. Probabilistic Quotient
#' Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures.
#' Application in H1 NMR Metabonomics. Anal. Chem. 78, 4281-4290 (2006).
#' @export

ms_normalize_pqn <- function(x, ref = c("median", "mean"), QC = NULL) {
  ref <- match.arg(ref, c("median", "mean"))
  if (inherits(x, what = "ms_alignment")){
    X <- x$tab
  } else if (class(x) %in% c("data.frame","matrix")){
    X <- x
  }
  if (!is.null(QC)){
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (ref == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (ref == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (ref == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (ref == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }

  # do the actual normalization
  X.norm <- t(apply(X, 1, function(Xi){
    Xi / median(as.numeric(Xi / mX), na.rm=TRUE)
  }))

  X.norm <- as.data.frame(X.norm)

  if (inherits(x, what = "ms_alignment")){
    x$tab <- X.norm
  } else{
    x <- X.norm
  }
  x
}

#' Total sum normalization
#'
#' Divides each row by the sum of the features in that row.
#'
#' @param x An \code{ms_alignment} object or matrix with rows as samples and features as columns.
#' @return A normalized \code{ms_alignment} object or \code{matrix},
#' according to the input.
#' @export

ms_normalize_tsn <- function(x) {
  if (inherits(x, what = "ms_alignment")){
    X <- x$tab
  } else if (class(x) %in% c("data.frame","matrix")){
    X <- x
  }

  # do the actual normalization
  X.norm <- t(apply(X, 1, function(Xi){
    Xi / sum(Xi, na.rm=TRUE)
  }))

  X.norm <- as.data.frame(X.norm)

  if (inherits(x, what = "ms_alignment")){
    x$tab <- X.norm
  } else{
    x <- X.norm
  }
  x
}

#' Internal standard normalization
#'
#' Normalize by internal standard.
#'
#' @param x An \code{ms_alignment} object or matrix with rows as samples and
#' features as columns.
#' @param idx Column index of internal standard.
#' @param plot_it Logical. Whether to plot ITSD against total peak area.
#' @importFrom graphics abline legend plot
#' @importFrom stats lm
#' @author Ethan Bass
#' @return A normalized \code{ms_alignment} object or \code{matrix},
#' according to the input.
#' @export
ms_normalize_itsd <- function(x, idx, plot_it = FALSE) {
  if (inherits(x, what = "ms_alignment")){
    X <- x$tab
  } else if (class(x) %in% c("data.frame","matrix")){
    X <- x
  }

  # do the actual normalization
  zeros <- which(X[,idx] == 0)
  if (length(zeros)>0){
    warning(paste("Internal standard is 0 in the following samples:",
                  paste(sQuote(row.names(x)[zeros]), collapse=", ")))
  }

  if (plot_it){
    m <- lm(rowSums(X) ~ X[,idx])
    plot(rowSums(X) ~ X[, idx], xlab = "[ITSD]", ylab = "total peak area")
    abline(m)
    r2 <- format(summary(m)$r.squared, digits = 3)
    legend("top",
           legend = bquote(R^2 == .(r2)),
           bty = "n")
  }

  X.norm <- t(apply(X, 1, function(Xi){
    Xi / Xi[idx]
  }))

  X.norm <- as.data.frame(X.norm)

  if (inherits(x, what = "ms_alignment")){
    x$tab <- X.norm
  } else{
    x <- X.norm
  }
  x
}

#' Subtract blanks
#' @param x A \code{ms_alignment} object.
#' @param blanks.idx Indices of blank samples
#' @param blanks.pattern A string that uniquely identifies blank samples by name
#' @param what Whether to subtract the mean or median value
#' @param drop Logical. Whether to drop columns containing only zeros. Defaults to TRUE.
#' @return A \code{ms_alignment} object with the mean or median of the blanks
#' subtracted from each peak.
#' @export

ms_subtract_blanks <- function(x, blanks.idx, blanks.pattern,
                            what=c("mean", "median"), drop = TRUE){
  if (missing(blanks.idx)){
    if (!missing(blanks.pattern)){
      blanks.idx <- grep(blanks.pattern, x$sample_meta$full.name)
      if (length(blanks.idx) == 0){
        stop("Blanks not found! Please double check the supplied pattern.")
      }
    } else{
      stop("Must define blanks by providing an index (`blanks.idx`) or pattern (`blanks.pattern`).")
    }
  }
  # For each column, subtract mean or median of blanks.
  what <- match.arg(what, c("mean","median"))
  fn <- switch(what, "mean" = mean,
               "median" = median)

  x.n <- apply(x$tab, 2, function(col){
    col - fn(col[blanks.idx])
  })

  # Round any negative nunbers up to 0
  x.n <- apply(x.n, c(1, 2), function(y) max(y,0))

  # drop 0 columns
  if (drop){
    zeros <- which(colMeans(x.n) == 0)
    if (length(zeros) > 0){
      x.n <- x.n[, -zeros]
    }
  }
  x.n <- as.data.frame(x.n)
  x$tab <- x.n
  x
}
