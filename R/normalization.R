#' @title Probabilistic Quotient Normalization
#'
#' @description Performs Probabilistic Quotient Normalization
#'
#' @param x An \code{msdial_alignment} object or matrix with rows as samples and features as columns.
#' @param n normalization reference: "mean" for using the overall average of variables as reference
#' or "median" (default) for using the overall median of variables as reference
#' @param QC vector of number(s) to specify samples which average to use as reference
#' (e.g. QC samples)
#'
#' @return A normalized \code{msdial_alignment} object or matrix, according to the input.
#'
#' @importFrom stats median
#'
#' @author E. Nevedomskaya
#' @author Rico Derks
#' @author Ethan Bass
#' @note Adapted from Rcpm package (https://github.com/ricoderks/Rcpm)
#' @references Dieterle, F., Ross, A., Schlotterbeck, G. & Senn, H. Probabilistic Quotient
#' Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures.
#' Application in H1 NMR Metabonomics. Anal. Chem. 78, 4281-4290 (2006).
#' @export

pqn <- function(x, n = "median", QC = NULL) {
  if (inherits(x, what = "msdial_alignment")){
    X <- x$tab
  } else if (class(x) %in% c("data.frame","matrix")){
    X <- x
  }
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)

  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }

  # do the actual normalization
  for (i in 1:nrow(X)) {
    X.norm[i, ] <- as.numeric(X[i, ] / median(as.numeric(X[i, ] / mX), na.rm=TRUE))
  }
  if (inherits(x, what = "msdial_alignment")){
    x$tab <- X.norm
  } else{x <- X.norm}
  x
}

#' Subtract blanks
#' @param x A \code{msdial_alignment} object.
#' @param blanks.idx Indices of blank samples
#' @param blanks.pattern A string that uniquely identifies blank samples by name
#' @param what Whether to subtract the mean or median value
#' @param drop Logical. Whether to drop columns containing only zeros. Defaults to TRUE.
#' @return A \code{msdial_alignment} object with the mean or median of the blanks
#' subtracted from each peak.

subtract_blanks <- function(x, blanks.idx, blanks.pattern,
                            what=c("mean","median"), drop = TRUE){
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
  x.n <- apply(x.n, c(1,2), function(y) max(y,0))

  # drop 0 columns
  if (drop){
    zeros <- which(colMeans(x.n) == 0)
    if (length(zeros) > 0){
      x.n <- x.n[,-zeros]
    }
  }

  x$tab <- x.n
  x
}

