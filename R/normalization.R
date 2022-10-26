#' @title Probabilistic Quotient Normalization
#'
#' @description Performs Probabilistic Quotient Normalization
#'
#' @param X matrix to normalize samples * variables (rows * columns)
#' @param n normalization reference: "mean" for using the overall average of variables as reference
#' or "median" (default) for using the overall median of variables as reference
#' @param QC vector of number(s) to specify samples which average to use as reference
#' (e.g. QC samples)
#'
#' @return Normalized table samples * variables (rows * columns)
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

pqn <- function(X, n = "median", QC = NULL) {
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

  X.norm
}


#' Subtract blanks
#' @param x A \code{msdial_alignment} object.
#' @param blanks.idx Indices of blank samples
#' @param blanks.pattern A string that uniquely identifies blank samples by name
#' @param what Whether to subtract the mean or median value
#' @return A \code{msdial_alignment} object with the mean or median of the blanks
#' subtracted from each peak.

subtract_blanks <- function(x, blanks.idx, blanks.pattern, what=c("mean","median")){
  if (missing(blanks.idx)){
    if (!missing(blanks.pattern)){
      blanks.idx <- grep(blanks.pattern, x$sample_meta$full.name)
    } else{
      stop("Must define blanks by providing an index (`blanks.idx`) or pattern (`blanks.pattern`).")
    }
  }
  # For each column, subtract mean or median of blanks.
  what <- match.arg(what, c("mean","median"))
  fn <- switch(what, "mean" = mean,
               "median" = median)
  x.n <- sapply(seq_along(x$tab), function(j){
    x$tab[,j] - fn(x$tab[blanks.idx, j])
  })

  # Round any negative nunbers up to 0
  x.n <- apply(x.n, c(1,2), function(y) max(y,0))

  # retrieve names
  dimnames(x.n) <-dimnames(x$tab)
  x$tab <- x.n
  # filter 0 columns
  x$tab <- x$tab[,-which(colMeans(x$tab)==0)]
  x
}
