#' Read MSDIAL alignment file
#' @param path Path to MS dial alignment file
#' @importFrom utils read.csv
#' @return Returns \code{msdial_alignment} object. A list of 3 data.frames,
#' containing peak data (\code{tab}), peak metadata (\code{peak_meta}) and
#' sample metadata (\code{sample_meta}).
#' @export
#' @author Ethan Bass

ms_read_alignment <- function(path){
  x<-read.csv(path, sep = "\t",skip = 4, header = TRUE)
  colnames(x)[(ncol(x)-1):ncol(x)] <- c("MEAN", "SD")

  # invert and convert to data.frame
  x1 <- as.data.frame(t(x))
  rownames(x1) <- colnames(x)

  # trim peak metadata
  meta.idx <- c(1:28, (ncol(x)-1):ncol(x))
  peak_meta <- as.data.frame(t(x1[meta.idx,]))

  x2 <- x1[-meta.idx,]
  x2 <- apply(x2, 2, as.numeric)
  x2 <- as.data.frame(x2)
  rownames(x2) <- rownames(x1)[-meta.idx]

  structure(.Data = list(tab = x2, peak_meta = peak_meta, sample_meta = data.frame(full.name = rownames(x2))),
            class="msdial_alignment")
}


#' Filter alignment by provided indices.
#' @param x An \code{msdial_alignment} object or matrix with rows as samples and features as columns.
#' @param idx Indices to be retained
#' @param what Which dimension to filter on. Either \code{rows} or columns (\code{cols}).
#' @param inverse Whether to retain (default) or remove the specified columns.
#' @author Ethan Bass
#' @export
ms_filter_alignment <- function(x, idx, what=c("rows","cols"), inverse = FALSE){
  what <- match.arg(what, c("rows","cols"))
  if (inverse){
    idx <- -idx
  }
  if (what == "rows"){
    x$tab <- x$tab[idx, ]
    x$sample_meta <- x$sample_meta[idx,]
  } else if (what == "cols"){
    x$tab <- x$tab[,idx]
    x$peak_meta <- x$peak_meta[idx,]
  }
  x
}

#' @importFrom utils head
#' @noRd
#' @export
head.msdial_alignment <- function(x,...){
  head(x$tab)
}

#' @importFrom utils tail
#' @noRd
#' @export
tail.msdial_alignment <- function(x,...){
  tail(x$tab)
}

#' @noRd
#' @export
print.msdial_alignment <- function(x, ...){
  print(x$tab)
}

#' @noRd
#' @export
dim.msdial_alignment <- function(x){
  dim(x$tab)
}

#' @noRd
#' @export
row.names.msdial_alignment <- function(x){
  row.names(x$tab)
}

