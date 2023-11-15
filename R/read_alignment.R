#' Read MSDIAL alignment file
#' @param path Path to mass spectrometry alignment file.
#' @param format The format of the provided alignment file. Currently, only
#' MS-DIAL '.txt' files are supported (\code{msdial}).
#' @importFrom utils read.csv
#' @return Returns \code{msdial_alignment} object. A list of 3 data.frames,
#' containing peak data (\code{tab}), peak metadata (\code{peak_meta}) and
#' sample metadata (\code{sample_meta}).
#' @export
#' @author Ethan Bass

ms_read_alignment <- function(path, format = c("msdial")){
  format <- match.arg(format, c("msdial"))
  x <- read.csv(path, sep = "\t", skip = 4, header = TRUE)
  colnames(x)[(ncol(x)-1):ncol(x)] <- c("MEAN", "SD")

  # get peak metadata
  meta.idx <- c(1:grep("EI.spectrum", colnames(x)), (ncol(x)-1):ncol(x))
  peak_meta <- x[,meta.idx]
  rownames(peak_meta) <- paste0("V", rownames(x))
  if (all(peak_meta$Average.RI==-1)){
    peak_meta$Average.RI <- NA
  }

  # get peak table
  tab <- as.data.frame(t(x[,-meta.idx]))

  structure(.Data = list(tab = tab, peak_meta = peak_meta,
                         sample_meta = data.frame(full.name = rownames(tab))),
            class = "msdial_alignment")
}


#' Filter alignment by provided indices.
#' @param x An \code{msdial_alignment} object or matrix with rows as samples and features as columns.
#' @param idx Indices to be retained or excluded according to the value of \code{inverse}.
#' @param what Which dimension to filter on. Either (\code{rows}) or columns
#' (\code{cols}).
#' @param inverse Whether to retain (default) or remove the specified columns.
#' @author Ethan Bass
#' @export
ms_filter_alignment <- function(x, idx, what=c("rows","cols"), inverse = FALSE){
  what <- match.arg(what, c("rows", "cols"))
  if (inverse){
    idx <- -idx
  }
  if (what == "rows"){
    x$tab <- x$tab[idx, ]
    x$sample_meta <- x$sample_meta[idx,]
  } else if (what == "cols"){
    x$tab <- x$tab[,idx]
    x$peak_meta <- x$peak_meta[idx,]
    if (!is.null(x$matches)){
      x$matches <- x$matches[idx]
    }
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

