#' Read MSDIAL alignment file
#' @param path Path to MS dial alignment file
#' @importFrom utils read.csv
#' @return Returns \code{msdial_alignment} object containing peak data and metadata.
#' @export
#' @author Ethan Bass

read_alignment <- function(path){
  x<-read.csv(path, sep="\t",skip = 4, header = TRUE)
  colnames(x)[(ncol(x)-1):ncol(x)] <- c("MEAN", "SD")

  # invert and convert to dataframe
  x1 <- as.data.frame(t(x))
  rownames(x1) <- colnames(x)

  # trim peak metadata
  meta.idx <- 1:28
  pk_meta <- as.data.frame(x1[meta.idx,])

  x2 <- x1[-meta.idx,]
  x2 <- apply(x2, 2, as.numeric)
  x2 <- as.data.frame(x2)
  rownames(x2) <- rownames(x1)[-meta.idx]

  structure(.Data = list(peaks = x2, peak_meta = x1, sample_meta = NA),
            class="msdial_alignment")
}
