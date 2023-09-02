#' Reshape peak table
#'
#' Convert peak table to tidy format for plotting.
#'
#' @param x An MS dial alignment object.
#' @param peaks A character vector specifying the peaks to include in tidy output.
#' If the character vector is named, the names of the vector elements will be
#' used in place of the original peak names.
#' @param metadata A character vector specifying the metadata to include in
#' the tidy output.
#' @param treatments This argument is deprecated as of version 0.3.2. It is
#' synonymous with the new metadata argument which should be used instead.
#' @importFrom dplyr select mutate any_of
#' @importFrom tidyr pivot_longer
#' @return If \code{export} is \code{TRUE}, returns spectrum as \code{data.frame}.
#' Otherwise, no return value.
#' @author Ethan Bass
#' @export

ms_reshape_peaktable <- function(x, peaks, metadata, treatments = NULL){
  if (!is.null(treatments)){
    .Deprecated("metadata", package="msdialreadr", old = "treatments",
                msg="The `treatments` argument is deprecated as of msdialreadr v0.3.2.
                \t Please use the `metadata` argument instead.")
    metadata <- treatments
  }
  if (!missing(peaks)){
    if (is.numeric(peaks)){
      peaks <- colnames(df)[peaks]
    }
    x$tab <- x$tab[,which(colnames(x$tab) %in% peaks), drop = FALSE]
    if (!is.null(names(peaks))){
      colnames(x$tab) <- names(peaks)[match(colnames(x$tab), peaks)]
      peaks <- colnames(x$tab)
    }
  }
  df <- as.data.frame(x$tab)
  df <- cbind(df, select(x$sample_meta, any_of(metadata)))
  df <- mutate(df, sample = row.names(x))
  df <- pivot_longer(df, cols = peaks, names_to = "peak")
  df
}

#' Reshape peak table
#'
#' Converts peak table to tidy format for plotting. This function is deprecated
#' as of version \code{0.3.3}. Please use \code{\link{ms_reshape_peaktable}} instead.
#'
#' @param x An MS dial alignment object.
#' @param peaks A character vector specifying the peaks to include in tidy output.
#' If the character vector is named, the names of the vector elements will be
#' used in place of the original peak names.
#' @param metadata A character vector specifying the metadata to include in
#' the tidy output.
#' @param treatments This argument is deprecated as of version 0.3.2. It is
#' synonymous with the new metadata argument which should be used instead.
#' @importFrom dplyr select mutate any_of
#' @importFrom tidyr pivot_longer
#' @return If \code{export} is \code{TRUE}, returns spectrum as \code{data.frame}.
#' Otherwise, no return value.
#' @author Ethan Bass
#' @export

ms_tidy_msdial <- function(x, peaks, metadata, treatments = NULL){
  .Deprecated("ms_reshape_peaktable", package = "tidy_msdial")
  ms_reshape_peaktable(x=x, peaks = peaks, metadata = metadata, treatments = treatments)
}
