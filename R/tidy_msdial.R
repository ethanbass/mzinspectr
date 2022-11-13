#' Convert peak to tidy format for plotting
#' @param x An MS dial alignment object
#' @param peaks Peaks to include in tidy output
#' @param treatments Treatments to include in tidy output
#' @importFrom dplyr select mutate any_of
#' @importFrom tidyr pivot_longer
#' @return If \code{export} is \code{TRUE}, returns spectrum as \code{data.frame}.
#' Otherwise, no return value.
#' @author Ethan Bass
#' @export

tidy_msdial <- function(x, peaks, treatments){
  df <- as.data.frame(x$tab)
  if (is.numeric(peaks)){
    peaks <- colnames(df)[peaks]
  }
  df <- select(df, any_of(peaks))
  df <- cbind(df, select(x$sample_meta, any_of(treatments)))
  df <- mutate(df, sample = row.names(x))
  df <- pivot_longer(df, cols=peaks, names_to="peak")
  df
}
