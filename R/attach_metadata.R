#' Attach experimental metadata
#'
#' Attaches experimental metadata to `MSDIAL_alignment` object. One of the columns in
#' the supplied metadata must match exactly the row names of the peak table.
#'
#' @aliases attach_metadata
#' @param x An `MSDIAL_alignment` object.
#' @param metadata A `data.frame` containing the sample metadata.
#' @param col The name of the column containing the sample names.
#' @return A \code{MSDIAL_alignment} object with attached metadata in the \code{
#' $sample_meta} slot.
#' @author Ethan Bass
#' @examples
#' @export attach_metadata

attach_metadata <- function(x, metadata, col){
  if (any(grepl("tbl", class(metadata)))){
    metadata <- as.data.frame(metadata)
  }
  if (!inherits(metadata, "data.frame")){
    stop("Please provide metadata as a `data.frame`")
  }
  if (!(col %in% colnames(metadata)))
    stop(paste0("Column, ", col, ", is not found."))
  if (sum((duplicated(metadata[,col]))) > 0)
    stop(paste("Sample names must be unique. Please check column", sQuote(col),
               "for duplicates."))
  if (!inherits(x,"msdial_alignment"))
    stop(paste("Provided peak table object must be of class 'x'."))
  meta <- data.frame(rownames(x$tab))
  names(meta) <- col
  metadata[, col] <- as.character(metadata[, col])
  missing_meta <- !(meta[, col] %in% metadata[, col])
  if (sum(missing_meta) > 0)
    warning("The supplied metadata does not include all samples.")
  meta <- keep_order(meta, merge, y = metadata, by = col,
                     all.x = TRUE, all.y = FALSE, sort = FALSE)
  if (any(!is.na(x$sample_meta))){
    x$sample_meta <- cbind(x$sample_meta, meta)
  } else{
    x$sample_meta <- meta
  }
  return(x)
}

keep_order <- function(data, fn, ...) {
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data)
  out <- fn(data, ...)
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function")
  out <- out[order(out[,col]),]
  out[,col] <- NULL
  out
}
