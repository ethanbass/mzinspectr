
#' Plot mass spectrum of peak given by \code{col}.
#' @param x MS dial alignment object
#' @param col Spectrum to plot
#' @param plot_labels Logical. Whether to plot labels or not.
#' @param lab.int Labels will be plotted above the specified proportion of the
#' largest ion.
#' @param type What kind of plot to produce. Either base R (\code{base}) or
#' plotly (\code{plotly})
#' @param width Width of bars.
#' @param digits How many figures to include on mz labels
#' @importFrom graphics text
#' @return No return value
#' @author Ethan Bass
#' @export

plot_spectrum <- function(x, col, plot_labels=TRUE, lab.int = 0.2,
                          type=c("plotly", "base"), width = 1, digits = 1){
  type <- match.arg(type, c("plotly", "base"))
  spec <- tidy_eispectrum(x$peak_meta["EI.spectrum", col])
  lab.idx <- which(spec$intensity > lab.int*max(spec$intensity))
  if (type == "base"){
    plot(spec, type="h")
    if (plot_labels){
      text(spec$mz[lab.idx],
           spec$intensity[lab.idx],
           round(spec$mz[lab.idx], digits), offset = .25, pos = 3, cex = 0.5)
    }
  } else if (type == "plotly"){
    check_for_pkg("plotly")
    layout <- list(
      # title = "Mass Spectrum for scan 11700",
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Ion Intensity")
    )
    p <- plotly::plot_ly()
    p <- plotly::add_trace(p, type="bar", x=spec$mz, y=spec$intensity, marker=list(line=list(width=width)))
    p <- plotly::layout(p, title=layout$title, xaxis=layout$xaxis, yaxis=layout$yaxis)
    if (plot_labels){
      p <- plotly::add_annotations(p, x = ~spec$mz[lab.idx], y = ~spec$intensity[lab.idx], text=~round(spec$mz[lab.idx], digits), yshift=10, xshift=1, showarrow=FALSE )
    }
    p
  }
}

#' @noRd
check_for_pkg <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(
      "Package", sQuote(pkg), "must be installed to perform this action:
          try", paste0("`install.packages('", pkg, "')`.")),
      call. = FALSE
    )
  }
}


#' @noRd
#' @importFrom graphics boxplot
#' @importFrom stats as.formula
#' @export
boxplot.msdial_alignment <- function(x, rows, cols, vars, ...){
  if (is.null(vars)){
    stop("Must provide independent variable(s) (`var`) to make a boxplot.")
  }
  if (missing(cols)){
    stop("A peak name must be provided to `loc` to make a boxplot.")
  }
  if (missing(rows)){
    rows <- seq_len(nrow(x[["tab"]]))
  }
  for (col in cols){
    boxplot(as.formula(paste("x[['tab']][rows,col]", vars, sep="~")), data=x$sample_meta[rows,],
            las=2, ylab='',xlab='',
            main=paste0("peak ", col, "; rt: ", x$peak_meta["Average.Rt.min.",col],
                        "; mz: ", x$peak_meta["Quant.mass",col]), ...)
  }
}
