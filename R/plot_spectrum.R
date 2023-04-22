#' Plot mass spectrum of peak given by \code{col}.
#' @param x MS dial alignment object
#' @param col Spectrum to plot
#' @param plot_labels Logical. Whether to plot labels or not.
#' @param lab_int Labels will be plotted above the specified proportion of the
#' largest ion.
#' @param type What kind of plot to produce. Either base R (\code{base}) or
#' plotly (\code{plotly})
#' @param scale Logical. Whether to scale mass spectrum. Defaults to FALSE.
#' @param bar_width Width of bars.
#' @param digits How many figures to include on mz labels
#' @param ... Additional arguments.
#' @importFrom graphics text
#' @return If \code{export} is \code{TRUE}, returns spectrum as \code{data.frame}.
#' Otherwise, no return value.
#' @author Ethan Bass
#' @export

ms_plot_spectrum <- function(x, col, plot_labels = TRUE, lab_int = 0.2,
                          type = c("plotly", "base"), scale = FALSE,
                          bar_width = 1, digits = 1, ...){
  if (inherits(x,"msdial_alignment")){
    xx <- x$peak_meta
  } else if (inherits(x, "data.frame") && "EI.spectrum" %in% colnames(x)){
    xx <- x
  }
  type <- match.arg(type, c("plotly", "base"))
  spec <- tidy_eispectrum(xx[col, "EI.spectrum"])
  if (scale){
    spec$intensity <- scales::rescale(spec$intensity)
  }
  lab.idx <- which(spec$intensity > lab_int*max(spec$intensity))
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
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Ion Intensity")
    )
    p <- plotly::plot_ly()
    p <- plotly::add_trace(p, type="bar", x=spec$mz, y=spec$intensity, marker=list(line=list(width=bar_width)))
    p <- plotly::layout(p, title=layout$title, xaxis=layout$xaxis, yaxis=layout$yaxis)
    if (plot_labels){
      p <- plotly::add_annotations(p, x = ~spec$mz[lab.idx], y = ~spec$intensity[lab.idx], text=~round(spec$mz[lab.idx], digits), yshift=10, xshift=1, showarrow=FALSE )
    }
    p
  }
}


#' Mirror plot function
#' @export
ms_mirror_plot <- function(x, ...){
  UseMethod("ms_mirror_plot", x)
}

#' Plot two spectra as a mirror plot.
#' @importFrom graphics lines
#' @param x Mass spectrum as data.frame with m/z values in column one and
#' ionization intensity in column two.
#' @param y Mass spectrum as data.frame with m/z values in column one and
#' ionization intensity in column two.
#' @param type What kind of plot to produce. Either base R (\code{base}) or
#' plotly (\code{plotly}).
#' @param scale Logical. Whether to scale mass spectrum. Defaults to TRUE.
#' @param plot_labels Logical. Whether to label m/z values on plot.
#' @param lab_int Labels will be plotted above the specified proportion of the
#' largest ion.
#' @param digits How many figures to include on m/z labels.
#' @param bar_width Width of bars.
#' @param match_score Logical. Whether to plot match score or not.
#' @param ... Additional arguments
#' @rdname ms_mirror_plot
#' @method ms_mirror_plot data.frame
#' @export

ms_mirror_plot.data.frame <- function(x, y, plot_labels = TRUE, type = c("plotly", "base"),
                        scale = TRUE, lab_int = 0.2, digits = 1, bar_width = 1,
                        match_score = TRUE, ...){
  type <- match.arg(type, c("plotly", "base"))
  colnames(x) <- c("mz","intensity")
  colnames(y) <- c("mz","intensity")
  if (scale){
    x$intensity <- scales::rescale(x$intensity)
    y$intensity <- -scales::rescale(y$intensity)
  }
  y <- y[which(y$mz > min(y$mz)),]
  lab1.idx <- which(x$intensity > lab_int*max(x$intensity))
  lab2.idx <- which(abs(y$intensity) > lab_int*max(abs(y$intensity)))
  if (type == "base"){
    max.int <- max(c(x$intensity,y$intensity))
    plot(x, type = "h", ylim = c(-max.int,max.int))
    lines(y, type = "h")
    if (plot_labels){
      text(x$mz[lab1.idx],
           x$intensity[lab1.idx],
           round(x$mz[lab1.idx], digits), offset = .25, pos = 3, cex = 0.5)

      text(y$mz[lab2.idx],
           y$intensity[lab2.idx],
           round(y$mz[lab2.idx], digits), offset = -.25, pos = 3, cex = 0.5)
    }
    if (match_score){
      y[,2] <- -y[,2]
      match <- try(OrgMassSpecR::SpectrumSimilarity(spec.top = x, spec.bottom = y,
                                                    print.graphic = FALSE))
      label <- paste0("Similarity score: ", format(round(match, 2), nsmall = 2), "%")
      legend("top", legend = label, cex = 0.7, bty = "n")
    }
  } else if (type == "plotly"){
    check_for_pkg("plotly")
    layout <- list(
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Ion Intensity")
    )
    p <- plotly::plot_ly()
    p <- plotly::add_trace(p, type="bar", x = x$mz, y = x$intensity,
                           marker = list(line = list(width = bar_width)),
                           name = "target")
    p <- plotly::add_trace(p, type="bar", x = y$mz, y = y$intensity,
                           marker = list(line = list(width = bar_width)),
                           name = "reference")
    p <- plotly::layout(p, title = layout$title,
                        xaxis = layout$xaxis, yaxis = layout$yaxis)
    if (plot_labels){
      p <- plotly::add_annotations(p, x = ~x$mz[lab1.idx],
                                   y = ~x$intensity[lab1.idx],
                                   text = ~round(x$mz[lab1.idx], digits),
                                   yshift = 10, xshift = 1, showarrow = FALSE,
                                   font = list(size = 8))
      p <- plotly::add_annotations(p, x = ~y$mz[lab2.idx],
                                   y = ~y$intensity[lab2.idx],
                                   text = ~round(y$mz[lab2.idx], digits),
                                   yshift = -10, xshift = 1, showarrow = FALSE,
                                   font = list(size = 8))
    }
    p
  }
}

#' Plot two spectra as a mirror plot.
#' @param x A \code{msdial_alignment} object.
#' @param cols One or more columns in the peak table \code{tab} to plot.
#' @param ref A row in the matches slot corresponding to the provided column.
#' @param type What kind of plot to produce. Either base R (\code{base}) or
#' plotly (\code{plotly}).
#' @param scale Logical. Whether to scale mass spectrum. Defaults to TRUE.
#' @param plot_labels Logical. Whether to label m/z values on plot.
#' @param lab_int Labels will be plotted above the specified proportion of the
#' largest ion.
#' @param digits How many figures to include on m/z labels.
#' @param bar_width Width of bars.
#' @param match_score Logical. Whether to plot match score or not.
#' @rdname ms_mirror_plot
#' @method ms_mirror_plot msdial_alignment
#' @export

ms_mirror_plot.msdial_alignment <- function(x, cols, ref, type=c("plotly", "base"),
                                         scale = TRUE, plot_labels = TRUE,
                                         lab_int = 0.2, digits = 1,
                                         bar_width = 1, match_score = TRUE, ...){
  type <- match.arg(type, c("plotly", "base"))
  if (length(cols) > 2){
    stop("Please provide a maximum of two columns.")
  }
  spec1 <- tidy_eispectrum(x$peak_meta[cols[1], "EI.spectrum"])
  if (length(cols) == 2){
    spec2 <- tidy_eispectrum(x$peak_meta[cols[2], "EI.spectrum"])
  } else if (!missing(ref) && !all(is.na(x$matches[[cols[1]]]))){
    spec2 <- tidy_eispectrum(x$matches[[cols[1]]][ref, "Spectra"])
  }
  ms_mirror_plot(x = spec1, y = spec2, plot_labels = plot_labels,
              lab_int = lab_int, type = type, scale = scale, bar_width = bar_width,
              digits = digits, match_score = match_score)
}
