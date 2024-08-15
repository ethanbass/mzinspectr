#' Plot spectrum from peak table
#'
#' Plots the trace and/or spectrum for a given peak in peak table.
#'
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Can also be used
#' to create simple box-plots to examine the distribution of a peak with respect
#' to variables defined in sample metadata.
#'
#' @importFrom graphics title boxplot
#' @importFrom stats as.formula
#' @param x A \code{ms.alignment} object.
#' @param col A vector specifying the peak(s) that you wish to plot.
#' @param plot_spectrum Logical. If TRUE, plots the mass spectrum of the chosen
#' peak. Defaults to TRUE.
#' @param box_plot Logical. If TRUE, plots box plot using factors
#' defined by \code{vars}.
#' @param vars Independent variables for boxplot. Righthand side of formula.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param engine Which plotting engine to use: either \code{base} or \code{plotly}.
#' @param ... Additional arguments to \code{\link[graphics]{boxplot}}.
#' @return No return value.
#' @section Side effects:
#'
#' If \code{plot_spectrum} is TRUE, plots the spectrum for the specified chromatogram
#' at the specified retention time. The spectrum is a single row from the chromatographic
#' matrix.
#'
#' If \code{box_plot} is TRUE, produces a \code{\link[graphics]{boxplot}} from the
#' specified peak with groups provided by \code{vars}.
#'
#' @author Ethan Bass
#' @rdname plot.peak_table
#' @concept Visualization
#' @export

plot.ms_alignment <- function(x, col, plot_spectrum = TRUE,
                            box_plot = FALSE,  vars = NULL,
                            spectrum_labels = TRUE,
                            engine = c("base", "plotly"), ...){
  engine <- match.arg(engine, c("base", "plotly"))
  for (col in col){
    if (plot_spectrum){
      ms_plot_spectrum(x = x, col = col, plot_labels = spectrum_labels,
                       type = engine)
    }
    if (box_plot){
      if (!is.data.frame(x$sample_meta)){
        stop("Attach metadata to `peak_table` to make a boxplot.")
      }
      if (is.null(vars)){
        stop("Must provide independent variable(s) (`var`) to make a boxplot.")
      }
      boxplot(as.formula(paste("x[['tab']][,col]", vars, sep="~")),
              data = x$sample_meta,
              main = paste(col, '\n', 'RT = ',
                           round(as.numeric(x$peak_meta[col, "Average.Rt.min."]), 2)),
              ylab = "abs", xlab = "", ...)
    }
  }
}

#' Make boxplot from MS peak table.
#'
#' The function can take multiple response variables on the left hand side of the
#' formula (separated by \code{+}). In this case, a separate boxplot will be
#' produced for each response variable.
#'
#' @param x A peak_table object
#' @param formula A formula object
#' @param ... Additional arguments to \code{\link[graphics]{boxplot}}
#' @importFrom stats reformulate terms
#' @importFrom graphics boxplot
#' @concept Visualization
#' @export

boxplot.ms_alignment <- function(x, formula, ...){
  if (missing(formula)){
    stop("Please provide a `formula` to make a boxplot.")
  }
  lhs <- get_lhs_vars(formula)
  for (li in lhs){
    response <- paste0("x[['tab']][,'", li, "']")
    form <- reformulate(labels(terms(formula)), response = response)

    idx <- ifelse(!is.na(x$peak_meta[col,"Average.RI"]),
                  paste0("RI: ", x$peak_meta[col,"Average.RI"]),
                  paste0("RT: ", x$peak_meta[col,"Average.Rt.min."])
    )
    title <- paste(li, idx, sep = "\n")
    boxplot(form,
            data = x$sample_meta,
            main = title,
            ylab = "", xlab = "", ...)
  }
}

#' Extract variables from the left-hand-side of a formula.
#' @param formula A \code{\link{formula}} object.
#' @importFrom Formula Formula
#' @noRd
#' @note Adapted from https://github.com/adibender/pammtools/blob/master/R/formula-utils.R
#' (c) Copyright © 2017 Andreas Bender and Fabian Scheipl under MIT license:
#' Permission is hereby granted, free of charge, to any person obtaining a copy of
#' this software and associated documentation files (the “Software”), to deal in
#' the Software without restriction, including without limitation the rights to use,
#' copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
#' Software, and to permit persons to whom the Software is furnished to do so,
#' subject to the following conditions:
#'
#' The above copyright notice and this permission notice shall be included in all
#' copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#' IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#' FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#' AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#' WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#' CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

get_lhs_vars <- function(formula) {
  if (is.character(formula) ) formula <- as.formula(formula)
  form <- formula(Formula::Formula(formula), lhs = TRUE, rhs = FALSE)
  all.vars(form)
}
