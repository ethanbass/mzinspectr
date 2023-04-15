% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attach_metadata.R
\name{ms_attach_metadata}
\alias{ms_attach_metadata}
\alias{attach_metadata}
\title{Attach experimental metadata}
\usage{
ms_attach_metadata(x, metadata, col)
}
\arguments{
\item{x}{An `MSDIAL_alignment` object.}

\item{metadata}{A `data.frame` containing the sample metadata.}

\item{col}{The name of the column containing the sample names.}
}
\value{
A \code{MSDIAL_alignment} object with attached metadata in the \code{
$sample_meta} slot.
}
\description{
Attaches experimental metadata to `MSDIAL_alignment` object. One of the columns in
the supplied metadata must match exactly the row names of the peak table.
}
\author{
Ethan Bass
}