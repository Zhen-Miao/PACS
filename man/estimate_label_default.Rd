% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_type_annotation.R
\name{estimate_label_default}
\alias{estimate_label_default}
\title{Estimate cell type label for new data}
\usage{
estimate_label_default(r_by_t, in_r_by_c, alpha = 1)
}
\arguments{
\item{r_by_t}{Region-by cell type matrix generated from a annotated dataset}

\item{in_r_by_c}{Input (unannotated) region by cell matrix}

\item{alpha}{Weight for negative peaks, default = 1}
}
\value{
A matrix of cell types by cell matrix, with elements representing
probability of being in that cell type
}
\description{
Estimate cell type label for new data
}
