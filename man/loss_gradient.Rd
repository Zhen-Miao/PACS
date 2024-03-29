% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param-estimate_logit_get_p_by_t_June.R
\name{loss_gradient}
\alias{loss_gradient}
\title{Gradient of loss function without Firth prior}
\usage{
loss_gradient(xdumm, p_bg, q_vec, y_vec)
}
\arguments{
\item{xdumm}{A design matrix, X with dummy variables}

\item{p_bg}{A value of open probability}

\item{q_vec}{A vector of capturing rates}

\item{y_vec}{A vector of observed open (1) or close (0) state for each cell,
must match the length of q}
}
\value{
The first order derivative of the loss function without Firth prior,
a vector
}
\description{
Gradient of loss function without Firth prior
}
