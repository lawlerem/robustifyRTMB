% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustifyRTMB.R
\name{robustify}
\alias{robustify}
\alias{robust_weight}
\title{Robustify likelihood components}
\usage{
robustify(x, robustness, robust_function = "id")

robust_weight(x, robustness, robust_function = "id")
}
\arguments{
\item{x}{A vector of log-likelihood values}

\item{robustness}{A positive number. Larger numbers downweight components more.}

\item{robust_function}{Either "id" for the identity (non-robust), "ll" for log-logistic, or
"ssh" for smooth semi-huber. "ssh" is currently an experimental function
and may quietly not work as intended.}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{The robustified log-likelihood values.
}\if{html}{\out{</div>}}

\if{html}{\out{<div class="sourceCode">}}\preformatted{A vector of robust weights.
}\if{html}{\out{</div>}}
}
\description{
Compute the robustified likelihood value
}
\section{Functions}{
\itemize{
\item \code{robust_weight()}: Compute the robust weight of log-likelihood values

}}
