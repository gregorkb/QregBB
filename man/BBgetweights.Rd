% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QregBB.R
\name{BBgetweights}
\alias{BBgetweights}
\title{Retrieve weights from from Monte-Carlo block samples for the MBB and ETBB}
\usage{
BBgetweights(n, l, B, seed = NA, c = 0.43)
}
\arguments{
\item{n}{length of observed time series.}

\item{l}{block length.}

\item{B}{the number of Monte-Carlo draws.}

\item{seed}{a random number generating seed (optional).}

\item{c}{a constant controlling the tapering of weight assignments over the blocks.}
}
\value{
a list containing \code{n} by \code{B} matrices containing the the mass placed on the \code{n} observations by the MBB and ETBB empirical distributions in each of the \code{B} Monte-Carlo draws, the sampled indices for the block starting points, and the scalar \eqn{m_l} corresponding to the weight function and the block length.
}
\description{
Retrieve weights from from Monte-Carlo block samples for the MBB and ETBB
}
