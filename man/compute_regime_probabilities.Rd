% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsvarTVPTOOLs.R
\name{compute_regime_probabilities}
\alias{compute_regime_probabilities}
\title{Computes posterior draws of regime probabilities}
\usage{
compute_regime_probabilities(
  posterior,
  Y,
  X,
  type = c("realized", "filtered", "forecasted", "smoothed")
)
}
\arguments{
\item{posterior}{posterior estimation outcome of regime-dependent heteroskedastic models 
- an object of either of the classes: PosteriorBSVARSVMSTVI, or PosteriorBSVARSVMS.}

\item{Y}{a \code{N x T} matrix of dependent variables}

\item{X}{a \code{K x T} matrix of lagged variables}

\item{type}{one of the values \code{"realized"}, \code{"filtered"}, \code{"forecasted"}, or \code{"smoothed"}
denoting the type of probabilities to be computed.}
}
\value{
An object of class PosteriorRegimePr, that is, an \code{MxTxS} array with attribute PosteriorRegimePr 
containing \code{S} draws of the regime probabilities.
}
\description{
Each of the draws from the posterior estimation of a model is transformed into
a draw from the posterior distribution of the regime probabilities. These represent either
the realisations of the regime indicators, when \code{type = "realized"}, filtered probabilities,
when \code{type = "filtered"}, forecasted regime probabilities, when \code{type = "forecasted"},
or the smoothed probabilities, when \code{type = "smoothed"}, .
}
\references{
Song, Y., and Woźniak, T., (2021) Markov Switching. \emph{Oxford Research Encyclopedia of Economics and Finance}, Oxford University Press, \doi{https://doi.org/10.1093/acrefore/9780190625979.013.174}.
}
\author{
Tomasz Woźniak \email{wozniak.tom@pm.me}
}
