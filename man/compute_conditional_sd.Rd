% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsvarTVPTOOLs.R
\name{compute_conditional_sd}
\alias{compute_conditional_sd}
\title{Computes posterior draws of structural shock conditional standard deviations}
\usage{
compute_conditional_sd(posterior)
}
\arguments{
\item{posterior}{posterior estimation outcome - an object of either of the classes: 
PosteriorBSVARSVMSTVI, PosteriorBSVARSVTVI, PosteriorBSVARSVMS, PosteriorBSVARMSTVI, 
PosteriorBSVARTVI, or PosteriorBSVARMS}
}
\value{
An object of class PosteriorSigma, that is, an \code{NxTxS} array with attribute PosteriorSigma 
containing \code{S} draws of the structural shock conditional standard deviations.
}
\description{
Each of the draws from the posterior estimation of a model is transformed into
a draw from the posterior distribution of the structural shock conditional standard deviations.
}
\author{
Tomasz Woźniak \email{wozniak.tom@pm.me}
}
