
#' @title US fiscal system for the period 1948 Q1 -- 2025 Q4
#'
#' @description A system used to identify the US fiscal policy shocks. 
#' Last data update was implemented on 2026-06-17.
#'
#' @usage data(us_fiscal)
#' 
#' @format A matrix and a \code{ts} object with time series of over three hundred observations on 3 variables:
#' \describe{
#'   \item{ttr}{quarterly US total tax revenue expressed in log, real, per person terms}
#'   \item{gs}{quarterly US total government spending expressed in log, real, per person terms}
#'   \item{gdp}{quarterly US gross domestic product expressed in log, real, per person terms}
#'   \item{iv}{quarterly narrative fiscal shocks from Romer & Romer (2010)}
#' }
#' 
#' The series are as described by Mertens & Ravn (2014) in footnote 3 and main body on page S3 of the paper. 
#' Differences with respect to Mertens & Ravn's data :
#' \itemize{
#' \item The sample period is from quarter 1 of 1948 to the last available observation,
#' \item The population variable is not from Francis & Ramey (2009) but from the FRED (with the same definition),
#' \item The original monthly population data is transformed to quarterly by taking monthly averages.
#' }
#' 
#' @references 
#' Francis, N., and Ramey, V.A. (2009) Measures of per capita Hours and Their 
#' Implications for the Technology‐hours Debate. \emph{Journal of Money, Credit and Banking}
#' , 41(6), 1071-1097, DOI: \doi{10.1111/j.1538-4616.2009.00247.x}.
#' 
#' Mertens, K., and Ravn, M.O. (2014) A Reconciliation of SVAR and Narrative 
#' Estimates of Tax Multipliers, \emph{Journal of Monetary Economics}, 68(S), 
#' S1–S19. DOI: \doi{10.1016/j.jmoneco.2013.04.004}.
#' 
#' Romer & Romer (2010) The Macroeconomic Effects of Tax Changes: Estimates 
#' Based on a New Measure of Fiscal Shocks, \emph{American Economic Review}, 
#' 100(3), 763-801. DOI: \doi{10.1257/aer.100.3.763}.
#' 
#' @source 
#' U.S. Bureau of Economic Analysis, National Income and Product Accounts, \url{https://www.bea.gov/}
#' 
#' FRED Economic Database, Federal Reserve Bank of St. Louis, \url{https://fred.stlouisfed.org/}
#' 
#' @examples 
#' plot(us_fiscal)
"us_fiscal"