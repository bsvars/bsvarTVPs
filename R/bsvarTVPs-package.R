#  #####################################################################################
#  R package bsvarTVPs by Tomasz Woźniak Copyright (C) 2022
#
#  This file is part of the R package bsvarTVPs: Bayesian Estimation of 
#  Heteroskedastic Structural Vector Autoregressions with Markov-Switching and 
#  Time-Varying Identification of the Structural Matrix
#
#  The R package bsvarTVPs is free software: you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 or
#  any later version of the License.
#
#  The R package bsvarTVPs is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the R package bsvarTVPs. If that is not the case, please
#  refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################
#
#' @title Bayesian Structural Vector Autoregressions with Time-Varying Identification
#'
#' @description Efficient algorithms for Bayesian estimation of Structural 
#' Vector Autoregressions (VARs) with Stochastic Volatility heteroskedasticity, 
#' Markov-switching and Time-Varying Identification of the Structural Matrix, 
#' and a three-level global-local hierarchical prior shrinkage for the 
#' structural and autoregressive matrices. The models were developed for a paper 
#' by Camehl & Woźniak (2025) <doi:10.48550/arXiv.2502.19659>. The 'bsvarTVPs' 
#' package is aligned regarding objects, workflows, and code structure with the 
#' R packages 'bsvars' by Woźniak (2024) <doi:10.32614/CRAN.package.bsvars> and 
#' 'bsvarSIGNs' by Wang & Woźniak (2024) <doi:10.32614/CRAN.package.bsvarSIGNs>, 
#' and they constitute an integrated toolset.
#' 
#' @details 
#' All the SVAR models in this package are specified by two equations, including 
#' the reduced form equation:
#' \deqn{Y = AX + E}
#' where \eqn{Y} is an \code{NxT} matrix of dependent variables, 
#' \eqn{X} is a \code{KxT} matrix of explanatory variables, 
#' \eqn{E} is an \code{NxT} matrix of reduced form error terms, 
#' and \eqn{A} is an \code{NxK} matrix of autoregressive slope coefficients and 
#' parameters on deterministic terms in \eqn{X}.
#' 
#' The structural equation is given by:
#' \deqn{B(s_t, kappa(s_t))E = U}
#' where \eqn{U} is an \code{NxT} matrix of structural form error terms, and
#' \eqn{B} is an \code{NxN} matrix of contemporaneous relationships.
#' 
#' Finally, all of the models share the following assumptions regarding the structural
#' shocks \code{U}, namely, joint conditional normality given the past observations collected
#' in matrix \code{X}, and temporal and contemporaneous independence. The latter implies 
#' zero correlations and autocorrelations. 
#' 
#' @name bsvarTVPs-package
#' @aliases bsvarTVPs-package bsvarTVPs
#' @useDynLib bsvarTVPs, .registration = TRUE
#' 
#' @importFrom bsvars specify_bsvar_sv specify_data_matrices estimate forecast compute_impulse_responses compute_fitted_values compute_historical_decompositions compute_structural_shocks compute_variance_decompositions compute_regime_probabilities compute_conditional_sd
#' @importFrom HDInterval hdi
#' @importFrom GIGrvg rgig
#' @importFrom R6 R6Class
#' @importFrom RcppTN rtn
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppProgress
#' 
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me} & Annika Camehl \email{camehl@ese.eur.nl}
#' 
#' @references
#' 
#' Camehl, A. & Woźniak, T. (2025) Time-Varying Identification of Structural Vector Autoregressions, <doi:10.48550/arXiv.2502.19659>.
#' 
#' Wang & Woźniak (2024) bsvarSIGNs: Bayesian SVARs with Sign, Zero, and Narrative Restrictions. 
#' R package version 1.0.1, <doi:10.32614/CRAN.package.bsvarSIGNs>.
#' 
#' Woźniak (2024) bsvars: Bayesian Estimation of Structural Vector Autoregressive Models. 
#' R package version 3.2, <doi:10.32614/CRAN.package.bsvars>.
#' 
#' @keywords package models ts
#' 
#' @examples
#' # upload data
NULL
