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
#' @title Bayesian Estimation of Heteroskedastic Structural Vector Autoregressions 
#' with Markov-Switching and Time-Varying Identification of the Structural Matrix
#'
#' @description Efficient algorithms for Bayesian estimation of Structural 
#' Vector Autoregressions with Stochastic Volatility heteroskedasticity, 
#' Markov-switching and Time-Varying Identification of the Structural Matrix, 
#' and a three-level global-local hierarchical prior shrinkage for the structural 
#' and autoregressive matrices.
#' The models were developed for a paper 
#' Camehl, Annika & Woźniak, Tomasz (2022) What do Data Say About 
#' Time-Variation in Monetary Policy Shock Identification?
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
#' @docType package
#' @useDynLib bsvarTVPs, .registration = TRUE
#' @importFrom bsvars specify_bsvar_sv
#' @importFrom HDInterval hdi
#' @importFrom GIGrvg rgig
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @import RcppProgress
#' @importFrom RcppTN rtn
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' @references
#' Camehl, A. & Woźniak, T. (2022) What do Data Say About Time-Variation in Monetary Policy Shock Identification?
#' @keywords package models ts
#' 
#' @examples
#' # upload data
NULL
