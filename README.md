
# bsvarTVPs

Bayesian Structural Vector Autoregressions with Time-Varying
Identification

Efficient algorithms for Bayesian estimation of Structural Vector
Autoregressions (VARs) with Stochastic Volatility heteroskedasticity,
Markov-switching and Time-Varying Identification of the Structural
Matrix, and a three-level global-local hierarchical prior shrinkage for
the structural and autoregressive matrices. The models were developed
for a paper by [Camehl & Woźniak
(2023)](https://doi.org/10.48550/arXiv.2311.05883). The ‘bsvarTVPs’
package is aligned regarding objects, workflows, and code structure with
the R packages ‘bsvars’ by [Woźniak
(2024)](https://doi.org/10.32614/CRAN.package.bsvars) and ‘bsvarSIGNs’
by [Wang & Woźniak
(2024)](https://doi.org/10.32614/CRAN.package.bsvarSIGNs), and they
constitute an integrated toolset.

# Installation

To install the **bsvarTVPs** package just type in **R**:

    devtools::install_github("bsvars/bsvarTVPs")

# Checks

[![R-CMD-check](https://github.com/bsvars/bsvarTVPs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bsvarTVPs/actions/workflows/R-CMD-check.yaml)
