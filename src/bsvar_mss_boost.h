
#ifndef _BSVAR_MSS_BOOST_H_
#define _BSVAR_MSS_BOOST_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

Rcpp::List bsvar_mss_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100,  // introduce thinning
    const bool                    hyper_boost = true
);

#endif  // _BSVAR_MSS_BOOST_H_