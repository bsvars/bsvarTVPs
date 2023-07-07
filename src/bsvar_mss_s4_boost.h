
#ifndef _BSVAR_MSS_S4_BOOST_H_
#define _BSVAR_MSS_S4_BOOST_H_

#include <RcppArmadillo.h>


Rcpp::List bsvar_mss_s4_boost_cpp (
    const int&                    SS,         // No. of posterior draws
    const arma::mat&              Y,          // NxT dependent variables
    const arma::mat&              X,          // KxT explanatory variables
    const Rcpp::List&             prior,      // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,        // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100  // introduce thinning
);


#endif  // _BSVAR_MSS_S4_BOOST_H_