
#ifndef _BSVAR_MSSA_TVI_SV_H_
#define _BSVAR_MSSA_TVI_SV_H_

#include <RcppArmadillo.h>

Rcpp::List bsvar_mssa_tvi_sv_cpp (
    const int&                    SS,                   // No. of posterior draws
    const arma::mat&              Y,                    // NxT dependent variables
    const arma::mat&              X,                    // KxT explanatory variables
    const Rcpp::List&             prior,                // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,                   // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100,           // introduce thinning
    const int                     sv_select = 1,        // {1 - non-centred, 2 - centred, 3 - homoskedastic};
    const int                     hyper_select = 1,     // {1 - horseshoe, 2 - boost, 3 - fixed}
    const bool                    finiteM = true,       // {true - stationary MS, false - overfitted};
    const bool                    studentt = false      // {true - normal, false - Student-t};
);

#endif  // _BSVAR_MSSA_TVI_SV_H_