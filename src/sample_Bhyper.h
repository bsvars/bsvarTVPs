
#ifndef _SAMPLE_BHYPER_H_
#define _SAMPLE_BHYPER_H_

#include <RcppArmadillo.h>


void sample_B_heterosk1_s4 (
    arma::mat&        aux_B,          // NxN
    arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const arma::mat&  aux_A,          // NxK
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL       // restrictions on B0 in S4 arrangement
);



void sample_hyperparameters_s4 (
    arma::vec&              aux_hyper,
    const arma::mat&        aux_B,
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
);


#endif  // _SAMPLE_BHYPER_H_