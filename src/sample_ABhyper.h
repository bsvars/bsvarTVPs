
#ifndef _SAMPLE_ABHYPER_TVP_H_
#define _SAMPLE_ABHYPER_TVP_H_

#include <RcppArmadillo.h>


void sample_B_mss (
    arma::cube&       aux_B,          // NxNxM
    const arma::mat&  aux_A,          // NxK
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
);


void sample_B_heterosk1_s4 (
    arma::mat&                    aux_B,          // NxN
    arma::ivec&                   aux_SL,         // Nx1 row-specific S4 indicators
    const arma::mat&              aux_A,          // NxK
    const arma::vec&              aux_hyper,      // NxM
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              Y,              // NxT dependent variables
    const arma::mat&              X,              // KxT dependent variables
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL       // restrictions on B0 in S4 arrangement
);


void sample_B_mss_s4 (
    arma::cube&       aux_B,          // NxNxM
    arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const arma::mat&  aux_A,          // NxK
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
);


void sample_A_heterosk1_mss (
    arma::mat&        aux_A,          // NxK
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior           // a list of priors - original dimensions
);


void sample_hyperparameters_s4 (
    arma::vec&              aux_hyper,
    const arma::mat&        aux_B,
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
);


void sample_hyperparameters_mss (
    arma::vec&              aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const Rcpp::List&       prior
);


void sample_hyperparameters_mss_s4 (
    arma::vec&              aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior
);

#endif  // _SAMPLE_ABHYPER_TVP_H_