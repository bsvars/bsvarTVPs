
#ifndef _SAMPLE_ABHYPER_TVP_H_
#define _SAMPLE_ABHYPER_TVP_H_

#include <RcppArmadillo.h>


double rig_inv1 (
    double alpha,
    double beta
);


arma::vec mvnrnd_prec_cond (
    arma::vec x,          // Nx1 vector to be filled with conditional normal draws when missing == 1
    arma::vec mu,         // Nx1 mean vector
    arma::mat precision,  // NxN precision matrix
    arma::vec to_sample   // Nx1 with 1 for missing observations
);


arma::mat orthogonal_complement_matrix_TW (const arma::mat& x);


double log_posterior_kernel_B (
    arma::mat                     aux_B,          // NxN
    arma::mat                     aux_Theta0_inv,
    arma::mat&                    shocks,         // NxT RF error terms
    arma::mat&                    aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::field<arma::mat>        prior_precision // (N)(N,N)
);


arma::mat sample_B_heterosk1_rown (
    int               n,              // row n of B to be sampled
    arma::mat         aux_B,          // NxN
    arma::mat&        shocks,         // NxT conditional STANDARD DEVIATIONS
    arma::mat&        aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::mat         prior_precision_n, // (N,N)
    const Rcpp::List& prior,          // a list of priors - original dimensions
    arma::mat&        VB_n        // restrictions on B0
);


Rcpp::List sample_B_heterosk1_s4 (
    arma::mat                     aux_B,          // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators
    arma::mat                     aux_Theta0,     // NxN
    arma::mat                     aux_Theta0_inv,
    arma::vec                     aux_SLlpr,
    arma::mat&                    shocks,         // NxT RF error terms
    arma::mat&                    aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::field<arma::mat>        prior_precision, // (N)(N,N)
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL             // restrictions on B0 in S4 arrangement
);


arma::mat sample_A_heterosk1_mss (
    arma::mat         aux_A,          // NxK
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    arma::field<arma::mat>  prior_precision, // (N)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VA  // restrictions on A
);


arma::cube sample_A_heterosk1_mssa (
    arma::cube        aux_A,          // NxKxM
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    arma::field<arma::mat>  prior_precision, // (N,M)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VA  // restrictions on A
);


Rcpp::List sample_hyperparameters_mss_s4_boost (
    Rcpp::List&             aux_hyper_list,      // (2 * N + 1) x 2
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
);


Rcpp::List sample_hyperparameters_mssa_s4_boost (
    Rcpp::List&             aux_hyper_list,      // (2 * N + 1) x 2
    const arma::cube&       aux_B,            // NxNxM
    const arma::cube&       aux_A,            // NxKxM
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
);


arma::field<arma::mat> hyper2precisionB_mss_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionA_msa_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionA_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


Rcpp::List construct_LR (
    const arma::mat& zeroB_R
);


double log_posterior_kernel_Theta0 (
    const int n,
    const arma::mat& aux_Theta0,
    const arma::mat& aux_B,
    const arma::mat& shocks,
    const arma::mat& aux_sigma,
    const arma::mat& prior_B0,
    const arma::mat& prior_VB0
);


arma::mat sample_Theta0_Hou24_heterosk1_coln (
    const int         n,
    arma::mat&        aux_Theta0,     // NxN
    const arma::mat&  shocks,         // NxT B(Y-AX)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::mat&  R_E_n,
    const arma::mat&  D_E_n
);


Rcpp::List sample_Theta0_Hou24_heterosk1_s4 (
    arma::mat                     aux_Theta0,     // NxN
    arma::mat                     aux_B,          // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators aux_SL.slice(1).col(m)
    arma::vec                     aux_SLlpr,       // N col os S4 indicators probs aux_SL.slice(1).col(m)
    const arma::mat&              shocks,         // NxT B(Y-AX)
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const Rcpp::List&             restrictions    // output of construct_LR
);


Rcpp::List sample_BTheta0_tvi (
    arma::cube                    aux_B,          // NxNxM
    arma::cube                    aux_Theta0,     // NxNxM
    arma::icube                   aux_SL,         // NxMx2 row-specific S4 indicators
    arma::cube                    aux_SLlpr,       // NxMx2 row-specific S4 indicators for prior 
    const arma::mat&              shocks,         // NxT shocks = Y - aux_A * X;
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              aux_xi,         // MxT
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    arma::field<arma::mat>        prior_precision, // (N,M)(N,N)
    const arma::field<arma::mat>& VB, // restrictions on B0
    const Rcpp::List&             VTheta0
);


#endif  // _SAMPLE_ABHYPER_TVP_H_