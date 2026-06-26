
#ifndef _SAMPLE_ABHYPER_TVP_H_
#define _SAMPLE_ABHYPER_TVP_H_

#include <RcppArmadillo.h>


arma::vec mvnrnd_prec_cond_gibbs (
    arma::vec x,          // Nx1 vector to be filled with conditional normal draws when missing == 1
    arma::vec mu,         // Nx1 mean vector
    arma::mat precision,  // NxN precision matrix
    arma::vec to_sample   // Nx1 with 1 for missing observations
);


arma::vec mvnrnd_prec_cond (
    arma::vec x,          // Nx1 vector to be filled with conditional normal draws when missing == 1
    arma::vec mu,         // Nx1 mean vector
    arma::mat precision,  // NxN precision matrix
    arma::vec to_sample   // Nx1 with 1 for missing observations
);


arma::mat orthogonal_complement_matrix_TW (const arma::mat& x);


arma::mat sample_B_heterosk1 (
    arma::mat         aux_B,          // NxN
    const arma::mat&  aux_A,          // NxK
    arma::field<arma::mat> prior_precision, // (N)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
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


double log_posterior_kernel_B (
    arma::mat                     aux_B,          // NxN
    arma::mat                     aux_Theta0_inv,
    arma::mat&                    shocks,         // NxT RF error terms
    arma::mat&                    aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::field<arma::mat>        prior_precision // (N)(N,N)
);


Rcpp::List sample_B_heterosk1_s4 (
    arma::mat                     aux_B,          // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators
    arma::mat                     aux_Theta0,     // NxN
    arma::mat                     aux_Theta0_inv,
    arma::mat&                    shocks,         // NxT RF error terms
    arma::mat&                    aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::field<arma::mat>        prior_precision, // (N)(N,N)
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL             // restrictions on B0 in S4 arrangement
);


Rcpp::List sample_B_mss_s4 (
    arma::cube        aux_B,          // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    arma::cube        aux_Theta0,
    arma::cube        aux_Theta0_inv,
    const arma::mat&  aux_A,          // NxK
    arma::field<arma::mat>  prior_precision, // (N,M)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
);


Rcpp::List sample_B_mssa_s4 (
    arma::cube        aux_B,          // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    arma::cube        aux_Theta0,
    arma::cube        aux_Theta0_inv,
    const arma::cube& aux_A,          // NxKxM
    arma::field<arma::mat>  prior_precision, // (N,M)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
);
  

arma::mat sample_A_heterosk1 (
      arma::mat         aux_A,          // NxK
      const arma::mat&  aux_B,          // NxN
      arma::field<arma::mat>  prior_precision, // (N)(N,N)
      const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
      const arma::mat&  Y,              // NxT dependent variables
      const arma::mat&  X,              // KxT dependent variables
      const Rcpp::List& prior           // a list of priors - original dimensions
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


Rcpp::List sample_hyperparameter_boost_s4 (
    Rcpp::List&             aux_hyper_list,      // (2 * N + 1) x 2
    const arma::mat&        aux_B,
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
);


Rcpp::List sample_hyperparameters_mss_boost (
    Rcpp::List&             aux_hyper_list,      // (2 * N + 1) x 2
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
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


double rig_inv1 (
    double alpha,
    double beta
);


Rcpp::List sample_hyperparameter_horseshoe (
    Rcpp::List&             aux_hyper_list,
    const arma::mat&        aux_B,                // (N, N)
    const arma::mat&        aux_A,                // (N, K)     
    const arma::field<arma::mat>& VB,             // (N)
    const arma::ivec&       aux_SL,               // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
);


Rcpp::List sample_hyperparameter_mss_horseshoe (
    Rcpp::List&             aux_hyper_list,
    const arma::cube&       aux_B,                // (N, N, M)
    const arma::mat&        aux_A,                // (N, K)     
    const arma::field<arma::mat>& VB,             // (N)
    const Rcpp::List&       prior
);


Rcpp::List sample_hyperparameter_mss_s4_horseshoe (
    Rcpp::List&             aux_hyper_list,
    const arma::cube&       aux_B,                // (N, N, M)
    const arma::mat&        aux_A,                // (N, K)     
    const arma::field<arma::mat>& VB,             // (R + 1)
    const arma::imat&       aux_SL,               // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
);


Rcpp::List sample_hyperparameter_mssa_s4_horseshoe (
    Rcpp::List&             aux_hyper_list,
    const arma::cube&       aux_B,                // (N, N, M)
    const arma::cube&       aux_A,                // (N, K, M)     
    const arma::field<arma::mat>& VB,             // (R + 1)
    const arma::imat&       aux_SL,               // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionB_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionA_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionB_mss_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionA_msa_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
);


arma::field<arma::mat> hyper2precisionB_horseshoe (
    Rcpp::List              aux_hyper_list
);


arma::field<arma::mat> hyper2precisionA_horseshoe (
    Rcpp::List              aux_hyper_list
);


arma::field<arma::mat> hyper2precisionB_mss_horseshoe (
    Rcpp::List              aux_hyper_list
);


arma::field<arma::mat> hyper2precisionA_msa_horseshoe (
    Rcpp::List              aux_hyper_list
);


#endif  // _SAMPLE_ABHYPER_TVP_H_