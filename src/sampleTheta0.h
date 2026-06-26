#ifndef _sampleTheta0_
#define _sampleTheta0_

#include <RcppArmadillo.h>


arma::vec Rodriguez_Yam_2004 (
    arma::vec& mu,
    arma::mat& Sig,
    arma::mat& A,
    arma::vec& l,
    arma::vec& u,
    arma::vec& x
);


double AN (
    const double& mu,
    const double& rho
);


double logf_w1 (
    const double& x,
    const double& P,
    const double& m1,
    const double& m2,
    const double& sig1,
    const double& sig2,
    const double& gam1, 
    const double& gam2,
    const int&    T
);


double logtar_w1 (
    const double& x,
    const double& gam1, 
    const double& gam2,
    const double& invVw10, 
    const double& w10,
    const int&    T
);


double exact_sampling_w1 (
    const double& gam1, 
    const double& gam2, 
    const int&    T, 
    const double& neww1, 
    const double& w1, 
    const double& invVw10, 
    const double& w10
);


Rcpp::List construct_LR (
    const arma::mat& zeroB_R,
    const arma::mat& signB_R
);


arma::mat sample_Theta0_Hou24_heterosk1 (
    arma::mat&        aux_Theta0,     // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const Rcpp::List& restrictions,   // output of construct_LR
    const bool        debug = false
);


arma::mat sample_Theta0_Hou24_heterosk1_coln (
    const int         n,
    arma::mat&        aux_Theta0,     // NxN
    const arma::mat&  shocks,         // NxT B(Y-AX)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::mat&  R_E_n,
    const arma::mat&  D_E_n,
    const Rcpp::List& restrictions,   // containing only R_E(n) and D_E(n)
    const bool        debug = false
);


Rcpp::List sample_Theta0_Hou24_heterosk1_s4 (
    arma::mat                     aux_Theta0,     // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators aux_SL.slice(1).col(m)
    const arma::mat&              shocks,         // NxT B(Y-AX)
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const Rcpp::List&             restrictions    // output of construct_LR
);


Rcpp::List sample_Theta0_mss_s4 (
    arma::cube        aux_Theta0,     // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  shocks,         // NxT
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const Rcpp::List& restrictions
);


Rcpp::List sample_BTheta0_tvi (
    arma::cube&                   aux_B,          // NxNxM
    arma::cube                    aux_Theta0,     // NxNxM
    arma::icube                   aux_SL,         // NxM row-specific S4 indicators
    const arma::mat&              shocks,         // NxT shocks = Y - aux_A * X;
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              aux_xi,         // MxT
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    arma::field<arma::mat>        prior_precision, // (N,M)(N,N)
    const arma::field<arma::mat>& VB, // restrictions on B0
    const Rcpp::List&             VTheta0
);


#endif
