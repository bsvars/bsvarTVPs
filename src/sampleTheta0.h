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


#endif
