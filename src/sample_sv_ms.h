
#ifndef _SAMPLE_SV_MS_H_
#define _SAMPLE_SV_MS_H_

#include <RcppArmadillo.h>


int csample_num1 (
    Rcpp::NumericVector x,
    Rcpp::NumericVector prob = Rcpp::NumericVector::create()
);


arma::vec find_mixture_indicator_cdf (
    const arma::vec& datanorm           // provide all that is conditionally normal
);


arma::uvec inverse_transform_sampling (
    const arma::vec&  mixprob,
    const int         T
);


double do_rgig1(
    double lambda, 
    double chi, 
    double psi
);


Rcpp::List cholesky_tridiagonal(
    const arma::vec&    omega_diag,
    const double&       omega_offdiag
);


arma::vec forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector
);


arma::vec backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp
);


arma::vec precision_sampler_ar1(
    const arma::vec&     precision_diag,
    const double&        precision_offdiag,
    const arma::vec&     location
);


Rcpp::List svar_nc1 (
    arma::rowvec    aux_h_n,                  // 1xT
    double          aux_rho_n,
    double          aux_omega_n,
    double          aux_sigma2v_n,
    double          aux_sigma2_omega_n,       // omega prior hyper-parameter 
    double          aux_s_n,                  // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec   aux_S_n,                  // 1xT
    const arma::rowvec&   u,                  // 1xT
    const Rcpp::List&     prior,
    bool            sample_s_ = true
);
  

Rcpp::List svar_nc1_mss (
    arma::rowvec&         aux_h_n,            // 1xT
    double&               aux_rho_n,
    arma::rowvec&         aux_omega_n,        // 1xM nth equation regime-dependent omegas
    double&               aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&               aux_s_n,            // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&        aux_S_n,            // 1xT
    const arma::mat&      aux_xi,             // MxT
    const arma::rowvec&   u,                  // 1xT
    const Rcpp::List&     prior,
    bool                  sample_s_ = true
);


arma::mat count_regime_transitions (
    const arma::mat& xi
);


arma::rowvec rDirichlet1 (
    const arma::rowvec&   alpha     // Kx1
);


arma::mat filtering (
    const arma::cube& Z,                  // NxTxM state-specific standardised residuals
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  pi_0                // Mx1
);


arma::mat smoothing (
    const arma::mat&  filtered,           // MxT
    const arma::mat&  aux_PR_TR           // MxM
);


arma::mat sample_Markov_process_mss (
    arma::mat         aux_xi,             // MxT
    const arma::mat&  E,                  // NxT
    const arma::cube& aux_B,              // NxNxM
    const arma::mat&  aux_sigma,          // NxM
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  aux_pi_0,           // Mx1
    const bool        finiteM = true
);


Rcpp::List sample_transition_probabilities (
    arma::mat           aux_PR_TR,    // MxM 
    arma::vec           aux_pi_0,     // Mx1
    const arma::mat&    aux_xi,       // MxT
    const Rcpp::List&   prior,         // a list of priors - original dimensions
    const bool          MSnotMIX = true
);


#endif  // _SAMPLE_SV_MS_H_