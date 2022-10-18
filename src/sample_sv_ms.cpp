
#include <RcppArmadillo.h>
#include "RcppTN.h"
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
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
) {
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_cols;
  const int     M = aux_xi.n_rows;
  const rowvec  U = log(square(u + ccc));
  const vec    Tm= sum(aux_xi, 1);
  
  rowvec        omega_aux(M);
  
  const double  prior_sv_a_ = prior["sv_a_"];
  const double  prior_sv_s_ = prior["sv_s_"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho_n;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  rowvec        omega_T(T);
  for (int t=0; t<T; t++) {
    omega_T(t)          = aux_omega_n(aux_xi.col(t).index_max());
  }
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = bsvars::find_mixture_indicator_cdf(trans(U - aux_h_n % omega_T));
  aux_S_n               = trans(bsvars::inverse_transform_sampling(mixprob, T));
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = alpha_s(aux_S_n(t));
    sigma_S_inv.col(t)  = 1/sigma_s(aux_S_n(t));
  }
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_n               = (prior_sv_s_ + 2 * aux_sigma2_omega_n)/R::rchisq(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega_n    = bsvars::do_rgig1( prior_sv_a_ - 0.5 * M, accu(square(aux_omega_n)), 2/aux_s_n );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  double    aux_rho_var = as_scalar(pow(hm1*hm1.t(), -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  double    upper_bound = sqrt(1-aux_sigma2_omega_n);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5), -upper_bound, upper_bound);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_omega
  for (int m=0; m<M; m++) {
    rowvec  aux_h_n_m(Tm(m));
    rowvec  U_m(Tm(m));
    rowvec  alpha_S_m(Tm(m));
    rowvec  sigma_S_inv_m(Tm(m));
    
    int ii = 0;
    for (int t=0; t<T; t++) {
      if (aux_xi(m,t)==1) {
        aux_h_n_m(ii)   = aux_h_n(t);
        U_m(ii)         = U(t);
        alpha_S_m(ii)   = alpha_S(t);
        sigma_S_inv_m(ii) = sigma_S_inv(t);
        ii++;
      }
    }
    double    V_omega_inv = 1/( as_scalar(aux_h_n_m * diagmat(sigma_S_inv_m) * aux_h_n_m.t()) + pow(aux_sigma2_omega_n, -1) );
    double    omega_bar   = as_scalar(aux_h_n_m * diagmat(sigma_S_inv_m) * (U_m - alpha_S_m).t());
    omega_aux(m)          = R::rnorm(V_omega_inv*omega_bar, sqrt(V_omega_inv) );
  } // END m loop
  
  for (int t=0; t<T; t++) {
    omega_T(t)          = omega_aux(aux_xi.col(t).index_max());
  }
  
  // sample aux_h
  mat       V_h         = diagmat(square(omega_T)) * diagmat(sigma_S_inv) + HH_rho;
  vec       h_bar       = diagmat(omega_T) * diagmat(sigma_S_inv) * (U - alpha_S).t();
  rowvec    h_aux       = trans(bsvars::precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar));
  
  // ASIS
  rowvec    aux_h_tilde = h_aux % omega_T;
  rowvec    aux_h_tilde_ar = trans(H_rho * aux_h_tilde.t());
  
  for (int m=0; m<M; m++) {
    double  hHHh = 0;
    for (int t=0; t<T; t++) {
      if (aux_xi(m,t)==1) {
        hHHh           += pow(aux_h_tilde(t), 2);
      }
    }
    double  sigma2_aux  = bsvars::do_rgig1( -0.5*(Tm(m)-1), hHHh, 1/aux_sigma2_omega_n );
    int     ss = 1;
    if (R::runif(0,1)<0.5) ss *= -1;
    aux_omega_n(m)      = ss * sqrt(sigma2_aux);
  } // END m loop
  
  for (int t=0; t<T; t++) {
    omega_T(t)          = aux_omega_n(aux_xi.col(t).index_max());
  }
  aux_h_n               = aux_h_tilde % pow(omega_T, -1);
  
  // resample aux_rho
  hm1                   = aux_h_n.cols(0,T-2);
  aux_rho_var           = as_scalar(pow(hm1*hm1.t(), -1));
  aux_rho_mean          = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  upper_bound           = sqrt(1-aux_sigma2_omega_n);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5), -upper_bound, upper_bound);
  
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END sv_nc1




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat filtering (
    const arma::cube& Z,                  // NxTxM state-specific standardised residuals
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  pi_0                // Mx1
) {
  
  // filtered probabilities for a model with MS structural matrix of SVAR-MSS-SV model
  
  const int   T = Z.n_cols;
  const int   N = Z.n_rows;
  const int   M = aux_PR_TR.n_rows;
  
  mat         eta_t(M, T);
  mat         xi_t_t(M, T);
  
  // This loop evaluates mvnormal pdf at Z being multivariate standard normal distribution
  for (int m=0; m<M; m++) {
    rowvec log_d    = -0.5 * sum(square(Z.slice(m)), 0);
    log_d          += -0.5 * N * log(2*M_PI);
    NumericVector   exp_log_d   = wrap(exp(log_d));
    exp_log_d[exp_log_d==0]     = 1e-300;
    eta_t.row(m)    = as<rowvec>(exp_log_d);
  } // END m loop
  
  vec xi_tm1_tm1    = pi_0;
  
  for (int t=0; t<T; t++) {
    vec     num     = eta_t.col(t) % (aux_PR_TR.t() * xi_tm1_tm1);
    double  den     = sum(num);
    xi_t_t.col(t)   = num/den;
    xi_tm1_tm1      = xi_t_t.col(t);
  } // END t loop
  
  return xi_t_t;
} // END filtering



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat smoothing (
    const arma::mat&  filtered,           // MxT
    const arma::mat&  aux_PR_TR           // MxM
) {
  // NOT the same as for msh (but could be the same if you get rid of arg U in the other one)
  const int   T = filtered.n_cols;
  const int   M = aux_PR_TR.n_rows;
  
  mat   smoothed(M, T);
  smoothed.col(T-1)   = filtered.col(T-1);
  
  for (int t=T-2; t>=0; --t) {
    smoothed.col(t)   = (aux_PR_TR * (smoothed.col(t+1)/(aux_PR_TR.t() * filtered.col(t)) )) % filtered.col(t);
  } // END t loop
  
  return smoothed;
} // smoothing










// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
void sample_Markov_process_mss (
    arma::mat&       aux_xi,             // MxT
    const arma::mat&  E,                  // NxT
    const arma::cube& aux_B,              // NxNxM
    const arma::mat&  aux_sigma,          // NxM
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  aux_pi_0,           // Mx1
    const bool        finiteM = true
) {
  // the function changes the value of aux_xi by reference (filling it with a new draw)
  
  int minimum_regime_occurrences = 10;
  int max_iterations = 50;
  if ( finiteM ) {
    minimum_regime_occurrences = 10;
    max_iterations = 50;
  }
  
  const int   T   = E.n_cols;
  const int   N   = E.n_rows;
  const int   M   = aux_PR_TR.n_rows;
  mat aux_xi_tmp = aux_xi;
  
  cube  Z(N, T, M);
  for (int m=0; m<M; m++) {
    Z.slice(m)    = pow(aux_sigma, -1) % (aux_B.slice(m) * E);
  }
  
  mat filtered    = filtering(Z, aux_PR_TR, aux_pi_0);
  mat smoothed    = smoothing(filtered, aux_PR_TR);
  mat    aj       = eye(M, M);
  
  mat xi(M, T);
  int draw        = bsvars::csample_num1(wrap(seq_len(M)), wrap(smoothed.col(T-1)));
  aux_xi_tmp.col(T-1)     = aj.col(draw-1);
  
  if ( minimum_regime_occurrences==0 ) {
    for (int t=T-2; t>=0; --t) {
      vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
      draw          = bsvars::csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
      aux_xi_tmp.col(t)   = aj.col(draw-1);
    }
    aux_xi = aux_xi_tmp;
  } else {
    int regime_occurrences  = 1;
    int iterations  = 1;
    while ( (regime_occurrences<minimum_regime_occurrences) & (iterations<max_iterations) ) {
      for (int t=T-2; t>=0; --t) {
        vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
        draw          = bsvars::csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
        aux_xi_tmp.col(t)   = aj.col(draw-1);
      }
      mat transitions       = bsvars::count_regime_transitions(aux_xi_tmp);
      regime_occurrences    = min(transitions.diag());
      iterations++;
    } // END while
    if ( iterations<max_iterations ) aux_xi = aux_xi_tmp;
  }
} // END sample_Markov_process_mss



