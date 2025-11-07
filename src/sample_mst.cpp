
#include <RcppArmadillo.h>
#include <RcppTN.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_lambda_ms (
    const arma::mat&    aux_df,     // NxM
    const arma::mat&    aux_xi      // MxT
) {
  const int N           = aux_df.n_rows;
  const int T           = aux_xi.n_cols;
  vec       Tm          = sum(aux_xi, 1);  
  
  mat       nu_lambda   = aux_df + 1;
  mat       s_lambda    = aux_df - 1;
  mat       aux_lambda(N, T);
  
  for (int n=0; n<N; n++) {
    for (int t=0; t<T; t++) {
  
      int     m         = aux_xi.col(t).index_max();
      double  draw      = as<double>(Rcpp::rchisq(1, nu_lambda(n, m)));
      aux_lambda(n, t)  = s_lambda(n, t) / draw;
  
    }
  }
  
  return aux_lambda;
} // END sample_lambda_ms


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_kernel_df_ms_nm (
    const double&         aux_df,
    const arma::rowvec&   aux_lambda,  // Tmx1
    const double&         prior_df_a
) {
  
  const int Tm  = aux_lambda.n_elem;
  double lk_df  = 0;
  if ( Tm != 0 ) {
    lk_df   -= Tm * lgamma(0.5 * aux_df);                        // lambda prior
    lk_df   += 0.5 * Tm * aux_df * log(0.5 * (aux_df - 2));      // lambda prior
    lk_df   -= 0.5 * (aux_df - 2) * accu(log(aux_lambda));       // lambda prior
    lk_df   -= 0.5 * (aux_df - 2) * accu(pow(aux_lambda, -1));  // lambda prior
  }
  lk_df   -= prior_df_a * (aux_df - 2);                       // df prior shifted exponential
  
  return lk_df;
} // END log_kernel_df_ms_nm


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_df_ms (
    arma::mat&        aux_df,             // NxM
    const arma::mat&  aux_lambda,         // NxT
    const arma::mat&  aux_xi,             // MxT
    const arma::mat&  U,                  // NxT
    const Rcpp::List& prior,              // hyper-parameter for exponential prior for aux_df
    const int&        s,                  // MCMC iteration
    arma::mat&        adaptive_scale,     // NxM
    const arma::vec&  adptive_alpha_gamma // 2x1 vector with target acceptance rate and step size
) {
  
  int N = aux_df.n_rows;
  int M = aux_df.n_cols;
  mat aux_df_star(N, M);
  mat alpha(N, M, fill::ones);
  
  double  prior_df_a          = as<double>(prior["df_a"]);
  
  // by sampling from truncated normal it is assumed that the asymmetry from truncation 
  // is negligible for alpha computation
  for (int n = 0; n < N; n++){
    for (int m = 0; m < M; m++){
      aux_df_star(n, m)       = RcppTN::rtn1( aux_df(n, m), adaptive_scale(n), 2, R_PosInf );
      
      uvec    aux_xi_m        = find(aux_xi.row(m));
      rowvec  lambda_tmp      = aux_lambda.row(n);
      rowvec  lambda          = lambda_tmp.cols(aux_xi_m);
      
      double  kernel_ratio     = exp( log_kernel_df_ms_nm(aux_df_star(n, m), lambda, prior_df_a) - log_kernel_df_ms_nm(aux_df(n, m), lambda, prior_df_a) );
      kernel_ratio            *= RcppTN::dtn1(aux_df_star(n, m), aux_df(n, m), adaptive_scale(n), 2, R_PosInf) / RcppTN::dtn1(aux_df(n, m), aux_df_star(n, m), adaptive_scale(n), 2, R_PosInf);
             
      if ( kernel_ratio < 1 ) alpha(n, m) = kernel_ratio;
      if ( R::runif(0, 1) < alpha(n, m) ) {
        aux_df(n, m) = aux_df_star(n, m);
      }
      if (s > 1) {
        adaptive_scale(n, m) = exp( log(adaptive_scale(n, m)) + 0.5 * log( 1 + pow(s, - adptive_alpha_gamma(1)) * (alpha(n, m) - adptive_alpha_gamma(0))) );
      }
    } // END m loop
  } // END n loop
  
  return List::create(
    _["aux_df"] = aux_df,
    _["adaptive_scale"] = adaptive_scale
  );
} // END sample_df_ms
