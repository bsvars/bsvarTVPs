
#include <RcppArmadillo.h>
#include <RcppTN.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_lambda_ms (
    const arma::mat&    aux_df,     // NxM
    const arma::mat&    aux_xi,     // MxT
    arma::mat&          U,          // NxT
    const arma::uvec    studentt    // Nx1, {FALSE - normal, TRUE - Student-t};
) {
  
  const int N           = aux_df.n_rows;
  const int T           = aux_xi.n_cols;
  
  U                     = square( U );
  U.each_col()         /= sum(U, 1) / T;        // normalisation E[u] = 1
  
  mat       nu_lambda   = aux_df + 1;
  mat       s_lambda    = U - 2; // + aux_df(n,m) done below
  mat       aux_lambda(N, T, fill::ones);
  
  for (int n=0; n<N; n++) {
    if (studentt(n) == 0) continue;         // skip normal errors
    for (int t=0; t<T; t++) {
      int     m         = aux_xi.col(t).index_max();
      double  draw      =  chi2rnd(nu_lambda(n, m));
      // lambda | U, df ~ IG((df+1)/2, (df-2+U^2)/2), with shape/scale parameters.
      aux_lambda(n, t)  = (aux_df(n, m) + s_lambda(n, t)) / draw;
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
    lk_df   -= 0.5 * (aux_df + 2) * accu(log(aux_lambda));       // lambda prior
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
    const arma::vec&  adptive_alpha_gamma,// 2x1 vector with target acceptance rate and step size
    const arma::uvec  studentt            // Nx1, {FALSE - normal, TRUE - Student-t};
) {
  
  int N = aux_df.n_rows;
  int M = aux_df.n_cols;
  mat aux_df_star(N, M);
  mat alpha(N, M, fill::ones);
  
  double  prior_df_a          = as<double>(prior["df_a"]);
  
  // Account for the asymmetry of the normal proposal truncated below at df=2.
  for (int n = 0; n < N; n++){
    if (studentt(n) == 0) continue;         // skip normal errors
    for (int m = 0; m < M; m++){
      
      const double proposal_sd = adaptive_scale(n, m);
      aux_df_star(n, m)       = RcppTN::rtn1( aux_df(n, m), proposal_sd, 2, R_PosInf );
      
      uvec    aux_xi_m        = find(aux_xi.row(m));
      rowvec  lambda_tmp      = aux_lambda.row(n);
      rowvec  lambda          = lambda_tmp.cols(aux_xi_m);
      
      double log_ratio = log_kernel_df_ms_nm(aux_df_star(n, m), lambda, prior_df_a)
        - log_kernel_df_ms_nm(aux_df(n, m), lambda, prior_df_a);
      // q(old | new) / q(new | old) = Phi((old-2)/sd) / Phi((new-2)/sd).
      log_ratio += R::pnorm5((aux_df(n, m)-2)/proposal_sd, 0, 1, true, true)
        - R::pnorm5((aux_df_star(n, m)-2)/proposal_sd, 0, 1, true, true);
      
      if (log_ratio < 0) alpha(n, m) = exp(log_ratio);
      if (std::log(R::runif(0, 1)) < log_ratio) {
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
