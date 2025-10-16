
#include <RcppArmadillo.h>
#include "progress.hpp"

#include "sample_ABhyper.h"
#include "sample_sv_ms.h"
#include "sample_mst.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_mss_tvi_sv_cpp (
    const int&                    SS,                   // No. of posterior draws
    const arma::mat&              Y,                    // NxT dependent variables
    const arma::mat&              X,                    // KxT explanatory variables
    const Rcpp::List&             prior,                // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,                   // restrictions on B0
    const Rcpp::List&             starting_values,
    const int                     thin = 100,           // introduce thinning
    const int                     sv_select = 1,        // {1 - non-centred, 2 - centred, 3 - homoskedastic};
    const int                     hyper_select = 1,     // {1 - horseshoe, 2 - boost, 3 - fixed}
    const bool                    finiteM = true,       // {true - stationary MS, false - overfitted};
    const bool                    studentt = false      // {true - normal, false - Student-t};
) {
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  Rcout << "**************************************************|" << endl;
  Rcout << "bsvarTVPs: Bayesian Structural VARs with          |" << endl;
  Rcout << "  Markov-Switching, Time-Varying Identification   |" << endl;
  Rcout << "  and Stochastic Volatility                       |" << endl;
  Rcout << "**************************************************|" << endl;
  Rcout << " Gibbs sampler for the SVAR-SV model              |" << endl;
  Rcout << "    with Markov-switching and regime-specific     |" << endl;
  Rcout << "    time-varying identification                   |" << endl;
  Rcout << "    for the structural matrix                     |" << endl;
  if ( studentt ) {
  Rcout << "    and Student-t structural shocks               |" << endl;
  }
  Rcout << "**************************************************|" << endl;
  Rcout << " Progress of the MCMC simulation for " << SS << " draws" << endl;
  Rcout << "    Every " << thin << "th draw is saved via MCMC thinning" << endl;
  Rcout << " Press Esc to interrupt the computations" << endl;
  Rcout << "**************************************************|" << endl;
  Progress p(50, true);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  const int   K     = X.n_rows;
  
  // aux_lambda - contains latent process
  // aux_lambda_sqrt = sqrt(aux_lambda)
  // aux_sigma - contains time-varying sds
  // aux_hetero - contains the diagonal of the covariance of conditional normal aux_hetero = aux_sigma * aux_lambda_sqrt
  
  cube  aux_B       = as<cube>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  List  aux_hyper   = as<List>(starting_values["hyper"]);  // (2*N+1)x2 (gamma_0, gamma_+, s_0, s_+, s_)
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  mat   aux_h       = as<mat>(starting_values["h"]);
  vec   aux_rho     = as<vec>(starting_values["rho"]);
  mat   aux_omega   = as<mat>(starting_values["omega"]);
  mat   aux_sigma2v = as<mat>(starting_values["sigma2v"]);
  umat  aux_S       = as<umat>(starting_values["S"]);
  vec   aux_sigma2_omega = as<vec>(starting_values["sigma2_omega"]);
  vec   aux_s_      = as<vec>(starting_values["s_"]);
  imat  aux_SL      = as<imat>(starting_values["S4_indicator"]) - 1;      // NxM S4 indicator matrix
  mat   aux_sigma   = as<mat>(starting_values["sigma"]);
  mat   aux_lambda  = as<mat>(starting_values["lambda"]);
  mat   aux_lambda_sqrt = sqrt(aux_lambda);
  mat   aux_df      = as<mat>(starting_values["df"]);
  mat   aux_hetero  = aux_sigma % aux_lambda_sqrt;
  
  const int M       = aux_PR_TR.n_cols;
  vec      Tm       = sum(aux_xi, 1);
  
  // parameters for adaptive sampling of degrees of freedom
  NumericVector aag = {0.44, 0.6};
  const vec     adptive_alpha_gamma = as<vec>(aag);
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for df evaluated at df = 30
  double  adaptive_scale_init = pow(R::psigamma(15, 1) - 29 * pow(28, -2), -1) / (T / M);
  mat     adaptive_scale(N, M, fill::value(adaptive_scale_init));
  
  const int   S     = floor(SS / thin);
  
  field<cube> posterior_B(S);
  cube  posterior_A(N, K, S);
  List  posterior_hyper(S);
  cube  posterior_PR_TR(M, M, S);
  mat   posterior_pi_0(M,S);
  cube  posterior_xi(M, T, S);
  cube  posterior_h(N, T, S);
  mat   posterior_rho(N, S);
  cube  posterior_omega(N, M, S);
  cube  posterior_sigma2v(N, M, S);
  ucube posterior_S(N, T, S);
  mat   posterior_sigma2_omega(N, S);
  mat   posterior_s_(N, S);
  icube posterior_SL(N, M, S);
  cube  posterior_lambda(N, T, S);
  cube  posterior_df(N, M, S);
  cube  posterior_sigma(N, T, S);
  
  vec   acceptance_count(4 + N);
  List  BSL = List::create(
    _["aux_B"]      = aux_B,
    _["aux_SL"]     = aux_SL
  );
  List  sv_n;
  List  PR_TR_tmp;
  
  field<mat> precisionB;
  field<mat> precisionA;
  
  mat aux_sigma_tmp_m(N, T, fill::ones);
  mat sigmaT(N, T, fill::ones);
  rowvec omega_T_n(T);
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    Rcout<<" s: "<<s<<endl;
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_lambda and aux_df
    Rcout<<" sample aux_lambda and aux_df"<<endl;
    mat E             = (Y - aux_A * X);
    if ( studentt ) {
      
      mat U           = E;
      for (int t=0; t<T; t++) {
        int   m       = aux_xi.col(t).index_max();
        U.col(t)      = aux_B.slice(m) * (Y.col(t) - aux_A * X.col(t));
      }
      U              /= aux_sigma;
      
      aux_lambda      = sample_lambda_ms(aux_df, aux_xi, U);
      aux_lambda_sqrt = sqrt(aux_lambda);
      aux_hetero      = aux_sigma % aux_lambda_sqrt;
      
      List aux_df_tmp = sample_df_ms (aux_df, aux_lambda, aux_xi, U, prior, ss, adaptive_scale, adptive_alpha_gamma);
      aux_df          = as<mat>(aux_df_tmp["aux_df"]);
      adaptive_scale  = as<mat>(aux_df_tmp["adaptive_scale"]);
    } // END studentt
    
    // sample aux_xi
    Rcout<<" sample aux_xi"<<endl;
    cube Z(N, T, M);
    for (int m=0; m<M; m++) {
      Z.slice(m)        = aux_B.slice(m) * (Y - aux_A * X);
      if ( sv_select != 3 ) {
        aux_sigma_tmp_m = exp(0.5 * diagmat(aux_omega.col(m)) * aux_h);
        Z.slice(m)     /= aux_sigma_tmp_m;
      }
    }
    if ( studentt ) {
      aux_xi          = sample_Markov_process_studentt(Z, aux_xi, aux_PR_TR, aux_pi_0, aux_df, finiteM);
    } else {
      aux_xi          = sample_Markov_process(Z, aux_xi, aux_PR_TR, aux_pi_0, finiteM);
    }    
    
    // sample aux_PR_TR and aux_pi_0
    Rcout<<" sample aux_PR_TR and aux_pi_0"<<endl;
    PR_TR_tmp         = sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior);
    aux_PR_TR         = as<mat>(PR_TR_tmp["aux_PR_TR"]);
    aux_pi_0          = as<vec>(PR_TR_tmp["aux_pi_0"]);
    
    // sample aux_hyper
    Rcout<<" sample aux_hyper"<<endl;
    if ( hyper_select == 1 ) {

      aux_hyper       = sample_hyperparameter_mss_s4_horseshoe(aux_hyper, aux_B, aux_A, VB, aux_SL, prior);
      precisionB      = hyper2precisionB_mss_horseshoe(aux_hyper);
      precisionA      = hyper2precisionA_horseshoe(aux_hyper);

    } else if ( hyper_select == 2 ) {

      aux_hyper       = sample_hyperparameters_mss_s4_boost( aux_hyper, aux_B, aux_A, VB, aux_SL, prior, true);
      precisionB      = hyper2precisionB_mss_boost(aux_hyper, prior);
      precisionA      = hyper2precisionA_boost(aux_hyper, prior);

    } else if ( hyper_select == 3 ) {

      aux_hyper       = sample_hyperparameters_mss_s4_boost( aux_hyper, aux_B, aux_A, VB, aux_SL, prior, false);
      precisionB      = hyper2precisionB_mss_boost(aux_hyper, prior);
      precisionA      = hyper2precisionA_boost(aux_hyper, prior);

    }
    
    // sample aux_B
    Rcout<<" sample aux_B"<<endl;
    BSL               = sample_B_mss_s4(aux_B, aux_SL, aux_A, precisionB, aux_hetero, aux_xi, Y, X, prior, VB);
    aux_B             = as<cube>(BSL["aux_B"]);
    aux_SL            = as<imat>(BSL["aux_SL"]);
    
    // sample aux_A
    Rcout<<" sample aux_A"<<endl;
    aux_A             = sample_A_heterosk1_mss(aux_A, aux_B, aux_xi, precisionA, aux_hetero, Y, X, prior);
    
    // sample aux_h, aux_omega and aux_S, aux_sigma2_omega
    Rcout<<" sample aux_h, aux_omega and aux_S, aux_sigma2_omega"<<endl;
    if ( sv_select != 3 ) {
      
      mat U(N, T);
      E = Y - aux_A * X;
      for (int m=0; m<M; m++) {
        for (int t=0; t<T; t++) {
          if (aux_xi(m,t)==1) {
            U.col(t)      = aux_B.slice(m) * E.col(t) / aux_lambda_sqrt.col(t);
          }
        }
      }
      
      for (int n=0; n<N; n++) {
        rowvec  h_tmp       = aux_h.row(n);
        double  rho_tmp     = aux_rho(n);
        rowvec  omega_tmp   = aux_omega.row(n);
        rowvec  sigma2v_tmp = square(aux_omega.row(n));
        urowvec S_tmp       = aux_S.row(n);
        rowvec  U_tmp       = U.row(n);
        double  s2o_tmp     = aux_sigma2_omega(n);
        double  s_n         = aux_s_(n);
        
        if ( sv_select == 2 ) {
          try {
            sv_n            = svar_ce1_mss( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, aux_xi, U_tmp, prior);
          } catch (std::runtime_error &e) {}
        } else if ( sv_select == 1 ) {
          sv_n              = svar_nc1_mss( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, aux_xi, U_tmp, prior);
        }
        
        aux_h.row(n)        = as<rowvec>(sv_n["aux_h_n"]);
        aux_rho(n)          = as<double>(sv_n["aux_rho_n"]);
        aux_omega.row(n)    = as<rowvec>(sv_n["aux_omega_n"]);
        aux_sigma2v.row(n)  = as<rowvec>(sv_n["aux_sigma2v_n"]);
        aux_S.row(n)        = as<urowvec>(sv_n["aux_S_n"]);
        aux_sigma2_omega(n) = as<double>(sv_n["aux_sigma2_omega_n"]);
        aux_s_(n)           = as<double>(sv_n["aux_s_n"]);
        
        for (int t=0; t<T; t++){
          omega_T_n(t)      = aux_omega(n, aux_xi.col(t).index_max());
        }
        aux_sigma.row(n)    = exp(0.5 * (aux_h.row(n) % omega_T_n));
        
      } // END n loop
      aux_hetero            = aux_sigma % aux_lambda_sqrt;
      
    } // END if( sv_select != 3 )
    
    Rcout<<" save in posterior"<<endl;
    if (ss % thin == 0) {
      posterior_B(s)                = aux_B;
      posterior_A.slice(s)          = aux_A;
      posterior_hyper(s)            = aux_hyper;
      posterior_PR_TR.slice(s)      = aux_PR_TR;
      posterior_pi_0.col(s)         = aux_pi_0;
      posterior_xi.slice(s)         = aux_xi;
      posterior_lambda.slice(s)     = aux_lambda;
      posterior_df.slice(s)         = aux_df;
      posterior_h.slice(s)          = aux_h;
      posterior_rho.col(s)          = aux_rho;
      posterior_omega.slice(s)      = aux_omega;
      posterior_sigma2v.slice(s)    = aux_sigma2v;
      posterior_S.slice(s)          = aux_S;
      posterior_sigma2_omega.col(s) = aux_sigma2_omega;
      posterior_s_.col(s)           = aux_s_;
      posterior_SL.slice(s)         = aux_SL;
      posterior_sigma.slice(s)      = aux_sigma;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["PR_TR"]    = aux_PR_TR,
      _["pi_0"]     = aux_pi_0,
      _["xi"]       = aux_xi,
      _["h"]        = aux_h,
      _["rho"]      = aux_rho,
      _["omega"]    = aux_omega,
      _["sigma2v"]  = aux_sigma2v,
      _["S"]        = aux_S,
      _["sigma2_omega"] = aux_sigma2_omega,
      _["s_"]       = aux_s_,
      _["S4_indicator"] = aux_SL + 1,
      _["sigma"]    = aux_sigma,
      _["lambda"]   = aux_lambda,
      _["df"]       = aux_df
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["PR_TR"]    = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi"]       = posterior_xi,
      _["h"]        = posterior_h,
      _["rho"]      = posterior_rho,
      _["omega"]    = posterior_omega,
      _["sigma2v"]  = posterior_sigma2v,
      _["S"]        = posterior_S,
      _["sigma2_omega"] = posterior_sigma2_omega,
      _["s_"]        = posterior_s_,
      _["S4_indicator"] = posterior_SL + 1,
      _["sigma"]    = posterior_sigma,
      _["lambda"]   = posterior_lambda,
      _["df"]       = posterior_df
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_mss_tvi_sv_cpp
