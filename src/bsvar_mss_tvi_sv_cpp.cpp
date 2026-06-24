
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "bsvars.h"

#include "sample_ABhyper.h"
#include "sampleTheta0.h"
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
    const Rcpp::List&             VTheta0,              // restrictions on Theta0
    const arma::field<arma::mat>& VA,                   // restrictions on A
    const Rcpp::List&             starting_values,
    const arma::uvec              sv_select,            // Nx1, for each equation: {1 - non-centred, 2 - centred, 3 - homoskedastic};
    const arma::uvec              studentt,             // Nx1, {FALSE - normal, TRUE - Student-t};
    const int                     thin = 100,           // introduce thinning
    const int                     hyper_select = 3,     // {1 - horseshoe, 2 - boost, 3 - fixed}
    const bool                    finiteM = true,       // {true - stationary MS, false - overfitted};
    const bool                    fixed_regime = false, // {true - don't estimate MS, false - estimate MS};
    const bool                    show_progress = true
) {
  
  const bool debug = false;
  if ( debug ) Rcout << " start!" << endl;
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, SS, 50));
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = bsvars::ordinal(thin) + " ";
  }
  
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << " bsvarTVPs: Bayesian Structural VARs              |" << endl;
    Rcout << "            with Time-Varying Identification      |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << SS << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, true);
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  const int   K     = X.n_rows;
  
  // aux_lambda - contains latent process
  // aux_lambda_sqrt = sqrt(aux_lambda)
  // aux_sigma - contains time-varying sds
  // aux_hetero - contains the diagonal of the covariance of conditional normal aux_hetero = aux_sigma * aux_lambda_sqrt
  cube  aux_B       = as<cube>(starting_values["B"]);
  cube  aux_Theta0  = as<cube>(starting_values["Theta0"]);
  cube  aux_Theta0_inv(size(aux_Theta0));
  cube  aux_struc(size(aux_B));
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
  icube aux_SL      = as<icube>(starting_values["S4_indicator"]) - 1;      // NxMx2 S4 indicator matrix
  mat   aux_sigma   = as<mat>(starting_values["sigma"]);
  mat   aux_lambda  = as<mat>(starting_values["lambda"]);
  mat   aux_lambda_sqrt = sqrt(aux_lambda);
  mat   aux_df      = as<mat>(starting_values["df"]);
  mat   aux_hetero  = aux_sigma % aux_lambda_sqrt;
  
  const int M       = aux_PR_TR.n_cols;
  vec      Tm       = sum(aux_xi, 1);
  for (int m=0; m<M; m++) {
    aux_Theta0_inv.slice(m) = inv(aux_Theta0.slice(m));
    aux_struc.slice(m)      = aux_Theta0.slice(m) * aux_B.slice(m);
  }
  
  // parameters for adaptive sampling of degrees of freedom
  NumericVector aag = {0.44, 0.6};
  const vec     adptive_alpha_gamma = as<vec>(aag);
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for df evaluated at df = 30
  double  adaptive_scale_init = abs(pow(R::psigamma(15, 1)/4 - 29 * pow(28, -2), -1) / (T / M));
  mat     adaptive_scale(N, M, fill::value(adaptive_scale_init));
  
  const int   S     = floor(SS / thin);
  
  field<cube> posterior_B(S);
  field<cube> posterior_Theta0(S);
  field<cube> posterior_struc(S);
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
  field<icube> posterior_SL(S);
  cube  posterior_lambda(N, T, S);
  cube  posterior_df(N, M, S);
  cube  posterior_sigma(N, T, S);
  
  vec   acceptance_count(4 + N);
  List  BSL = List::create(
    _["aux_B"]      = aux_B,
    _["aux_SL"]     = aux_SL
  );
  List  TSL = List::create(
    _["aux_Theta0"] = aux_Theta0,
    _["aux_SL"]     = aux_SL
  );
  List  sv_n = List::create(
    _["aux_h_n"]              = aux_h.row(0),
    _["aux_rho_n"]            = aux_rho(0),
    _["aux_omega_n"]          = aux_omega.row(0),
    _["aux_sigma2v_n"]        = square(aux_omega.row(0)),
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega(0),
    _["aux_s_n"]              = aux_s_(0),
    _["aux_S_n"]              = aux_S.row(0)
  );
  List  PR_TR_tmp;
  
  field<mat> precisionB;
  field<mat> precisionA;
  
  mat aux_sigma_tmp_m(N, T, fill::ones);
  mat sigmaT(N, T, fill::ones);
  rowvec omega_T_n(T);
  List aux_df_tmp;
  
  int   s = 0;
  
  for (int ss=0; ss<SS; ss++) {
    
    if ( debug ) Rcout<<" s: "<< s <<endl;
    
    // Increment progress bar
    if (any(prog_rep_points == ss)) p.increment();
    // Check for user interrupts
    if (ss % 200 == 0) checkUserInterrupt();
    
    // sample aux_lambda and aux_df
    if ( debug ) Rcout<<" sample aux_lambda and aux_df"<<endl;
    mat E             = (Y - aux_A * X);
    
    mat U           = E;
    for (int t=0; t<T; t++) {
      int   m       = aux_xi.col(t).index_max();
      U.col(t)      = aux_struc.slice(m) * (Y.col(t) - aux_A * X.col(t));
    }
    U              /= aux_sigma;
    
    try {
      aux_lambda      = sample_lambda_ms(aux_df, aux_xi, U, studentt);
      aux_lambda_sqrt = sqrt(aux_lambda);
      aux_hetero      = aux_sigma % aux_lambda_sqrt;
    } catch (std::runtime_error &e) {}
    
    try {
      aux_df_tmp      = sample_df_ms (aux_df, aux_lambda, aux_xi, U, prior, ss, adaptive_scale, adptive_alpha_gamma, studentt);
      aux_df          = as<mat>(aux_df_tmp["aux_df"]);
      adaptive_scale  = as<mat>(aux_df_tmp["adaptive_scale"]);
    } catch (std::runtime_error &e) {}
    
    // sample aux_xi
    if ( debug ) Rcout<<" sample aux_xi"<<endl;
    if (!fixed_regime) {
      cube Z(N, T, M);
      for (int m=0; m<M; m++) {
        Z.slice(m)     = aux_struc.slice(m) * (Y - aux_A * X);
        for (int n=0; n<N; n++) {
          if ( sv_select(n) != 3 ) {
            aux_sigma_tmp_m.row(n) = exp(0.5 * aux_omega(n,m) * aux_h.row(n));
            Z.slice(m).row(n)     /= aux_sigma_tmp_m.row(n);
          }
        }
      } // END loop m
      if ( all(studentt == 1) ) {
        try {
          aux_xi        = sample_Markov_process_studentt(Z, aux_xi, aux_PR_TR, aux_pi_0, aux_df, finiteM);
        } catch (std::runtime_error &e) {
          Rcout << "   sample_Markov_process_studentt failure " << endl;
        }
      } else {
        try {
          aux_xi        = sample_Markov_process(Z, aux_xi, aux_PR_TR, aux_pi_0, finiteM);
        } catch (std::runtime_error &e) {
          Rcout << "   sample_Markov_process failure " << endl;
        }
      } // END if ( studentt )
    } // END (!fixed_regime)
    
    // sample aux_PR_TR and aux_pi_0
    if ( debug ) Rcout<<" sample aux_PR_TR and aux_pi_0"<<endl;
    if (!fixed_regime) {
      try {
        PR_TR_tmp       = sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior);
        aux_PR_TR       = as<mat>(PR_TR_tmp["aux_PR_TR"]);
        aux_pi_0        = as<vec>(PR_TR_tmp["aux_pi_0"]);
      } catch (std::runtime_error &e) {}
    } // END (!fixed_regime)
    
    // sample aux_hyper
    if ( debug ) Rcout<<" sample aux_hyper"<<endl;
    if ( hyper_select == 1 ) {

      try {
        aux_hyper     = sample_hyperparameter_mss_s4_horseshoe(aux_hyper, aux_B, aux_A, VB, aux_SL.slice(0), prior);
        precisionB    = hyper2precisionB_mss_horseshoe(aux_hyper);
        precisionA    = hyper2precisionA_horseshoe(aux_hyper);
      } catch (std::logic_error &e) {
        Rcout << "   sample_hyperparameter_mss_s4_horseshoe failure " << endl;
      }

    } else if ( hyper_select == 2 ) {

      aux_hyper       = sample_hyperparameters_mss_s4_boost( aux_hyper, aux_B, aux_A, VB, aux_SL.slice(0), prior, true);
      precisionB      = hyper2precisionB_mss_boost(aux_hyper, prior);
      precisionA      = hyper2precisionA_boost(aux_hyper, prior);

    } else if ( hyper_select == 3 ) {

      aux_hyper       = sample_hyperparameters_mss_s4_boost( aux_hyper, aux_B, aux_A, VB, aux_SL.slice(0), prior, false);
      precisionB      = hyper2precisionB_mss_boost(aux_hyper, prior);
      precisionA      = hyper2precisionA_boost(aux_hyper, prior);

    }
    
    // sample aux_h, aux_omega and aux_S, aux_sigma2_omega
    if ( debug ) Rcout<<" sample aux_h, aux_omega and aux_S, aux_sigma2_omega"<<endl;
    mat UU(N, T);
    E = Y - aux_A * X;
    for (int m=0; m<M; m++) {
      for (int t=0; t<T; t++) {
        if (aux_xi(m,t)==1) {
          UU.col(t)      = aux_struc.slice(m) * E.col(t) / aux_lambda_sqrt.col(t);
        }
      }
    }
    
    for (int n=0; n<N; n++) {
      if ( sv_select(n) != 3 ) {
      if ( debug ) Rcout<<" sv n: "<< n << endl;
        rowvec  h_tmp       = aux_h.row(n);
        double  rho_tmp     = aux_rho(n);
        rowvec  omega_tmp   = aux_omega.row(n);
        rowvec  sigma2v_tmp = square(aux_omega.row(n));
        urowvec S_tmp       = aux_S.row(n);
        rowvec  U_tmp       = UU.row(n);
        double  s2o_tmp     = aux_sigma2_omega(n);
        double  s_n         = aux_s_(n);
        
        if ( sv_select(n) == 2 ) {
          try {
            sv_n            = svar_ce1_mss( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, aux_xi, U_tmp, prior);
          } 
          catch (std::runtime_error &e) {}
          catch (std::logic_error &e) {}
        } else if ( sv_select(n) == 1 ) {
          try {
            sv_n            = svar_nc1_mss( h_tmp, rho_tmp, omega_tmp, sigma2v_tmp, s2o_tmp, s_n, S_tmp, aux_xi, U_tmp, prior);
          } 
          catch (std::runtime_error &e) {}
          catch (std::logic_error &e) {}
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
        
      } // END if( sv_select(n) != 3 )
    } // END n loop
    aux_hetero            = aux_sigma % aux_lambda_sqrt;
    
    
    // sample aux_A
    if ( debug ) Rcout<<" sample aux_A"<<endl;
    try {
      aux_A           = sample_A_heterosk1_mss(aux_A, aux_struc, aux_xi, precisionA, aux_hetero, Y, X, prior, VA);
    } catch (std::runtime_error &e) {}
    
    
    // sample aux_B, aux_SL
    if ( debug ) Rcout<<" sample aux_B, aux_SL"<< endl;
    try {
      
      BSL             = sample_B_mss_s4 (aux_B, aux_SL.slice(0), aux_Theta0, aux_Theta0_inv, aux_A,
                                         precisionB, aux_hetero, aux_xi, Y, X, prior, VB);
      aux_B           = as<cube>(BSL["aux_B"]);
      aux_SL.slice(0) = as<imat>(BSL["aux_SL"]);
    } catch (std::runtime_error &e) {}
    
    // sample aux_Theta0, aux_SL
    if ( debug ) Rcout<<" sample aux_Theta0, aux_SL"<< endl;
    try {
      mat shocks      = Y - aux_A * X;
      TSL             = sample_Theta0_mss_s4( aux_Theta0, aux_SL.slice(1), aux_B, shocks, aux_hetero, aux_xi, prior, VTheta0 );
      aux_Theta0      = as<cube>(TSL["aux_Theta0"]);
      aux_SL.slice(1) = as<imat>(TSL["aux_SL"]);
      
      for (int m=0; m<M; m++) {
        aux_Theta0_inv.slice(m) = inv(aux_Theta0.slice(m));
        aux_struc.slice(m)      = aux_Theta0_inv.slice(m) * aux_B.slice(m);
      }
    } catch (std::runtime_error &e) {}
  
    
    if ( debug ) Rcout<<" save in posterior"<<endl;
    if (ss % thin == 0) {
      posterior_B(s)                = aux_B;
      posterior_Theta0(s)           = aux_Theta0;
      posterior_struc(s)            = aux_struc;  
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
      posterior_SL(s)               = aux_SL + 1;
      posterior_sigma.slice(s)      = aux_sigma;
      s++;
    }
  } // END ss loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["Theta0"]   = aux_Theta0,
      _["structural"] = aux_struc,
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
      _["B_cpp"]    = posterior_B,
      _["Theta0_cpp"] = posterior_Theta0,
      _["structural_cpp"] = posterior_struc,
      _["A_cpp"]    = posterior_A,
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
      _["S4_indicator"] = posterior_SL,
      _["sigma"]    = posterior_sigma,
      _["lambda"]   = posterior_lambda,
      _["df"]       = posterior_df
    ),
    _["acceptance_rate"] = 1- acceptance_count/SS
  );
} // END bsvar_mss_tvi_sv_cpp
