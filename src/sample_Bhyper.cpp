
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

#include "bsvars.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
void sample_B_heterosk1_s4 (
    arma::mat&                    aux_B,          // NxN
    arma::ivec&                   aux_SL,         // Nx1 row-specific S4 indicators
    const arma::mat&              aux_A,          // NxK
    const arma::vec&              aux_hyper,      // NxM
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              Y,              // NxT dependent variables
    const arma::mat&              X,              // KxT dependent variables
    const Rcpp::List&              prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL       // restrictions on B0 in S4 arrangement
) {
  // this function replaces sample_B_heterosk1
  // the function changes the value of aux_B and aux_SL by reference
  
  const int N           = aux_B.n_rows;
  const int T           = Y.n_cols;
  
  int         Ltmp      = VBL.n_elem - 1;
  vec         Lm        = VBL(Ltmp);
  double      L         = accu(Lm);
  field<mat>  VB        = VBL.rows(0, L-1);
  
  const int posterior_nu    = T + N;
  mat prior_SS_inv          = pow(aux_hyper(0), -1) * as<mat>(prior["B_V_inv"]);
  mat shocks                = Y - aux_A * X;
  
  for (int n=0; n<N; n++) {
    mat aux_B_nL(Lm(n), N);
    vec log_posterior_kernel_nL(Lm(n));
    mat aux_B_tmp           = aux_B;
    
    for (int l=0; l<Lm(n); l++) {
      int ll = 0;
      if (n == 0) {
        ll                    = l;
      } else {
        vec Lm_cs             = cumsum(Lm);
        ll                    = Lm_cs(n-1) + l;
      }
      
      // set scale matrix
      mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
      mat posterior_SS_inv    = prior_SS_inv + shocks_sigma * shocks_sigma.t();
      mat posterior_S_inv     = VB(ll) * posterior_SS_inv * VB(ll).t();
      posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
      
      // sample B
      mat Un                  = chol(posterior_nu * inv_sympd(posterior_S_inv));
      mat B_tmp               = aux_B;
      B_tmp.shed_row(n);
      rowvec w                = trans(bsvars::orthogonal_complement_matrix_TW(B_tmp.t()));
      vec w1_tmp              = trans(w * VB(ll).t() * Un.t());
      double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
      mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
      mat Wn;
      int rn                  = VB(ll).n_rows;
      if (rn==1) {
        Wn                    = w1;
      } else {
        Wn                    = join_rows(w1.t(), bsvars::orthogonal_complement_matrix_TW(w1.t()));
      }
      
      vec   alpha(rn);
      vec   u                 = rnorm(posterior_nu+1, 0, pow(posterior_nu, -0.5));
      alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
      if (R::runif(0,1)<0.5) {
        alpha(0)             *= -1;
      }
      if (rn>1){
        vec nn                = Rcpp::rnorm(rn-1, 0, pow(posterior_nu, -0.5));
        alpha.rows(1,rn-1)    = nn;
      }
      rowvec b0n              = alpha.t() * Wn * Un;
      aux_B_nL.row(l)         = b0n * VB(ll);
      
      // posterior kernel
      aux_B_tmp.row(n)            = aux_B_nL.row(l);
      mat std_shocks              = (aux_B_tmp * shocks) / aux_sigma;
      log_posterior_kernel_nL(l) -= 0.5 * accu(square(std_shocks));                   // likelihood kernel exp part
      double abs_log_det_B;
      double sign_;
      log_det(abs_log_det_B, sign_, aux_B_tmp); 
      log_posterior_kernel_nL(l) += T * abs_log_det_B;                                // likelihood kernel det part
      log_posterior_kernel_nL(l) -= 0.5 * as_scalar(b0n * posterior_S_inv * b0n.t()); // prior B kernel
      log_posterior_kernel_nL(l) += log(1/Lm(n));                                     // prior multinomial - it's flat, so this does not matter
    } // END loop l
    
    // Sample S4 indicator
    int     index_s4          = 0;
    if (Lm(n) > 1) {
      // Compute S4 components probabilities
      log_posterior_kernel_nL -= log_posterior_kernel_nL.max();
      vec     pr_s4           = exp(log_posterior_kernel_nL)/accu(exp(log_posterior_kernel_nL));
      
      // Sample S4 indicator
      NumericVector seq_1S    = wrap(seq_len(Lm(n)) - 1);
      NumericVector indi_tmp  = Rcpp::RcppArmadillo::sample(seq_1S, 1, false, wrap(pr_s4));
      index_s4                = indi_tmp(0);
    }
    aux_SL(n)                 = index_s4;
    aux_B.row(n)              = aux_B_nL.row(index_s4);
    
  } // END n loop
} // END sample_B_heterosk1_s4




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
void sample_hyperparameters_s4 (
    arma::vec&              aux_hyper,
    const arma::mat&        aux_B,
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
) {
  // the function changes the value of aux_hyper by reference (filling it with a new draw)
  const int N = aux_B.n_rows;
  const int K = aux_A.n_cols;
  
  int         Ltmp = VB.n_elem - 1;
  vec         Lm  = VB(Ltmp);
  vec         Lm_cs = cumsum(Lm);
  int rn=0;
  for (int n=0; n<N; n++) {
    int ll        = aux_SL(n);
    if (n>0) {
      ll         += Lm_cs(n-1);
    }
    rn           += VB(ll).n_rows;
  }
  
  // aux_hyper(4)    = ( as<double>(prior["hyper_S"]) + aux_hyper(2) + aux_hyper(3) ) / R::rchisq(as<double>(prior["hyper_V"]) + 4*as<double>(prior["hyper_a"]));
  // aux_hyper(3)    = R::rgamma( as<double>(prior["hyper_a"]) + 0.5*as<double>(prior["hyper_nu"]) ,
  //           1/((1/aux_hyper(4)) + (1/(2*aux_hyper(1)))));
  // aux_hyper(2)    = R::rgamma( as<double>(prior["hyper_a"]) + 0.5*as<double>(prior["hyper_nu"]) ,
  //           1/((1/aux_hyper(4)) + (1/(2*aux_hyper(0)))) );
  aux_hyper(1)    = ( aux_hyper(3) + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
    R::rchisq( as<double>(prior["hyper_nu"]) + N * K );
  aux_hyper(0)    = ( aux_hyper(2) + trace(aux_B * as<mat>(prior["B_V_inv"]) * trans(aux_B) )) /
    R::rchisq( as<double>(prior["hyper_nu"]) + rn );
} // END sample_hyperparameters_s4

