
#include <RcppArmadillo.h>
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarTVPs_ir (
    arma::field<arma::cube>&  posterior_B,        // (S)(N, N, M)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p
) {
  
  const int       N = posterior_B(0).n_rows;
  const int       M = posterior_B(0).n_slices;
  const int       S = posterior_B.n_elem;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S, M);
  
  for (int s=0; s<S; s++) {
    Rcout << "Iteration: " << s << endl;
    for (int m=0; m<M; m++) {
      Rcout << "  Regime: " << m << endl;
      aux_irfs            = bsvars::bsvars_ir1( posterior_B(s).slice(m), posterior_A.slice(s), horizon, p );
      irfs(s, m)          = aux_irfs;
    }
  } // END s loop
  
  return irfs;
} // END bsvarTVPs_ir

