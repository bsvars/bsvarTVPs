
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"

#include "bsvars.h"

using namespace Rcpp;
using namespace arma;



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_B_heterosk1 (
    arma::mat         aux_B,          // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B0 and aux_Bplus by reference (filling it with a new draw)
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  mat prior_SS_inv          = pow(aux_hyper(0), -1) * as<mat>(prior["B_V_inv"]);
  mat shocks                = Y - aux_A * X;
  
  
  for (int n=0; n<N; n++) {
    
    // set scale matrix
    mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
    mat posterior_SS_inv    = prior_SS_inv + shocks_sigma * shocks_sigma.t();
    mat posterior_S_inv     = VB(n) * posterior_SS_inv * VB(n).t();
    posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
    
    // sample B
    mat Un                  = chol(posterior_nu * inv_sympd(posterior_S_inv));
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(bsvars::orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), bsvars::orthogonal_complement_matrix_TW(w1.t()));
    }
    
    vec   alpha(rn);
    vec   u(posterior_nu+1, fill::randn);
    u                      *= pow(posterior_nu, -0.5);
    alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn(rn-1, fill::randn);
      nn                   *= pow(posterior_nu, -0.5);
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)           = b0n * VB(n);
  } // END n loop
  
  return aux_B;
} // END sample_B_heterosk1



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube sample_B_mss (
    arma::cube        aux_B,          // NxNxM
    const arma::mat&  aux_A,          // NxK
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  const int T     = Y.n_cols;
  const int N     = Y.n_rows;
  const int K     = X.n_rows;
  const int M     = aux_B.n_slices;
  const vec T_m  = sum(aux_xi,1);
  
  // aux_sigma, Y, X
  
  for (int m=0; m<M; m++) {
    int ii = 0;
    mat aux_sigma_m(N, T_m(m));
    mat Y_m(N, T_m(m));
    mat X_m(K, T_m(m));
    
    for (int t=0; t<T; t++){
      if (aux_xi(m,t)==1) {
        aux_sigma_m.col(ii) = aux_sigma.col(t);
        Y_m.col(ii)         = Y.col(t);
        X_m.col(ii)         = X.col(t);
        ii++;
      }
    }
    aux_B.slice(m)    = sample_B_heterosk1(aux_B.slice(m), aux_A, aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
  }
  
  return aux_B;
} // END sample_B_mss



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_B_heterosk1_s4 (
    arma::mat                     aux_B,          // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators
    const arma::mat&              aux_A,          // NxK
    const arma::vec&              aux_hyper,      // NxM
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              Y,              // NxT dependent variables
    const arma::mat&              X,              // KxT dependent variables
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL       // restrictions on B0 in S4 arrangement
) {
  // the function draws new values of aux_B and aux_SL
  
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
  
  return List::create(
    _["aux_B"]    = aux_B,
    _["aux_SL"]   = aux_SL
  );
} // END sample_B_heterosk1_s4





/*______________________function sample_B_mss_s4______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_B_mss_s4 (
    arma::cube        aux_B,          // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    const arma::mat&  aux_A,          // NxK
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  const int T     = Y.n_cols;
  const int N     = Y.n_rows;
  const int K     = X.n_rows;
  const int M     = aux_B.n_slices;
  const vec T_m  = sum(aux_xi,1);
  
  // aux_sigma, Y, X
  
  for (int m=0; m<M; m++) {
    int ii = 0;
    mat aux_sigma_m(N, T_m(m));
    mat Y_m(N, T_m(m));
    mat X_m(K, T_m(m));
    
    for (int t=0; t<T; t++){
      if (aux_xi(m,t)==1) {
        aux_sigma_m.col(ii) = aux_sigma.col(t);
        Y_m.col(ii)         = Y.col(t);
        X_m.col(ii)         = X.col(t);
        ii++;
      }
    }
    // ivec aux_SL_m           = aux_SL.col(m);
    List BSL_m              = sample_B_heterosk1_s4(aux_B.slice(m), aux_SL.col(m), aux_A, aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
    aux_B.slice(m)          = as<mat>(BSL_m["aux_B"]);
    aux_SL.col(m)           = as<ivec>(BSL_m["aux_SL"]);
  }
  
  return List::create(
      _["aux_B"]    = aux_B,
      _["aux_SL"]   = aux_SL
  );
} // END sample_B_mss_s4





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_heterosk1 (
    arma::mat         aux_A,          // NxK
    const arma::mat&  aux_B,          // NxN
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior           // a list of priors - original dimensions
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_Vinv    = pow(aux_hyper(1), -1) * as<mat>(prior["A_V_inv"]);
  rowvec    zerosA(K);
  vec sigma_vectorised= vectorise(aux_sigma);
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B * (Y - A0 * X) );
    mat   zn_sigma    = zn / sigma_vectorised;
    mat   Wn          = kron( trans(X), aux_B.col(n) );
    mat   Wn_sigma    = Wn.each_col() / sigma_vectorised;
    
    mat     precision = prior_A_Vinv + trans(Wn_sigma) * Wn_sigma;
    precision         = 0.5 * (precision + precision.t());
    rowvec  location  = prior_A_mean.row(n) * prior_A_Vinv + trans(zn_sigma) * Wn_sigma;
    
    mat     precision_chol = trimatu(chol(precision));
    vec     xx(K, fill::randn);
    vec     draw      = solve(precision_chol, 
                              solve(trans(precision_chol), trans(location)) + xx);
    aux_A.row(n)      = trans(draw);
  } // END n loop
  
  return aux_A;
} // END sample_A_heterosk1





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_heterosk1_mss (
    arma::mat         aux_A,          // NxK
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    const arma::vec&  aux_hyper,      // NxM
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior           // a list of priors - original dimensions
) {
  // the function draws the value of aux_A
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  const int T         = Y.n_cols;
  const int M         = aux_xi.n_rows;
  const vec Tm       = sum(aux_xi, 1);
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_Vinv    = pow(aux_hyper(1), -1) * as<mat>(prior["A_V_inv"]);
  rowvec    zerosA(K);
  
  field<mat> YM(M);
  field<mat> XM(M);
  field<mat> SM(M);
  for (uword m=0; m<M; m++) {
    mat Y_m(N, Tm(m));
    mat X_m(K, Tm(m));
    mat S_m(N, Tm(m));
    
    YM(m)             = Y_m;
    XM(m)             = X_m;
    SM(m)             = S_m;
  }
  vec m_iter(M);
  
  for (uword t=0; t<T; t++) {
    int mm            = aux_xi.col(t).index_max();
    YM(mm).col(m_iter(mm))    = Y.col(t);
    XM(mm).col(m_iter(mm))    = X.col(t);
    SM(mm).col(m_iter(mm))    = aux_sigma.col(t);
    m_iter += aux_xi.col(t);
  }
  
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B.slice(0) * (YM(0) - A0 * XM(0)) );
    mat   Wn          = kron( trans(XM(0)), aux_B.slice(0).col(n) );
    vec   Sn          = vectorise(SM(0));
    
    for (uword m=1; m<M; m++) {
      zn          = join_cols(zn, vectorise(aux_B.slice(m) * (YM(m) - A0 * XM(m))));
      Wn          = join_cols(Wn, kron(trans(XM(m)), aux_B.slice(m).col(n)));
      Sn          = join_cols(Sn, vectorise(SM(m)));
    }
    
    mat   zn_sigma    = zn / Sn;
    mat   Wn_sigma    = Wn.each_col() / Sn;
    
    mat     precision = prior_A_Vinv + trans(Wn_sigma) * Wn_sigma;
    precision         = 0.5 * (precision + precision.t());
    rowvec  location  = prior_A_mean.row(n) * prior_A_Vinv + trans(zn_sigma) * Wn_sigma;
    
    mat     precision_chol = trimatu(chol(precision));
    vec     draw      = solve(precision_chol, 
                              solve(trans(precision_chol), trans(location)) + as<vec>(rnorm(K)));
    aux_A.row(n)      = trans(draw);
  } // END n loop
  
  return aux_A;
} // END sample_A_heterosk1_mss




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec sample_hyperparameters_s4 (
    arma::vec               aux_hyper,
    const arma::mat&        aux_B,
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior
) {
  // the function draws the value of aux_hyper
  const int N = aux_B.n_rows;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu     = as<double>(prior["hyper_nu"]);
  double prior_hyper_s      = as<double>(prior["hyper_s"]);
  
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
  
  aux_hyper(1)    = ( prior_hyper_s + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
    R::rchisq( prior_hyper_nu + N * K );
  aux_hyper(0)    = ( prior_hyper_s + trace(aux_B * as<mat>(prior["B_V_inv"]) * trans(aux_B) )) /
    R::rchisq( prior_hyper_nu + rn );
  
  return aux_hyper;
} // END sample_hyperparameters_s4



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec sample_hyperparameters_mss (
    arma::vec               aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const Rcpp::List&       prior
) {
  // the function changes the value of aux_hyper by reference (filling it with a new draw)
  const int N = aux_B.n_rows;
  const int M = aux_B.n_slices;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu     = as<double>(prior["hyper_nu"]);
  double prior_hyper_s      = as<double>(prior["hyper_s"]);
  
  int rn=0;
  for (int n=0; n<N; n++) {
    rn       += VB(n).n_rows;
  }
  
  aux_hyper(1)    = ( prior_hyper_s + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
    R::rchisq( prior_hyper_nu + N * K );
  double BVB      = 0;
  for (int m=0; m<M; m++) {
    BVB          += trace(aux_B.slice(m) * as<mat>(prior["B_V_inv"]) * trans(aux_B.slice(m)) );
  }
  aux_hyper(0)    = ( prior_hyper_s + BVB) /
    R::rchisq( prior_hyper_nu + M * rn );
  
  return aux_hyper;
} // END sample_hyperparameters_mss



/*______________________function sample_hyperparameters_mss_s4______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyperparameters_mss_s4 (
    arma::vec               aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior
) {
  // the function changes the value of aux_hyper by reference
  const int N = aux_B.n_rows;
  const int M = aux_B.n_slices;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu     = as<double>(prior["hyper_nu"]);
  double prior_hyper_s      = as<double>(prior["hyper_s"]);
  
  int         Ltmp  = VB.n_elem - 1;
  vec         Lm    = VB(Ltmp);
  vec         Lm_cs = cumsum(Lm);
  int rn=0;
  for (int m=0; m<M; m++) {
    for (int n=0; n<N; n++) {
      int ll        = aux_SL(n,m);
      if (n>0) {
        ll         += Lm_cs(n-1);
      }
      rn           += VB(ll).n_rows;
    }
  }
  
  aux_hyper(1)    = ( prior_hyper_s + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
    R::rchisq( prior_hyper_nu + N * K );
  double BVB      = 0;
  for (int m=0; m<M; m++) {
    BVB          += trace(aux_B.slice(m) * as<mat>(prior["B_V_inv"]) * trans(aux_B.slice(m)) );
  }
  aux_hyper(0)    = ( prior_hyper_s + BVB) /
    R::rchisq( prior_hyper_nu + rn );
  
  return aux_hyper;
} // END sample_hyperparameters_mss_s4


