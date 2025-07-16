
#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"

#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;


// 
// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::mat sample_B_heterosk1 (
//     arma::mat         aux_B,          // NxN
//     const arma::mat&  aux_A,          // NxK
//     const arma::mat&  aux_hyper,      //
//     const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
//     const arma::mat&  Y,              // NxT dependent variables
//     const arma::mat&  X,              // KxT dependent variables
//     const Rcpp::List& prior,          // a list of priors - original dimensions
//     const arma::field<arma::mat>& VB        // restrictions on B0
// ) {
//   // the function changes the value of aux_B0 and aux_Bplus by reference (filling it with a new draw)
//   const int N               = aux_B.n_rows;
//   const int T               = Y.n_cols;
//   
//   const int posterior_nu    = T + as<int>(prior["B_nu"]);
//   mat prior_SS_inv          = pow(aux_hyper(0), -1) * as<mat>(prior["B_V_inv"]);
//   mat shocks                = Y - aux_A * X;
//   
//   
//   for (int n=0; n<N; n++) {
//     
//     // set scale matrix
//     mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
//     mat posterior_SS_inv    = prior_SS_inv + shocks_sigma * shocks_sigma.t();
//     mat posterior_S_inv     = VB(n) * posterior_SS_inv * VB(n).t();
//     posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
//     
//     // sample B
//     mat Un                  = chol(posterior_nu * inv_sympd(posterior_S_inv));
//     mat B_tmp               = aux_B;
//     B_tmp.shed_row(n);
//     rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
//     vec w1_tmp              = trans(w * VB(n).t() * Un.t());
//     double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
//     mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
//     mat Wn;
//     const int rn            = VB(n).n_rows;
//     if (rn==1) {
//       Wn                    = w1;
//     } else {
//       Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
//     }
//     
//     vec   alpha(rn);
//     vec   u(posterior_nu+1, fill::randn);
//     u                      *= pow(posterior_nu, -0.5);
//     alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
//     if (R::runif(0,1)<0.5) {
//       alpha(0)       *= -1;
//     }
//     if (rn>1){
//       vec nn(rn-1, fill::randn);
//       nn                   *= pow(posterior_nu, -0.5);
//       alpha.rows(1,rn-1)    = nn;
//     }
//     rowvec b0n              = alpha.t() * Wn * Un;
//     aux_B.row(n)           = b0n * VB(n);
//   } // END n loop
//   
//   return aux_B;
// } // END sample_B_heterosk1



// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::cube sample_B_mss (
//     arma::cube        aux_B,          // NxNxM
//     const arma::mat&  aux_A,          // NxK
//     const arma::vec&  aux_hyper,      // NxM
//     const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
//     const arma::mat&  aux_xi,         // MxT
//     const arma::mat&  Y,              // NxT dependent variables
//     const arma::mat&  X,              // KxT dependent variables
//     const Rcpp::List& prior,          // a list of priors - original dimensions
//     const arma::field<arma::mat>& VB        // restrictions on B0
// ) {
//   const int T     = Y.n_cols;
//   const int N     = Y.n_rows;
//   const int K     = X.n_rows;
//   const int M     = aux_B.n_slices;
//   const vec T_m  = sum(aux_xi,1);
//   
//   // aux_sigma, Y, X
//   
//   for (int m=0; m<M; m++) {
//     int ii = 0;
//     mat aux_sigma_m(N, T_m(m));
//     mat Y_m(N, T_m(m));
//     mat X_m(K, T_m(m));
//     
//     for (int t=0; t<T; t++){
//       if (aux_xi(m,t)==1) {
//         aux_sigma_m.col(ii) = aux_sigma.col(t);
//         Y_m.col(ii)         = Y.col(t);
//         X_m.col(ii)         = X.col(t);
//         ii++;
//       }
//     }
//     aux_B.slice(m)    = sample_B_heterosk1(aux_B.slice(m), aux_A, aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
//   }
//   
//   return aux_B;
// } // END sample_B_mss








// /*______________________function sample_B_mss_s4______________________*/
// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// Rcpp::List sample_B_mss_s4 (
//     arma::cube        aux_B,          // NxNxM
//     arma::imat        aux_SL,         // NxM row-specific S4 indicators
//     const arma::mat&  aux_A,          // NxK
//     const arma::vec&  aux_hyper,      // NxM
//     const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
//     const arma::mat&  aux_xi,         // MxT
//     const arma::mat&  Y,              // NxT dependent variables
//     const arma::mat&  X,              // KxT dependent variables
//     const Rcpp::List& prior,          // a list of priors - original dimensions
//     const arma::field<arma::mat>& VB        // restrictions on B0
// ) {
//   const int T     = Y.n_cols;
//   const int N     = Y.n_rows;
//   const int K     = X.n_rows;
//   const int M     = aux_B.n_slices;
//   const vec T_m  = sum(aux_xi,1);
//   
//   // aux_sigma, Y, X
//   
//   for (int m=0; m<M; m++) {
//     int ii = 0;
//     mat aux_sigma_m(N, T_m(m));
//     mat Y_m(N, T_m(m));
//     mat X_m(K, T_m(m));
//     
//     for (int t=0; t<T; t++){
//       if (aux_xi(m,t)==1) {
//         aux_sigma_m.col(ii) = aux_sigma.col(t);
//         Y_m.col(ii)         = Y.col(t);
//         X_m.col(ii)         = X.col(t);
//         ii++;
//       }
//     }
//     // ivec aux_SL_m           = aux_SL.col(m);
//     List BSL_m              = sample_B_heterosk1_s4(aux_B.slice(m), aux_SL.col(m), aux_A, aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
//     aux_B.slice(m)          = as<mat>(BSL_m["aux_B"]);
//     aux_SL.col(m)           = as<ivec>(BSL_m["aux_SL"]);
//   }
//   
//   return List::create(
//       _["aux_B"]    = aux_B,
//       _["aux_SL"]   = aux_SL
//   );
// } // END sample_B_mss_s4



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_B_heterosk1_boost (
    arma::mat         aux_B,          // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_hyper,      // (2*N+1)x2
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
  mat prior_SS_inv          = as<mat>(prior["B_V_inv"]);
  mat shocks                = Y - aux_A * X;
  
  for (int n=0; n<N; n++) {
    
    Rcout << " structural equation: " << n + 1 << endl;
    // set scale matrix
    mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
    mat posterior_SS_inv    = pow(aux_hyper(n, 0), -1) * prior_SS_inv + shocks_sigma * shocks_sigma.t();
    mat posterior_S_inv     = VB(n) * posterior_SS_inv * VB(n).t();
    posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
    
    // sample B
    mat posterior_S(size(posterior_S_inv));
    mat Un(size(posterior_S_inv));
    
    try {
      posterior_S           = inv_sympd(posterior_S_inv);
    } catch (std::runtime_error &e) {
      Rcout << "   B Inversion failure " << endl;
      continue; 
    }
    
    try {
      Un                    = chol(posterior_nu * posterior_S);
    } catch (std::runtime_error &e) {
      Rcout << "   B Cholesky failure " << endl;
      continue; 
    }
    
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
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
} // END sample_B_heterosk1_boost



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube sample_B_mss_boost (
    arma::cube        aux_B,          // NxNxM
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_hyper,      // (2*N+1)x2
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
    
    aux_B.slice(m)    = sample_B_heterosk1_boost(aux_B.slice(m), aux_A, aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
  }
  
  return aux_B;
} // END sample_B_mss_boost








// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_B_mss_s4_boost (
    arma::cube        aux_B,          // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_hyper,      // (2*N+1x2)
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
    
    List BSL_m              = sample_B_heterosk1_s4_boost(aux_B.slice(m), aux_SL.col(m), aux_A, aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
    aux_B.slice(m)          = as<mat>(BSL_m["aux_B"]);
    aux_SL.col(m)           = as<ivec>(BSL_m["aux_SL"]);
  }
  
  return List::create(
    _["aux_B"]    = aux_B,
    _["aux_SL"]   = aux_SL
  );
} // END sample_B_mss_s4_boost





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_B_mssa_s4_boost (
    arma::cube        aux_B,          // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    const arma::cube& aux_A,          // NxKxM
    const arma::mat&  aux_hyper,      // (2*N+1x2)
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
    
    List BSL_m              = sample_B_heterosk1_s4_boost(aux_B.slice(m), aux_SL.col(m), aux_A.slice(m), aux_hyper, aux_sigma_m, Y_m, X_m, prior, VB);
    aux_B.slice(m)          = as<mat>(BSL_m["aux_B"]);
    aux_SL.col(m)           = as<ivec>(BSL_m["aux_SL"]);
  }
  
  return List::create(
    _["aux_B"]    = aux_B,
    _["aux_SL"]   = aux_SL
  );
} // END sample_B_mssa_s4_boost


















// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::mat sample_A_heterosk1 (
//     arma::mat         aux_A,          // NxK
//     const arma::mat&  aux_B,          // NxN
//     const arma::vec&  aux_hyper,      // NxM
//     const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
//     const arma::mat&  Y,              // NxT dependent variables
//     const arma::mat&  X,              // KxT dependent variables
//     const Rcpp::List& prior           // a list of priors - original dimensions
// ) {
//   // the function changes the value of aux_A by reference
//   const int N         = aux_A.n_rows;
//   const int K         = aux_A.n_cols;
//   
//   mat prior_A_mean    = as<mat>(prior["A"]);
//   mat prior_A_Vinv    = pow(aux_hyper(1), -1) * as<mat>(prior["A_V_inv"]);
//   rowvec    zerosA(K);
//   vec sigma_vectorised= vectorise(aux_sigma);
//   
//   for (int n=0; n<N; n++) {
//     mat   A0          = aux_A;
//     A0.row(n)         = zerosA;
//     vec   zn          = vectorise( aux_B * (Y - A0 * X) );
//     mat   zn_sigma    = zn / sigma_vectorised;
//     mat   Wn          = kron( trans(X), aux_B.col(n) );
//     mat   Wn_sigma    = Wn.each_col() / sigma_vectorised;
//     
//     mat     precision = prior_A_Vinv + trans(Wn_sigma) * Wn_sigma;
//     precision         = 0.5 * (precision + precision.t());
//     rowvec  location  = prior_A_mean.row(n) * prior_A_Vinv + trans(zn_sigma) * Wn_sigma;
//     
//     mat     precision_chol = trimatu(chol(precision));
//     vec     xx(K, fill::randn);
//     vec     draw      = solve(precision_chol, 
//                               solve(trans(precision_chol), trans(location)) + xx);
//     aux_A.row(n)      = trans(draw);
//   } // END n loop
//   
//   return aux_A;
// } // END sample_A_heterosk1





// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::mat sample_A_heterosk1_mss (
//     arma::mat         aux_A,          // NxK
//     const arma::cube& aux_B,          // NxNxM
//     const arma::mat&  aux_xi,         // MxT
//     const arma::vec&  aux_hyper,      // NxM
//     const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
//     const arma::mat&  Y,              // NxT dependent variables
//     const arma::mat&  X,              // KxT dependent variables
//     const Rcpp::List& prior           // a list of priors - original dimensions
// ) {
//   // the function draws the value of aux_A
//   const int N         = aux_A.n_rows;
//   const int K         = aux_A.n_cols;
//   const int T         = Y.n_cols;
//   const int M         = aux_xi.n_rows;
//   const vec Tm       = sum(aux_xi, 1);
//   
//   mat prior_A_mean    = as<mat>(prior["A"]);
//   mat prior_A_Vinv    = pow(aux_hyper(1), -1) * as<mat>(prior["A_V_inv"]);
//   rowvec    zerosA(K);
//   
//   field<mat> YM(M);
//   field<mat> XM(M);
//   field<mat> SM(M);
//   for (uword m=0; m<M; m++) {
//     mat Y_m(N, Tm(m));
//     mat X_m(K, Tm(m));
//     mat S_m(N, Tm(m));
//     
//     YM(m)             = Y_m;
//     XM(m)             = X_m;
//     SM(m)             = S_m;
//   }
//   vec m_iter(M);
//   
//   for (uword t=0; t<T; t++) {
//     int mm            = aux_xi.col(t).index_max();
//     YM(mm).col(m_iter(mm))    = Y.col(t);
//     XM(mm).col(m_iter(mm))    = X.col(t);
//     SM(mm).col(m_iter(mm))    = aux_sigma.col(t);
//     m_iter += aux_xi.col(t);
//   }
//   
//   
//   for (int n=0; n<N; n++) {
//     mat   A0          = aux_A;
//     A0.row(n)         = zerosA;
//     vec   zn          = vectorise( aux_B.slice(0) * (YM(0) - A0 * XM(0)) );
//     mat   Wn          = kron( trans(XM(0)), aux_B.slice(0).col(n) );
//     vec   Sn          = vectorise(SM(0));
//     
//     for (uword m=1; m<M; m++) {
//       zn          = join_cols(zn, vectorise(aux_B.slice(m) * (YM(m) - A0 * XM(m))));
//       Wn          = join_cols(Wn, kron(trans(XM(m)), aux_B.slice(m).col(n)));
//       Sn          = join_cols(Sn, vectorise(SM(m)));
//     }
//     
//     mat   zn_sigma    = zn / Sn;
//     mat   Wn_sigma    = Wn.each_col() / Sn;
//     
//     mat     precision = prior_A_Vinv + trans(Wn_sigma) * Wn_sigma;
//     precision         = 0.5 * (precision + precision.t());
//     rowvec  location  = prior_A_mean.row(n) * prior_A_Vinv + trans(zn_sigma) * Wn_sigma;
//     
//     mat     precision_chol = trimatu(chol(precision));
//     vec     draw      = solve(precision_chol, 
//                               solve(trans(precision_chol), trans(location)) + as<vec>(rnorm(K)));
//     aux_A.row(n)      = trans(draw);
//   } // END n loop
//   
//   return aux_A;
// } // END sample_A_heterosk1_mss



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_heterosk1_boost (
    arma::mat         aux_A,          // NxK
    const arma::mat&  aux_B,          // NxN
    const arma::mat&  aux_hyper,      // (2*N+1)x2
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior           // a list of priors - original dimensions
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_V_inv   = as<mat>(prior["A_V_inv"]);
  
  rowvec    zerosA(K);
  vec sigma_vectorised= vectorise(aux_sigma);
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B * (Y - A0 * X) );
    mat   zn_sigma    = zn / sigma_vectorised;
    mat   Wn          = kron( trans(X), aux_B.col(n) );
    mat   Wn_sigma    = Wn.each_col() / sigma_vectorised;
    
    mat     precision = pow(aux_hyper(n, 1), -1) * prior_A_V_inv + trans(Wn_sigma) * Wn_sigma;
    precision         = 0.5 * (precision + precision.t());
    rowvec  location  = prior_A_mean.row(n) * pow(aux_hyper(n, 1), -1) * prior_A_V_inv + trans(zn_sigma) * Wn_sigma;
    
    mat     precision_chol = trimatu(chol(precision));
    vec     xx(K, fill::randn);
    vec     draw      = solve(precision_chol, 
                              solve(trans(precision_chol), trans(location)) + xx);
    aux_A.row(n)      = trans(draw);
  } // END n loop
  
  return aux_A;
} // END sample_A_heterosk1_boost



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_heterosk1_mss_boost (
    arma::mat         aux_A,          // NxK
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  aux_hyper,      // (2*N+1)x2
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
  mat prior_A_V_inv    = as<mat>(prior["A_V_inv"]);
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
    Rcout << " AR equation: " << n + 1 << endl;
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
    
    // Rcout << " pow(aux_hyper(n, 1), -1): " << pow(aux_hyper(n, 1), -1) << endl;
    mat     precision = pow(aux_hyper(n, 1), -1) * prior_A_V_inv + trans(Wn_sigma) * Wn_sigma;
    precision         = 0.5 * (precision + precision.t());
    rowvec  location  = prior_A_mean.row(n) * pow(aux_hyper(n, 1), -1) * prior_A_V_inv + trans(zn_sigma) * Wn_sigma;
    
    mat     precision_chol(size(precision));
    try {
      precision_chol = trimatu(chol(precision));
    } catch (std::runtime_error &e) {
      Rcout << "   Cholesky failure " << endl;
      continue; 
    }
    vec     draw      = solve(precision_chol, 
                              solve(trans(precision_chol), trans(location)) + as<vec>(rnorm(K)));
    aux_A.row(n)      = trans(draw);
  } // END n loop
  
  return aux_A;
} // END sample_A_heterosk1_mss_boost





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube sample_A_heterosk1_mssa_boost (
    arma::cube        aux_A,          // NxKxM
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    const arma::mat&  aux_hyper,      // (2*N+1)x2
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior           // a list of priors - original dimensions
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  const int T         = Y.n_cols;
  const int M         = aux_xi.n_rows;
  const vec Tm        = sum(aux_xi, 1);
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_V_inv   = as<mat>(prior["A_V_inv"]);
  
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
  
  for (int m=0; m<M; m++) {

    vec sigma_vectorised= vectorise(SM(m));
    
    for (int n=0; n<N; n++) {
      mat   A0          = aux_A.slice(m);
      A0.row(n)         = zerosA;
      vec   zn          = vectorise( aux_B.slice(m) * (YM(m) - A0 * XM(m)) );
      mat   zn_sigma    = zn / sigma_vectorised;
      mat   Wn          = kron( trans(XM(m)), aux_B.slice(m).col(n) );
      mat   Wn_sigma    = Wn.each_col() / sigma_vectorised;
      
      mat     precision = pow(aux_hyper(n, 1), -1) * prior_A_V_inv + trans(Wn_sigma) * Wn_sigma;
      precision         = 0.5 * (precision + precision.t());
      rowvec  location  = prior_A_mean.row(n) * pow(aux_hyper(n, 1), -1) * prior_A_V_inv + trans(zn_sigma) * Wn_sigma;
      
      mat     precision_chol = trimatu(chol(precision));
      vec     xx(K, fill::randn);
      vec     draw          = solve(precision_chol, 
                                solve(trans(precision_chol), trans(location)) + xx);
      aux_A.slice(m).row(n) = trans(draw);
    } // END n loop
  } // END m loop
  
  return aux_A;
} // END sample_A_heterosk1_mssa_boost









// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::vec sample_hyperparameters_s4 (
//     arma::vec               aux_hyper,
//     const arma::mat&        aux_B,
//     const arma::mat&        aux_A,
//     const arma::field<arma::mat>& VB,
//     const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
//     const Rcpp::List&       prior
// ) {
//   // the function draws the value of aux_hyper
//   const int N = aux_B.n_rows;
//   const int K = aux_A.n_cols;
//   
//   double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
//   double prior_hyper_s_B      = as<double>(prior["hyper_s_B"]);
//   double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
//   double prior_hyper_s_A      = as<double>(prior["hyper_s_A"]);
//   // double prior_hyper_a      = as<double>(prior["hyper_a"]);
//   
//   int         Ltmp = VB.n_elem - 1;
//   vec         Lm  = VB(Ltmp);
//   vec         Lm_cs = cumsum(Lm);
//   int rn=0;
//   for (int n=0; n<N; n++) {
//     int ll        = aux_SL(n);
//     if (n>0) {
//       ll         += Lm_cs(n-1);
//     }
//     rn           += VB(ll).n_rows;
//   }
//   
//   // aux_hyper(4)    = ( as<double>(prior["hyper_S"]) + aux_hyper(2) + aux_hyper(3) ) / R::rchisq(as<double>(prior["hyper_V"]) + 4 * prior_hyper_a);
//   // aux_hyper(3)    = R::rgamma( prior_hyper_a + 0.5 * prior_hyper_nu ,
//   //           1/((1/aux_hyper(4)) + (1/(2*aux_hyper(1)))));
//   // aux_hyper(2)    = R::rgamma( prior_hyper_a + 0.5 * prior_hyper_nu ,
//   //           1/((1/100) + (1/(2*aux_hyper(0)))) );
//   aux_hyper(1)    = ( prior_hyper_s_A + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
//     R::rchisq( prior_hyper_nu_A + N * K );
//   aux_hyper(0)    = ( prior_hyper_s_B + trace(aux_B * as<mat>(prior["B_V_inv"]) * trans(aux_B) )) /
//     R::rchisq( prior_hyper_nu_B + rn );
//   
//   return aux_hyper;
// } // END sample_hyperparameters_s4



// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::vec sample_hyperparameters_mss (
//     arma::vec               aux_hyper,
//     const arma::cube&       aux_B,            // NxNxM
//     const arma::mat&        aux_A,
//     const arma::field<arma::mat>& VB,
//     const Rcpp::List&       prior
// ) {
//   // the function changes the value of aux_hyper by reference (filling it with a new draw)
//   const int N = aux_B.n_rows;
//   const int M = aux_B.n_slices;
//   const int K = aux_A.n_cols;
//   
//   double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
//   double prior_hyper_s_B      = as<double>(prior["hyper_s_B"]);
//   double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
//   double prior_hyper_s_A      = as<double>(prior["hyper_s_A"]);
//   
//   int rn=0;
//   for (int n=0; n<N; n++) {
//     rn       += VB(n).n_rows;
//   }
//   
//   // aux_hyper(4)    = ( as<double>(prior["hyper_S"]) + aux_hyper(2) + aux_hyper(3) ) / R::rchisq(as<double>(prior["hyper_V"]) + 4 * prior_hyper_a);
//   // aux_hyper(3)    = R::rgamma( prior_hyper_a + 0.5 * prior_hyper_nu ,
//   //           1/((1/aux_hyper(4)) + (1/(2*aux_hyper(1)))));
//   // aux_hyper(2)    = R::rgamma( prior_hyper_a + 0.5 * prior_hyper_nu ,
//   //           1/((1/100) + (1/(2*aux_hyper(0)))) );
//   aux_hyper(1)    = ( prior_hyper_s_A + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
//     R::rchisq( prior_hyper_nu_A + N * K );
//   double BVB      = 0;
//   for (int m=0; m<M; m++) {
//     BVB          += trace(aux_B.slice(m) * as<mat>(prior["B_V_inv"]) * trans(aux_B.slice(m)) );
//   }
//   aux_hyper(0)    = ( prior_hyper_s_B + BVB) /
//     R::rchisq( prior_hyper_nu_B + M * rn );
//   
//   return aux_hyper;
// } // END sample_hyperparameters_mss



// /*______________________function sample_hyperparameters_mss_s4______________________*/
// // [[Rcpp::interfaces(cpp)]]
// // [[Rcpp::export]]
// arma::mat sample_hyperparameters_mss_s4 (
//     arma::vec               aux_hyper,
//     const arma::cube&       aux_B,            // NxNxM
//     const arma::mat&        aux_A,
//     const arma::field<arma::mat>& VB,
//     const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
//     const Rcpp::List&       prior
// ) {
//   // the function changes the value of aux_hyper by reference
//   const int N = aux_B.n_rows;
//   const int M = aux_B.n_slices;
//   const int K = aux_A.n_cols;
//   
//   double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
//   double prior_hyper_s_B      = as<double>(prior["hyper_s_B"]);
//   double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
//   double prior_hyper_s_A      = as<double>(prior["hyper_s_A"]);
//   
//   int         Ltmp  = VB.n_elem - 1;
//   vec         Lm    = VB(Ltmp);
//   vec         Lm_cs = cumsum(Lm);
//   int rn=0;
//   for (int m=0; m<M; m++) {
//     for (int n=0; n<N; n++) {
//       int ll        = aux_SL(n,m);
//       if (n>0) {
//         ll         += Lm_cs(n-1);
//       }
//       rn           += VB(ll).n_rows;
//     }
//   }
//   
//   // aux_hyper(4)    = ( as<double>(prior["hyper_S"]) + aux_hyper(2) + aux_hyper(3) ) / R::rchisq(as<double>(prior["hyper_V"]) + 4 * prior_hyper_a);
//   // aux_hyper(3)    = R::rgamma( prior_hyper_a + 0.5 * prior_hyper_nu ,
//   //           1/((1/aux_hyper(4)) + (1/(2*aux_hyper(1)))));
//   // aux_hyper(2)    = R::rgamma( prior_hyper_a + 0.5 * prior_hyper_nu ,
//   //           1/((1/100) + (1/(2*aux_hyper(0)))) );
//   aux_hyper(1)    = ( prior_hyper_s_A + trace((aux_A - as<mat>(prior["A"])) * as<mat>(prior["A_V_inv"]) * trans(aux_A - as<mat>(prior["A"]))) ) /
//     R::rchisq( prior_hyper_nu_A + N * K );
//   double BVB      = 0;
//   for (int m=0; m<M; m++) {
//     BVB          += trace(aux_B.slice(m) * as<mat>(prior["B_V_inv"]) * trans(aux_B.slice(m)) );
//   }
//   aux_hyper(0)    = ( prior_hyper_s_B + BVB) /
//     R::rchisq( prior_hyper_nu_B + rn );
//   
//   return aux_hyper;
// } // END sample_hyperparameters_mss_s4





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyperparameter_boost_s4 (
    arma::mat               aux_hyper,      // (2 * N + 1) x 2
    const arma::mat&        aux_B,
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::ivec&       aux_SL,         // Nx1 row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
) {
  // the function draws the value of aux_hyper
  const int N = aux_B.n_rows;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  mat   prior_A               = as<mat>(prior["A"]);
  mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  int         Ltmp = VB.n_elem - 1;
  vec         Lm  = VB(Ltmp);
  vec         Lm_cs = cumsum(Lm);
  ivec        rn(N);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  BB_scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  BB_shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 0) = BB_scale_tmp / R::rchisq(BB_shape_tmp);
  }
  
  // aux_A - related hyper-parameters 
  ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  double AA_scale_tmp = prior_hyper_s_AA + 2 * sum(ss_tmp);
  double AA_shape_tmp           = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 1) = AA_scale_tmp / R::rchisq(AA_shape_tmp);
  }
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int ll        = aux_SL(n);
    if (n>0) {
      ll         += Lm_cs(n-1);
    }
    rn(n)         = VB(ll).n_rows;
    
    // aux_B - related hyper-parameters 
    if ( hyper_boost ) {
      BB_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
      BB_shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
      aux_hyper(N + n, 0) = R::rgamma(BB_shape_tmp, BB_scale_tmp);
      
      BB_scale_tmp         = aux_hyper(N + n, 0) + 
        as_scalar(aux_B.row(n) * prior_B_V_inv * trans(aux_B.row(n)));
      BB_shape_tmp         = prior_hyper_nu_B + rn(n);
      aux_hyper(n, 0)   = BB_scale_tmp / R::rchisq(BB_shape_tmp);
    }
    
    // aux_A - related hyper-parameters 
    if ( hyper_boost ) {
      AA_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
      AA_shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
      aux_hyper(N + n, 1) = R::rgamma(AA_shape_tmp, AA_scale_tmp);
      
      AA_scale_tmp         = aux_hyper(N + n, 1) + 
        as_scalar((aux_A.row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.row(n) - prior_A.row(n)));
      AA_shape_tmp         = prior_hyper_nu_A + K;
      aux_hyper(n, 1)   = AA_scale_tmp / R::rchisq(AA_shape_tmp);
    }
  }
  
  return aux_hyper;
} // END sample_hyperparameter_boost_s4



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyperparameters_mss_boost (
    arma::mat               aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
) {
  // the function returns aux_hyper
  
  const int N = aux_B.n_rows;
  const int M = aux_B.n_slices;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  mat   prior_A               = as<mat>(prior["A"]);
  mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  BB_scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  BB_shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 0) = BB_scale_tmp / R::rchisq(BB_shape_tmp);
  }
  
  // aux_A - related hyper-parameters 
  ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  double AA_scale_tmp = prior_hyper_s_AA + 2 * sum(ss_tmp);
  double AA_shape_tmp = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 1) = AA_scale_tmp / R::rchisq(AA_shape_tmp);
  }
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int rn            = VB(n).n_rows;
    
    // aux_B - related hyper-parameters 
    if ( hyper_boost ) {
      BB_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
      BB_shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
      aux_hyper(N + n, 0) = R::rgamma(BB_shape_tmp, BB_scale_tmp);
      
      double BVB        = 0;
      for (int m=0; m<M; m++) {
        rowvec  Bnm     = aux_B.subcube(n, 0, m, n, N - 1, m);
        BVB            += as_scalar(Bnm * prior_B_V_inv * Bnm.t());
      }
      BB_scale_tmp         = aux_hyper(N + n, 0) + BVB;
      BB_shape_tmp         = prior_hyper_nu_B + M * rn;
      aux_hyper(n, 0)   = BB_scale_tmp / R::rchisq(BB_shape_tmp);
    }
    
    // aux_A - related hyper-parameters 
    if ( hyper_boost ) {
      AA_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
      AA_shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
      aux_hyper(N + n, 1) = R::rgamma(AA_shape_tmp, AA_scale_tmp);
      
      AA_scale_tmp         = aux_hyper(N + n, 1) +
        as_scalar((aux_A.row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.row(n) - prior_A.row(n)));
      AA_shape_tmp         = prior_hyper_nu_A + K;
      aux_hyper(n, 1)   = AA_scale_tmp / R::rchisq(AA_shape_tmp);
    }
  }
  
  return aux_hyper;
} // END sample_hyperparameters_mss_boost



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyperparameters_mss_s4_boost (
    arma::mat               aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
) {
  // the function changes the value of aux_hyper by reference
  const int N = aux_B.n_rows;
  const int M = aux_B.n_slices;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  mat   prior_A               = as<mat>(prior["A"]);
  mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  int   Ltmp                  = VB.n_elem - 1;
  vec   Lm                    = VB(Ltmp);
  vec   Lm_cs                 = cumsum(Lm);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  BB_scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  BB_shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 0) = BB_scale_tmp / R::rchisq(BB_shape_tmp);
  }
  
  // aux_A - related hyper-parameters 
  ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  double AA_scale_tmp           = prior_hyper_s_AA + 2 * sum(ss_tmp);
  double AA_shape_tmp           = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 1) = AA_scale_tmp / R::rchisq(AA_shape_tmp);
  }
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int rn=0;
    for (int m=0; m<M; m++) {
      int ll        = aux_SL(n,m);
      if (n>0) {
        ll         += Lm_cs(n-1);
      }
      rn           += VB(ll).n_rows;
    }
    
    // aux_B - related hyper-parameters 
    if ( hyper_boost ) {
      BB_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
      BB_shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
      aux_hyper(N + n, 0) = R::rgamma(BB_shape_tmp, BB_scale_tmp);
      
      double BVB        = 0;
      for (int m=0; m<M; m++) {
        rowvec  Bnm     = aux_B.subcube(n, 0, m, n, N - 1, m);
        BVB            += as_scalar(Bnm * prior_B_V_inv * Bnm.t());
      }
      BB_scale_tmp         = aux_hyper(N + n, 0) + BVB;
      BB_shape_tmp         = prior_hyper_nu_B + rn;
      aux_hyper(n, 0)   = BB_scale_tmp / R::rchisq(BB_shape_tmp);
    }
    
    // aux_A - related hyper-parameters 
    if ( hyper_boost ) {
      AA_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
      AA_shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
      aux_hyper(N + n, 1) = R::rgamma(AA_shape_tmp, AA_scale_tmp);
      
      AA_scale_tmp         = aux_hyper(N + n, 1) + 
        as_scalar((aux_A.row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.row(n) - prior_A.row(n)));
      AA_shape_tmp         = prior_hyper_nu_A + K;
      aux_hyper(n, 1)   = AA_scale_tmp / R::rchisq(AA_shape_tmp);
    }
  }
  
  return aux_hyper;
} // END sample_hyperparameters_mss_s4_boost



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyperparameters_mssa_s4_boost (
    arma::mat               aux_hyper,
    const arma::cube&       aux_B,            // NxNxM
    const arma::cube&       aux_A,            // NxKxM
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
) {
  // the function changes the value of aux_hyper by reference
  const int N = aux_B.n_rows;
  const int M = aux_B.n_slices;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  mat   prior_A               = as<mat>(prior["A"]);
  mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  int   Ltmp                  = VB.n_elem - 1;
  vec   Lm                    = VB(Ltmp);
  vec   Lm_cs                 = cumsum(Lm);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  BB_scale_tmp  = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  BB_shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 0) = BB_scale_tmp / R::rchisq(BB_shape_tmp);
  }
  
  // aux_A - related hyper-parameters 
  ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  double AA_scale_tmp = prior_hyper_s_AA + 2 * sum(ss_tmp);
  double AA_shape_tmp = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  if ( hyper_boost ) {
    aux_hyper(2 * N, 1) = AA_scale_tmp / R::rchisq(AA_shape_tmp);
  }
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int rn=0;
    for (int m=0; m<M; m++) {
      int ll        = aux_SL(n,m);
      if (n>0) {
        ll         += Lm_cs(n-1);
      }
      rn           += VB(ll).n_rows;
    }
    
    // aux_B - related hyper-parameters 
    if ( hyper_boost ) {
      BB_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
      BB_shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
      aux_hyper(N + n, 0) = R::rgamma(BB_shape_tmp, BB_scale_tmp);
      
      double BVB        = 0;
      for (int m=0; m<M; m++) {
        rowvec  Bnm     = aux_B.subcube(n, 0, m, n, N - 1, m);
        BVB            += as_scalar(Bnm * prior_B_V_inv * Bnm.t());
      }
      BB_scale_tmp         = aux_hyper(N + n, 0) + BVB;
      BB_shape_tmp         = prior_hyper_nu_B + rn;
      aux_hyper(n, 0)   = BB_scale_tmp / R::rchisq(BB_shape_tmp);
    }
    
    // aux_A - related hyper-parameters 
    if ( hyper_boost ) {
      AA_scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
      AA_shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
      aux_hyper(N + n, 1) = R::rgamma(AA_shape_tmp, AA_scale_tmp);
      
      AA_scale_tmp         = aux_hyper(N + n, 1);
      for (int m=0; m<M; m++) {
        AA_scale_tmp      += as_scalar((aux_A.slice(m).row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.slice(m).row(n) - prior_A.row(n)));
      }
      AA_shape_tmp         = prior_hyper_nu_A + N * K;
      aux_hyper(n, 1)   = AA_scale_tmp / R::rchisq(AA_shape_tmp);
    }
  }
  
  return aux_hyper;
} // END sample_hyperparameters_mssa_s4_boost
