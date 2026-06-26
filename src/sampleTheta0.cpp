
#include <RcppArmadillo.h>
#include "sample_sv_ms.h"
#include "sample_ABhyper.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec Rodriguez_Yam_2004 (
    arma::vec& mu,
    arma::mat& Sig,
    arma::mat& A,
    arma::vec& l,
    arma::vec& u,
    arma::vec& x
) {
  
  int n = mu.n_elem;
  mat C = trimatu( chol(Sig, "lower") );
  mat Atilde = A * C;
  vec ltilde = l - A * mu;
  vec utilde = u - A * mu;
  vec z = solve(C, x - mu);
  uvec idx(n);
  mat lutilde_s;
  vec ltildei, utildei;
  
  for (int i = 0; i < n; ++i) {
    
    idx = regspace<uvec>(0, n - 1);
    idx.shed_row(i);
    
    ltildei = (ltilde - Atilde.cols(idx) * z.elem(idx)) / Atilde.col(i);
    utildei = (utilde - Atilde.cols(idx) * z.elem(idx)) / Atilde.col(i);
    
    mat lutilde(ltildei.n_elem, 2);
    lutilde.col(0) = ltildei;
    lutilde.col(1) = utildei;
    
    lutilde_s = arma::sort(lutilde, "ascend", 1);
    
    double maxli = max(lutilde_s.col(0));
    double minui = min(lutilde_s.col(1));
    
    double u_rand = randu();
    double zi = R::qnorm(
      R::pnorm(maxli, 0.0, 1.0, 1, 0)
      + u_rand * (R::pnorm(minui, 0.0, 1.0, 1, 0) - R::pnorm(maxli, 0.0, 1.0, 1, 0)),
      0.0, 1.0, 1, 0
    );
    
    z(i) = zi;
    
  }
  
  vec out = mu + C * z;
  return out;
} // END Rodriguez_Yam_2004



// [[Rcpp::export]]
double AN (
    const double& mu,
    const double& rho
) {
  // a sampler for the absolute normal distribution by Villani (2004, JAE)
  
  double mu1 = 0.5 * mu - 0.5 * sqrt(mu * mu + 4.0);
  double mu2 = 0.5 * mu + 0.5 * sqrt(mu * mu + 4.0);
  double sig1 = mu1 * mu1 / (1.0 + mu1 * mu1) * rho;
  double sig2 = mu2 * mu2 / (1.0 + mu2 * mu2) * rho;
  double w = 1.0 / (1.0 + exp(2.0 * mu / rho));
  double x = randn();
  if ( randu() < w ) {
    x *= sqrt(sig1);
    x += mu1;
  } else {
    x *= sqrt(sig2);
    x += mu2;
  }
  
  return x;
} // END AN





// [[Rcpp::export]]
double logf_w1 (
    const double& x,
    const double& P,
    const double& m1,
    const double& m2,
    const double& sig1,
    const double& sig2,
    const double& gam1, 
    const double& gam2,
    const int&    T
) {
  double out = 0;
  out += -2.0 * log( fabs(x) );
  out += log( P * pow(sig1, -0.5) * exp( -0.5 / sig1 * ( sqrt( gam2 / (T - 2) ) / x - m1 ) * ( sqrt( gam2 / (T - 2) ) / x - m1 ) )
                + (1.0 - P) * pow(sig2, -0.5) * exp( -0.5 / sig2 * ( sqrt( gam2 / (T - 2) ) / x - m2 ) * ( sqrt( gam2 / (T - 2) ) / x - m2 ) ) );
  return out;
} // END logf_w1



// [[Rcpp::export]]
double logtar_w1 (
    const double& x,
    const double& gam1, 
    const double& gam2,
    const double& invVw10, 
    const double& w10,
    const int&    T
) {
  double out = 0;
  out -= T * log( fabs(x) );
  out -= 0.5 * ( gam1 / x + gam2 / (x * x) );
  out -= 0.5 * invVw10 * (x - w10) * (x - w10);
  return out;
} // END logtar_w1


// [[Rcpp::export]]
double exact_sampling_w1 (
    const double& gam1, 
    const double& gam2, 
    const int&    T, 
    const double& neww1, 
    const double& w1, 
    const double& invVw10, 
    const double& w10
) {
  
  double a = - 0.5 * gam1 / sqrt( (T - 2) * gam2 );
  double b = 1.0 / (T - 2);
  double m1 = 0.5 * a - 0.5 * sqrt( a * a + 4.0 );
  double m2 = 0.5 * a + 0.5 * sqrt( a * a + 4.0 );
  double sig1 = m1 * m1 / (1.0 + m1 * m1) * b;
  double sig2 = m2 * m2 / (1.0 + m2 * m2) * b;
  double P = 1.0 / (1.0 + exp( 2.0 * a / b ));
  
  double logf_w1_neww1 = logf_w1( neww1, P, m1, m2, sig1, sig2, gam1, gam2, T );
  double logf_w1_w1    = logf_w1( w1,   P, m1, m2, sig1, sig2, gam1, gam2 , T );
  double logtar_w1_neww1 = logtar_w1( neww1, gam1, gam2, invVw10, w10, T );
  double logtar_w1_w1    = logtar_w1( w1,   gam1, gam2, invVw10, w10, T );
  double out = logtar_w1_neww1 - logtar_w1_w1 + logf_w1_w1 - logf_w1_neww1;
  
  return out;
} // END exact_sampling_w1




// [[Rcpp::export]]
Rcpp::List construct_LR (
    const arma::mat& zeroB_R,
    const arma::mat& signB_R
) {
  
  const int n = zeroB_R.n_cols;
  
  field<mat> R_E(n);
  field<mat> D_E(n);
  for (int i = 0; i < n; ++i) {
    uvec zeroB_Ri = find(zeroB_R.col(i) == 0);
    mat SS = eye<mat>(n, n);
    R_E(i) = SS.rows(zeroB_Ri);
    D_E(i) = zeros<mat>(zeroB_Ri.n_elem, 1);
  }
  field<mat> R_I(n);
  field<mat> U_I(n);
  field<mat> L_I(n);
  for (int i = 0; i < n; ++i) {
    // R_I
    uvec id = find(signB_R.col(i) != 0);
    mat R_Ii = eye<mat>(n, n);
    R_I(i) = R_Ii.rows(id);
    // U_I, L_I
    mat R_Ii_mat = R_I(i);
    if (R_Ii_mat.n_rows > 0) {
      mat signRI_i = R_Ii_mat * signB_R.col(i);
      vec ui(signRI_i.n_rows);
      vec li(signRI_i.n_rows);
      
      for (size_t ii = 0; ii < signRI_i.n_rows; ++ii) {
        if (signRI_i(ii) == 1) {
          ui(ii) = datum::inf;
          li(ii) = 0;
        } else if (signRI_i(ii) == -1) {
          ui(ii) = 0;
          li(ii) = -datum::inf;
        } else if (signRI_i(ii) == 0) {
          ui(ii) = datum::inf;
          li(ii) = -datum::inf;
        }
      } 
      U_I(i) = ui;
      L_I(i) = li;
    } else {
      U_I(i) = mat(0, 1);
      L_I(i) = mat(0, 1);
    }
  }
  
  List out = List::create(
    _["R_E"] = R_E,
    _["D_E"] = D_E,
    _["R_I"] = R_I,
    _["U_I"] = U_I,
    _["L_I"] = L_I
  );
  return out;
} // END construct_LR




// [[Rcpp::export]]
arma::mat sample_Theta0_Hou24_heterosk1 (
    arma::mat&        aux_Theta0,     // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const Rcpp::List& restrictions,   // output of construct_LR
    const bool        debug = false
) {
  
  if (debug) {Rcout << "START" << endl;}
  const int N           = aux_Theta0.n_rows;
  const int TT          = Y.n_cols;
  const double T        = TT;
  
  mat prior_VB0         = as<mat>(prior["VB0"]);
  mat prior_B0          = as<mat>(prior["B0"]);
  
  field<mat> R_E        = as<field<mat>>(restrictions["R_E"]);
  field<mat> D_E        = as<field<mat>>(restrictions["D_E"]);
  field<mat> R_I        = as<field<mat>>(restrictions["R_I"]);
  field<mat> U_I        = as<field<mat>>(restrictions["U_I"]);
  field<mat> L_I        = as<field<mat>>(restrictions["L_I"]);
  
  mat shocks            = Y - aux_A * X;
  mat Y_XA;
  
  for (int n=0; n<N; n++) {
    if (debug) {Rcout << " iteration: " << n << endl;}
    
    Y_XA                = trans(shocks.each_row() / aux_sigma.row(n));
    
    mat Bi              = aux_Theta0;
    mat Sci             = eye(N, N);
    
    if (n > 0) {
    
      if (n == N - 1) {
        Sci             = join_horiz(Sci.col(n), Sci.cols(0, n - 1));
      } else {
        Sci             = join_horiz(Sci.col(n), Sci.cols(0, n - 1), Sci.cols(n + 1, N - 1));
      }
      Bi                = aux_Theta0 * Sci;
    }
    
    int nri             = D_E(n).n_rows;
    
    // prior for b1 
    if (debug) {Rcout << " prior for b1 " << endl;}
    mat Vb10            = diagmat(prior_VB0.col(n));
    vec b10             = prior_B0.col(n);
  
    // construct T, Btilde, Ytilde, w0 and invVw0, invBtilde22
    if (debug) {Rcout << " construct T, Btilde, ..." << endl;}
    mat nullBi_1        = trans( null( trans(Bi.cols(1, N - 1)) ) );
    mat Ti;
    
    if (R_E(n).n_rows == 0) {
      Ti                = join_vert( nullBi_1, trans(Bi.cols(1, N - 1)) );
    } else {
      mat nBR_tmp       = join_vert( nullBi_1, R_E(n) );
      Ti                = join_vert( nullBi_1, trans(null(nBR_tmp)), R_E(n));
    }
    mat Btilde          = Ti * Bi;
    mat Ytilde          = Y_XA * Ti.t();
    mat CStilde         = kron(Ytilde, eye(N, N));
    vec w0              = Ti * b10;
    mat invVw0          = inv_sympd(Ti * Vb10 * Ti.t());
    mat invBtilde22     = inv(Btilde.submat(1, 1, N - 1, N - 1));
 
    // inequality constraints
    if (debug) {Rcout << " inequality constraints" << endl;}
    vec l_d, u_d, Rw1_I;
    mat Rw_11_I;
    if (R_I(n).n_rows != 0) {
      mat Rw_I          = R_I(n) * inv(Ti);
      vec li            = L_I(n);
      vec ui            = U_I(n);
      Rw1_I             = Rw_I.col(0);
      if (N - nri - 1 >= 1) {
        Rw_11_I           = Rw_I.cols(1, N - nri - 1);
      }
      if (N - nri <= N - 1) {
        mat Rw_12_I     = Rw_I.cols(N - nri, N - 1);
        l_d             = li - Rw_12_I * D_E(n);
        u_d             = ui - Rw_12_I * D_E(n);
      } else {
        l_d             = li;
        u_d             = ui;
      }
    }
    
    // initialize w (w1 and w_1)
    if (debug) {Rcout << " initialize w (w1 and w_1)" << endl;}
    vec     w           = Btilde.col(0);
    double  w1          = w(0);
    vec     w_1         = w.subvec(1, N - 1);
    
    // sample w_1
    // conditional prior
    if (debug) {Rcout << " sample w_1 conditional prior" << endl;}
    vec w_10            = w0.subvec(1, N - 1) - invVw0.submat(1, 1, N - 1, N - 1) * invVw0.submat(1, 0, N - 1, 0) * (w1 - w0(0));
    mat invVw_10        = invVw0.submat(1, 1, N - 1, N - 1);
    vec w_110; mat invVw_110;
    if (R_E(n).n_rows != 0 && N - nri - 2 >= 0) {
      w_110             = w_10.subvec(0, N - nri - 2) - invVw_10.submat(0, 0, N - nri - 2, N - nri - 2) * invVw_10.submat(0, N - nri - 1, N - nri - 2, N - 2) * (D_E(n) - w_10.subvec(N - nri - 1, N - 2));
      invVw_110         = invVw_10.submat(0, 0, N - nri - 2, N - nri - 2);
    } else {
      w_110             = w_10;
      invVw_110         = invVw_10;
    }
    
    // construct q and Z
    if (debug) {Rcout << " construct q and Z" << endl;}
    vec w1w1            = vec(1, fill::value(1 / w1));
    mat q_w_1           = CStilde * join_vert( w1w1, zeros(N - 1, 1), vectorise( join_vert( zeros(1, N - 1), invBtilde22 ) ) );
    mat Z_w_1           = CStilde * join_vert( zeros(1, N - 1), invBtilde22 / w1, zeros(N * (N - 1), N - 1) );
    
    
    // construct qtilde and Ztilde
    if (debug) {Rcout << " construct qtilde and Ztilde" << endl;}
    mat Z_w11; vec q_w11;
    if (R_E(n).n_rows != 0 && N - nri - 2 >= 0) {
      
      mat S             = eye(N - 1, N - 1);
      mat S1            = S.cols(0, N - nri - 2);     // others
      mat S2            = S.cols(N - nri - 1, N - 2); // equality constraints
      Z_w11             = Z_w_1 * S1;
      q_w11             = q_w_1 - Z_w_1 * S2 * D_E(n);
      
    } else {
      Z_w11             = Z_w_1;
      q_w11             = q_w_1;
    }
    
    // sample w_11 (or w_1)    
    if (debug) {Rcout << " sample w_11 (or w_1)" << endl;}
    mat Dw_11           = inv_sympd(Z_w11.t() * Z_w11 + invVw_110);
    Dw_11               = 0.5 * (Dw_11 + Dw_11.t());
    vec w_11hat         = Dw_11 * (Z_w11.t() * q_w11 + invVw_110 * w_110);
    vec w_11, neww_11;
    
    if (R_I(n).n_rows != 0 && N - nri - 2 >= 0) {
    
      // inequality constaints for w11 conditional on w1, i.e., lw_11 < Rw_11_I*w_11 < uw_11 
      vec lw_11         = l_d - Rw1_I * w1;
      vec uw_11         = u_d - Rw1_I * w1;
      w_11              = w_1.subvec(0, N - nri - 2);   
      
      // This is quite some change! Not using Botev's algo at all!
      w_11              = Rodriguez_Yam_2004(w_11hat, Dw_11, Rw_11_I, lw_11, uw_11, w_11);
      
    } else {
      w_11              = w_11hat + chol(Dw_11, "lower") * randn(w_11hat.n_elem);
    }
    if (N - nri - 2 >= 0) {
      w_1                 = join_vert(w_11, D_E(n));
    } else {
      w_1                 = D_E(n);
    }
    
    // sample w1
    // conditional prior 
    if (debug) {Rcout << " conditional prior" << endl;}
    double invVw10      = invVw0(0, 0);
    double w10          = w0(0) - as_scalar(invVw10 * invVw0.submat(0, 1, 0, N - 1) * (w_1 - w0.subvec(1, N - 1)));
    
    // construct gam1 and gam2
    if (debug) {Rcout << " construct gam1 and gam2" << endl;}
    vec q_w1           = CStilde * join_vert( zeros(N, 1), vectorise( join_vert( zeros(1, N - 1), invBtilde22 ) ) );
    vec Z_w1           = CStilde * join_vert( vec(1, fill::ones), -invBtilde22 * w_1, zeros(N * ( N - 1), 1) );
    double gam1        = 2 * as_scalar(Z_w1.t() * q_w1);
    double gam2        = as_scalar(Z_w1.t() * Z_w1);
    
    // sample w1
    if (debug) {Rcout << " sample w1" << endl;}
    double w1tilde      = AN ( -0.5 * gam1 / sqrt((T - 2) * gam2), 1/(T - 2) );
    double neww1        = sqrt(gam2 / (T - 2)) / w1tilde;
    double logAP        = exact_sampling_w1 ( gam1, gam2, TT, neww1, w1, invVw10, w10 );
    
    if (R_I(n).n_rows != 0 && N - nri - 2 >= 0) {
      // inequality constraint for w1 conditional on w_1, i.e., lw1 < Rw1_I*w1 < uw1
      vec lw1           = l_d - Rw_11_I * w_11;
      vec uw1           = u_d - Rw_11_I * w_11;
      if (log(randu()) < logAP && prod( Rw1_I * neww1 < uw1 ) == 1 && prod( lw1 < Rw1_I * neww1 ) == 1 ) {
        w1              = neww1;
      }
    } else if (log(randu()) < logAP) {
      w1                = neww1;
    }
    vec w1_vec(1); w1_vec(0) = w1;
    w                   = join_vert(w1_vec, w_1);
    aux_Theta0.col(n)   = solve(Ti, w);
    
  } // END n loop
  
  return aux_Theta0;
} // END sample_Theta0_Hou24_heterosk1



// [[Rcpp::export]]
arma::mat sample_Theta0_Hou24_heterosk1_coln (
    const int         n,
    arma::mat&        aux_Theta0,     // NxN
    const arma::mat&  shocks,         // NxT B(Y-AX)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::mat&  R_E_n,
    const arma::mat&  D_E_n,
    const Rcpp::List& restrictions,   // containing only R_E(n) and D_E(n)
    const bool        debug = false
) {
  
  if (debug) {Rcout << "START" << endl;}
  const int N           = aux_Theta0.n_rows;
  const int TT          = shocks.n_cols;
  const double T        = TT;
  
  mat prior_VB0         = as<mat>(prior["VB0"]);
  mat prior_B0          = as<mat>(prior["B0"]);
  
  field<mat> R_I        = as<field<mat>>(restrictions["R_I"]);
  field<mat> U_I        = as<field<mat>>(restrictions["U_I"]);
  field<mat> L_I        = as<field<mat>>(restrictions["L_I"]);
  
  mat Y_XA;
  
  if (debug) {Rcout << " iteration: " << n << endl;}
  
  Y_XA                = trans(shocks.each_row() / aux_sigma.row(n));
  
  mat Bi              = aux_Theta0;
  mat Sci             = eye(N, N);
  
  if (n > 0) {
    
    if (n == N - 1) {
      Sci             = join_horiz(Sci.col(n), Sci.cols(0, n - 1));
    } else {
      Sci             = join_horiz(Sci.col(n), Sci.cols(0, n - 1), Sci.cols(n + 1, N - 1));
    }
    Bi                = aux_Theta0 * Sci;
  }
  
  int nri             = D_E_n.n_rows;
  
  // prior for b1 
  if (debug) {Rcout << " prior for b1 " << endl;}
  mat Vb10            = diagmat(prior_VB0.col(n));
  vec b10             = prior_B0.col(n);
  
  // construct T, Btilde, Ytilde, w0 and invVw0, invBtilde22
  if (debug) {Rcout << " construct T, Btilde, ..." << endl;}
  mat nullBi_1        = trans( null( trans(Bi.cols(1, N - 1)) ) );
  mat Ti;
  
  if (R_E_n.n_rows == 0) {
    Ti                = join_vert( nullBi_1, trans(Bi.cols(1, N - 1)) );
  } else {
    mat nBR_tmp       = join_vert( nullBi_1, R_E_n );
    Ti                = join_vert( nullBi_1, trans(null(nBR_tmp)), R_E_n);
  }
  mat Btilde          = Ti * Bi;
  mat Ytilde          = Y_XA * Ti.t();
  mat CStilde         = kron(Ytilde, eye(N, N));
  vec w0              = Ti * b10;
  mat invVw0          = inv_sympd(Ti * Vb10 * Ti.t());
  mat invBtilde22     = inv(Btilde.submat(1, 1, N - 1, N - 1));
  
  // inequality constraints
  if (debug) {Rcout << " inequality constraints" << endl;}
  vec l_d, u_d, Rw1_I;
  mat Rw_11_I;
  if (R_I(n).n_rows != 0) {
    mat Rw_I          = R_I(n) * inv(Ti);
    vec li            = L_I(n);
    vec ui            = U_I(n);
    Rw1_I             = Rw_I.col(0);
    if (N - nri - 1 >= 1) {
      Rw_11_I           = Rw_I.cols(1, N - nri - 1);
    }
    if (N - nri <= N - 1) {
      mat Rw_12_I     = Rw_I.cols(N - nri, N - 1);
      l_d             = li - Rw_12_I * D_E_n;
      u_d             = ui - Rw_12_I * D_E_n;
    } else {
      l_d             = li;
      u_d             = ui;
    }
  }
  
  // initialize w (w1 and w_1)
  if (debug) {Rcout << " initialize w (w1 and w_1)" << endl;}
  vec     w           = Btilde.col(0);
  double  w1          = w(0);
  vec     w_1         = w.subvec(1, N - 1);
  
  // sample w_1
  // conditional prior
  if (debug) {Rcout << " sample w_1 conditional prior" << endl;}
  vec w_10            = w0.subvec(1, N - 1) - invVw0.submat(1, 1, N - 1, N - 1) * invVw0.submat(1, 0, N - 1, 0) * (w1 - w0(0));
  mat invVw_10        = invVw0.submat(1, 1, N - 1, N - 1);
  vec w_110; mat invVw_110;
  if (R_E_n.n_rows != 0 && N - nri - 2 >= 0) {
    w_110             = w_10.subvec(0, N - nri - 2) - invVw_10.submat(0, 0, N - nri - 2, N - nri - 2) * invVw_10.submat(0, N - nri - 1, N - nri - 2, N - 2) * (D_E_n - w_10.subvec(N - nri - 1, N - 2));
    invVw_110         = invVw_10.submat(0, 0, N - nri - 2, N - nri - 2);
  } else {
    w_110             = w_10;
    invVw_110         = invVw_10;
  }
  
  // construct q and Z
  if (debug) {Rcout << " construct q and Z" << endl;}
  vec w1w1            = vec(1, fill::value(1 / w1));
  mat q_w_1           = CStilde * join_vert( w1w1, zeros(N - 1, 1), vectorise( join_vert( zeros(1, N - 1), invBtilde22 ) ) );
  mat Z_w_1           = CStilde * join_vert( zeros(1, N - 1), invBtilde22 / w1, zeros(N * (N - 1), N - 1) );
  
  
  // construct qtilde and Ztilde
  if (debug) {Rcout << " construct qtilde and Ztilde" << endl;}
  mat Z_w11; vec q_w11;
  if (R_E_n.n_rows != 0 && N - nri - 2 >= 0) {
    
    mat S             = eye(N - 1, N - 1);
    mat S1            = S.cols(0, N - nri - 2);     // others
    mat S2            = S.cols(N - nri - 1, N - 2); // equality constraints
    Z_w11             = Z_w_1 * S1;
    q_w11             = q_w_1 - Z_w_1 * S2 * D_E_n;
    
  } else {
    Z_w11             = Z_w_1;
    q_w11             = q_w_1;
  }
  
  // sample w_11 (or w_1)    
  if (debug) {Rcout << " sample w_11 (or w_1)" << endl;}
  mat Dw_11           = inv_sympd(Z_w11.t() * Z_w11 + invVw_110);
  Dw_11               = 0.5 * (Dw_11 + Dw_11.t());
  vec w_11hat         = Dw_11 * (Z_w11.t() * q_w11 + invVw_110 * w_110);
  vec w_11, neww_11;
  
  if (R_I(n).n_rows != 0 && N - nri - 2 >= 0) {
    
    // inequality constaints for w11 conditional on w1, i.e., lw_11 < Rw_11_I*w_11 < uw_11 
    vec lw_11         = l_d - Rw1_I * w1;
    vec uw_11         = u_d - Rw1_I * w1;
    w_11              = w_1.subvec(0, N - nri - 2);   
    
    // This is quite some change! Not using Botev's algo at all!
    w_11              = Rodriguez_Yam_2004(w_11hat, Dw_11, Rw_11_I, lw_11, uw_11, w_11);
    
  } else {
    w_11              = w_11hat + chol(Dw_11, "lower") * randn(w_11hat.n_elem);
  }
  if (N - nri - 2 >= 0) {
    w_1                 = join_vert(w_11, D_E_n);
  } else {
    w_1                 = D_E_n;
  }
  
  // sample w1
  // conditional prior 
  if (debug) {Rcout << " conditional prior" << endl;}
  double invVw10      = invVw0(0, 0);
  double w10          = w0(0) - as_scalar(invVw10 * invVw0.submat(0, 1, 0, N - 1) * (w_1 - w0.subvec(1, N - 1)));
  
  // construct gam1 and gam2
  if (debug) {Rcout << " construct gam1 and gam2" << endl;}
  vec q_w1           = CStilde * join_vert( zeros(N, 1), vectorise( join_vert( zeros(1, N - 1), invBtilde22 ) ) );
  vec Z_w1           = CStilde * join_vert( vec(1, fill::ones), -invBtilde22 * w_1, zeros(N * ( N - 1), 1) );
  double gam1        = 2 * as_scalar(Z_w1.t() * q_w1);
  double gam2        = as_scalar(Z_w1.t() * Z_w1);
  
  // sample w1
  if (debug) {Rcout << " sample w1" << endl;}
  double w1tilde      = AN ( -0.5 * gam1 / sqrt((T - 2) * gam2), 1/(T - 2) );
  double neww1        = sqrt(gam2 / (T - 2)) / w1tilde;
  double logAP        = exact_sampling_w1 ( gam1, gam2, TT, neww1, w1, invVw10, w10 );
  
  if (R_I(n).n_rows != 0 && N - nri - 2 >= 0) {
    // inequality constraint for w1 conditional on w_1, i.e., lw1 < Rw1_I*w1 < uw1
    vec lw1           = l_d - Rw_11_I * w_11;
    vec uw1           = u_d - Rw_11_I * w_11;
    if (log(randu()) < logAP && prod( Rw1_I * neww1 < uw1 ) == 1 && prod( lw1 < Rw1_I * neww1 ) == 1 ) {
      w1              = neww1;
    }
  } else if (log(randu()) < logAP) {
    w1                = neww1;
  }
  vec w1_vec(1); w1_vec(0) = w1;
  w                   = join_vert(w1_vec, w_1);
  vec aux_Theta0n     = solve(Ti, w);
  aux_Theta0.col(n)   = aux_Theta0n / aux_Theta0n(n);
  
  return aux_Theta0;
} // END sample_Theta0_Hou24_heterosk1_coln







// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_Theta0_Hou24_heterosk1_s4 (
    arma::mat                     aux_Theta0,     // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators aux_SL.slice(1).col(m)
    const arma::mat&              shocks,         // NxT B(Y-AX)
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const Rcpp::List&             restrictions    // output of construct_LR
) {
  // the function draws new values of aux_Theta0 and aux_SL
  
  field<mat> R_E        = as<field<mat>>(restrictions["R_E"]);
  field<mat> D_E        = as<field<mat>>(restrictions["D_E"]);
  field<mat> R_I        = as<field<mat>>(restrictions["R_I"]);
  field<mat> U_I        = as<field<mat>>(restrictions["U_I"]);
  field<mat> L_I        = as<field<mat>>(restrictions["L_I"]);
  
  mat prior_VB0         = as<mat>(prior["VB0"]);
  mat prior_B0          = as<mat>(prior["B0"]);
  
  const int N           = aux_Theta0.n_rows;
  const int T           = shocks.n_cols;
  int       Ltmp        = R_E.n_elem - 1;
  vec       Lm          = R_E(Ltmp);
  double    L           = accu(Lm);
  field<mat> RE         = R_E.rows(0, L-1);
  mat aux_Theta0_tmp_inv;
  
  for (int n=0; n<N; n++) {
    
    mat aux_Theta0_nL(N, Lm(n));
    vec log_posterior_kernel_nL(Lm(n));
    mat aux_Theta0_tmp      = aux_Theta0;
    
    for (int l=0; l<Lm(n); l++) {
      int ll = 0;
      if (n == 0) {
        ll                    = l;
      } else {
        vec Lm_cs             = cumsum(Lm);
        ll                    = Lm_cs(n-1) + l;
      }
    
      aux_Theta0_tmp          = sample_Theta0_Hou24_heterosk1_coln( n, aux_Theta0_tmp, shocks, aux_sigma, prior, RE(ll), D_E(ll), restrictions );
      aux_Theta0_nL.col(l)    = aux_Theta0_tmp.col(n);
      
      // posterior kernel
      aux_Theta0_tmp_inv          = inv(aux_Theta0_tmp);
      log_posterior_kernel_nL(l)  = accu(square((aux_Theta0_tmp_inv * shocks) / aux_sigma));                  // likelihood
      log_posterior_kernel_nL(l) += accu(square(aux_Theta0_tmp.col(n) - prior_B0.col(n)) / prior_VB0.col(n)); // prior
    } // END loop l
    
    // Sample S4 indicator
    int     index_s4          = 0;
    if (Lm(n) > 1) {
      // Compute S4 components probabilities
      log_posterior_kernel_nL -= log_posterior_kernel_nL.max();
      vec     pr_s4           = exp(log_posterior_kernel_nL)/accu(exp(log_posterior_kernel_nL));
      
      // Sample S4 indicator
      NumericVector seq_1S    = wrap(seq_len(Lm(n)) - 1);
      index_s4                = csample_num1(seq_1S, wrap(pr_s4));
    }
    
    aux_SL(n)                 = index_s4;
    aux_Theta0.col(n)         = aux_Theta0_nL.col(index_s4);
    
  } // END n loop
  
  return List::create(
    _["aux_Theta0"] = aux_Theta0,
    _["aux_SL"]     = aux_SL
  );
} // END sample_Theta0_Hou24_heterosk1_s4





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_Theta0_mss_s4 (
    arma::cube        aux_Theta0,     // NxNxM
    arma::imat        aux_SL,         // NxM row-specific S4 indicators
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  shocks,         // NxT
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  aux_xi,         // MxT
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const Rcpp::List& restrictions
) {
  const int T     = shocks.n_cols;
  const int N     = shocks.n_rows;
  const int M     = aux_Theta0.n_slices;
  const vec T_m  = sum(aux_xi,1);
  
  for (int m=0; m<M; m++) {
    int ii = 0;
    mat aux_sigma_m(N, T_m(m));
    mat shocks_m(N, T_m(m));
    
    for (int t=0; t<T; t++){
      if (aux_xi(m,t)==1) {
        aux_sigma_m.col(ii) = aux_sigma.col(t);
        shocks_m.col(ii)    = aux_B.slice(m) * shocks.col(t);
        ii++;
      }
    }
    
    List BSL_m              = sample_Theta0_Hou24_heterosk1_s4( aux_Theta0.slice(m), aux_SL.col(m), shocks_m, aux_sigma_m, prior, restrictions );
    aux_Theta0.slice(m)     = as<mat>(BSL_m["aux_Theta0"]);
    aux_SL.col(m)           = as<ivec>(BSL_m["aux_SL"]);
  }
  
  return List::create(
    _["aux_Theta0"] = aux_Theta0,
    _["aux_SL"]     = aux_SL
  );
} // END sample_B_mss_s4




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_BTheta0_tvi (
    arma::cube&                   aux_B,          // NxNxM
    arma::cube                    aux_Theta0,     // NxNxM
    arma::icube                   aux_SL,         // NxMx2 row-specific S4 indicators
    const arma::mat&              shocks,         // NxT shocks = Y - aux_A * X;
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              aux_xi,         // MxT
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    arma::field<arma::mat>        prior_precision, // (N,M)(N,N)
    const arma::field<arma::mat>& VB, // restrictions on B0
    const Rcpp::List&             VTheta0
) {
  const int T         = shocks.n_cols;
  const int N         = shocks.n_rows;
  const int M         = aux_Theta0.n_slices;
  const vec T_m       = sum(aux_xi,1);
  imat aux_SL_B       = aux_SL.slice(0);
  imat aux_SL_Theta0  = aux_SL.slice(1);
  cube aux_Theta0_inv(N, N, M);
  cube aux_struc(N, N, M);
  
  for (int m=0; m<M; m++) {
    int ii = 0;
    mat aux_sigma_m(N, T_m(m));
    mat shocks_m(N, T_m(m));
    
    for (int t=0; t<T; t++){
      if (aux_xi(m,t)==1) {
        aux_sigma_m.col(ii) = aux_sigma.col(t);
        shocks_m.col(ii)    = aux_B.slice(m) * shocks.col(t);
        ii++;
      }
    }
    
    List Theta0SL_m         = sample_Theta0_Hou24_heterosk1_s4( aux_Theta0.slice(m), aux_SL_Theta0.col(m), shocks_m, aux_sigma_m, prior, VTheta0 );
    aux_Theta0.slice(m)     = as<mat>(Theta0SL_m["aux_Theta0"]);
    aux_SL_Theta0.col(m)    = as<ivec>(Theta0SL_m["aux_SL"]);
    aux_Theta0_inv.slice(m) = inv(aux_Theta0.slice(m));
    
    List BSL_m              = sample_B_heterosk1_s4( aux_B.slice(m), aux_SL_B.col(m), aux_Theta0.slice(m), aux_Theta0_inv.slice(m), shocks_m, aux_sigma_m, prior_precision.col(m), prior, VB );
    aux_B.slice(m)          = as<mat>(BSL_m["aux_B"]);
    aux_SL_B.col(m)         = as<ivec>(BSL_m["aux_SL"]);
    
    aux_struc.slice(m)      = aux_Theta0_inv.slice(m) * aux_B.slice(m);
    
  }
  aux_SL.slice(0)       = aux_SL_B;
  aux_SL.slice(1)       = aux_SL_Theta0;
  
  return List::create(
    _["aux_B"]          = aux_B,
    _["aux_Theta0"]     = aux_Theta0,
    _["aux_struc"]      = aux_struc,
    _["aux_SL"]         = aux_SL
  );
} // END sample_BTheta0_tvi
