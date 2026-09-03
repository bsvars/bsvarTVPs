
#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"

#include "sample_sv_ms.h"

using namespace Rcpp;
using namespace arma;



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double rig_inv1 (
    double alpha,
    double beta
) {
  double out = randg( distr_param(alpha, pow(beta, -1)) );
  return out;
} // END rig_inv1




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec mvnrnd_prec_cond (
    arma::vec x,          // Nx1 vector to be filled with conditional normal draws when missing == 1
    arma::vec mu,         // Nx1 mean vector
    arma::mat precision,  // NxN precision matrix
    arma::vec to_sample   // Nx1 with 1 for missing observations
) {
  
  int   N         = x.n_elem;
  uvec  ind       = find( to_sample == 0 );
  uvec  ind_miss  = find(to_sample);
  mat   aj        = eye(N, N);
  
  vec   x2        = x(ind); 
  
  vec   mu1       = mu(ind_miss);
  vec   mu2       = mu(ind);
  mat   prec11    = precision(ind_miss, ind_miss);
  prec11          = 0.5 * ( prec11 + prec11.t() );
  mat   prec12    = precision(ind_miss, ind);
   
  mat   prec11_inv = inv_sympd(prec11);
  vec   mu_cond   = mu1 + prec11_inv * prec12 * (x2 - mu2);
  vec   draw      = mvnrnd( mu_cond, prec11_inv );
  
  vec   out       = aj.cols(ind_miss) * draw + aj.cols(ind) * x2;
  return out;
} // END mvnrnd_prec_cond





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat orthogonal_complement_matrix_TW (const arma::mat& x) {
  // # x is a mxn matrix and m>n
  // # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
  // # t(x)%*%out = 0 and det(cbind(x,out))!=0
  int n_nrow     = x.n_rows;
  int n_ncol     = x.n_cols;
  mat Q;
  mat R;
  qr(Q, R, x);
  mat ocm = Q.tail_cols(n_nrow-n_ncol);
  return ocm;
} // END orthogonal_complement_matrix_TW






// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_posterior_kernel_B (
    arma::mat                     aux_B,          // NxN
    arma::mat                     aux_Theta0_inv,
    arma::mat&                    shocks,         // NxT RF error terms
    arma::mat&                    aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::field<arma::mat>        prior_precision // (N)(N,N)
) {
  
  double abs_log_det_B;
  double sign_;
  log_det(abs_log_det_B, sign_, aux_B); 
  
  int       T     = aux_sigma.n_cols;
  const int N     = aux_B.n_rows;
  mat bdiag_prior_precision(N*N, N*N);
  
  for (int n=0; n<N; n++) {
    mat pp  = prior_precision(n);
    int first = 0+n*N;
    int last  = N-1+n*N;
    bdiag_prior_precision.submat(first, first, last, last) = pp;
  }
  mat S           = kron(aux_Theta0_inv, shocks.t());
  S               = S.each_col() / vectorise(aux_sigma.t());
  vec aux_Bt_vec  = vectorise(aux_B.t());
  double out      = - 0.5 * as_scalar( aux_Bt_vec.t() *(S.t() * S + bdiag_prior_precision) * aux_Bt_vec );
  out            += T * abs_log_det_B;
  
  return out;
} // END log_posterior_kernel_B





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_B_heterosk1_rown (
    int               n,              // row n of B to be sampled
    arma::mat         aux_B,          // NxN
    arma::mat&        shocks,         // NxT conditional STANDARD DEVIATIONS
    arma::mat&        aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::mat         prior_precision_n, // (N,N)
    const Rcpp::List& prior,          // a list of priors - original dimensions
    arma::mat&        VB_n        // restrictions on B0
) {
  
  bool debug = false;
  
  if ( debug ) Rcout<<" sample Bn pre "<<endl;
  const int N               = aux_B.n_rows;
  const int T               = shocks.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
    
  // Rcout << " structural equation: " << n + 1 << endl;
  // set scale matrix
  if ( debug ) Rcout<<" sample Bn sho "<<endl;
  mat shocks_sigma(N, T);
  mat ss(N, N);
  if ( T != 0 ) {
    shocks_sigma       += shocks.each_row() / aux_sigma.row(n);
    ss                 += shocks_sigma * shocks_sigma.t();
  }
  if ( debug ) Rcout<<" sample Bn SSi "<<endl;
  mat posterior_SS_inv    = prior_precision_n + ss;
  mat posterior_S_inv     = VB_n * posterior_SS_inv * VB_n.t();
  posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
  
  // sample B
  mat posterior_S(size(posterior_S_inv));
  mat Un(size(posterior_S_inv));
  
  if ( debug ) Rcout<<" sample Bn inv "<<endl;
  posterior_S           = inv_sympd(posterior_S_inv);
  
  if ( debug ) Rcout<<" sample Bn cho "<<endl;
  Un                    = chol(posterior_nu * posterior_S);
  
  mat B_tmp               = aux_B;
  B_tmp.shed_row(n);
  
  if ( debug ) Rcout<<" sample Bn ocm "<<endl;
  mat w1;
  rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
  vec w1_tmp              = trans(w * VB_n.t() * Un.t());
  double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
  w1                      = w1_tmp.t()/sqrt(w1w1_tmp);
  
  mat Wn;
  const int rn            = VB_n.n_rows;
  if (rn==1) {
    Wn                    = w1;
  } else {
    Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
  }

    if ( debug ) Rcout<<" sample Bn fin "<<endl;
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
  rowvec Brown            = b0n * VB_n;
  if (!Brown.has_nan()) {
    aux_B.row(n)           = Brown; 
  }
  
  return aux_B;
} // END sample_B_heterosk1_rown


    
    
    

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_B_heterosk1_s4 (
    arma::mat                     aux_B,          // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators
    arma::mat                     aux_Theta0,     // NxN
    arma::mat                     aux_Theta0_inv,
    arma::vec                     aux_SLlpr,
    arma::mat&                    shocks,         // NxT RF error terms
    arma::mat&                    aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    arma::field<arma::mat>        prior_precision, // (N)(N,N)
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VBL             // restrictions on B0 in S4 arrangement
) {
  
  bool debug = false;
  
  // the function draws new values of aux_B and aux_SL
  if ( debug ) Rcout<<" sample B pre "<<endl;
  const int N           = aux_B.n_rows;
  int         Ltmp      = VBL.n_elem - 1;
  vec         Lm        = VBL(Ltmp);
  double      L         = accu(Lm);
  field<mat>  VB        = VBL.rows(0, L-1);
  vec         Sdiag     = diagvec(aux_Theta0 * aux_Theta0.t()); // to be used for approximate covariance rescaling
  double      lpk_new, lpk_old;
  
  for (int n=0; n<N; n++) {
    
    if ( debug ) Rcout<<" sample B n: "<<n<<endl;
    
    mat aux_B_nL(Lm(n), N);
    mat aux_B_tmp           = aux_B;
    vec pr_s4(Lm(n), fill::value(1/Lm(n)));
    // mat shocks_rescaled     = shocks / Sdiag(n);
    vec lpks_new(Lm(n));
    if ( debug ) Rcout<<" sample B lpk "<<endl;
    lpk_old                 = aux_SLlpr(n) + log_posterior_kernel_B( aux_B, aux_Theta0_inv, shocks, aux_sigma, prior_precision ); 
  
    for (int l=0; l<Lm(n); l++) {
      if ( debug ) Rcout<<" sample B l: "<<l<<endl;
      int ll = 0;
      if (n == 0) {
        ll                    = l;
      } else {
        vec Lm_cs             = cumsum(Lm);
        ll                    = Lm_cs(n-1) + l;
      }
      try {
        if ( debug ) Rcout<<" sample B before: "<<l<<endl;
        aux_B_tmp               = sample_B_heterosk1_rown( n, aux_B, shocks, aux_sigma, prior_precision(n), prior, VB(ll) );
        if ( debug ) Rcout<<" sample B after: "<<l<<endl;
      } catch (std::runtime_error &e) {}
      
      aux_B_nL.row(l)         = aux_B_tmp.row(n);
      lpks_new(l)             = log_posterior_kernel_B( aux_B, aux_Theta0_inv, shocks, aux_sigma, prior_precision );
    } // END loop l
    
    // Sample S4 indicator
    if ( debug ) Rcout<<" sample B S4 "<<endl;
    int index_s4              = 0;
    if (Lm(n) > 1) {
      // Compute S4 components probabilities
      vec lpks_norm           = lpks_new - lpks_new.max();
      pr_s4                   = exp(lpks_norm)/accu(exp(lpks_norm));
    
      // Sample S4 indicator
      NumericVector seq_1S    = wrap(seq_len(Lm(n)) - 1);
      index_s4                = csample_num1(seq_1S, wrap(pr_s4));
    }
    
    lpk_new                   = lpks_new(index_s4);
    double lprs               = log(pr_s4(index_s4)); 
    // MH gate
    if (exp((lpk_new + lprs) - lpk_old) > randu()) {
      aux_SL(n)                 = index_s4;
      aux_SLlpr(n)              = lprs;
      aux_B.row(n)              = aux_B_nL.row(index_s4);
    }
    
  } // END n loop
  
  return List::create(
    _["aux_B"]    = aux_B,
    _["aux_SL"]   = aux_SL,
    _["aux_SLlpr"] = aux_SLlpr
  );
} // END sample_B_heterosk1_s4




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_A_heterosk1_mss (
    arma::mat         aux_A,          // NxK
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    arma::field<arma::mat>  prior_precision, // (N)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VA  // restrictions on A
) {
  
  bool debug = false;
  
  // the function draws the value of aux_A
  if ( debug ) Rcout << " 1" << endl;
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  const int T         = Y.n_cols;
  const int M         = aux_xi.n_rows;
  const vec Tm       = sum(aux_xi, 1);
  
  mat prior_A_mean    = as<mat>(prior["A"]);
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
  
  if ( debug ) Rcout << " 2" << endl;
  for (uword t=0; t<T; t++) {
    int mm            = aux_xi.col(t).index_max();
    YM(mm).col(m_iter(mm))    = Y.col(t);
    XM(mm).col(m_iter(mm))    = X.col(t);
    SM(mm).col(m_iter(mm))    = aux_sigma.col(t);
    m_iter += aux_xi.col(t);
  }
  
  if ( debug ) Rcout << " 3" << endl;
  for (int n=0; n<N; n++) {
    if ( debug ) Rcout << " n: " << n << endl;
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B.slice(0) * (YM(0) - A0 * XM(0)) );
    mat   Wn          = kron( trans(XM(0)), aux_B.slice(0).col(n) );
    vec   Sn          = vectorise(SM(0));
    
    if ( debug ) Rcout << " 4" << endl;
    for (uword m=1; m<M; m++) {
      zn          = join_cols(zn, vectorise(aux_B.slice(m) * (YM(m) - A0 * XM(m))));
      Wn          = join_cols(Wn, kron(trans(XM(m)), aux_B.slice(m).col(n)));
      Sn          = join_cols(Sn, vectorise(SM(m)));
    }
    
    mat   zn_sigma    = zn / Sn;
    mat   Wn_sigma    = Wn.each_col() / Sn;
    
    if ( debug ) Rcout << " 5" << endl;
    mat     precision = prior_precision(n) + trans(Wn_sigma) * Wn_sigma;
    precision         = 0.5 * (precision + precision.t());
    rowvec  location  = prior_A_mean.row(n) * prior_precision(n) + trans(zn_sigma) * Wn_sigma;
    vec     mean      = solve(precision, location.t(), solve_opts::likely_sympd);
    
    vec     draw      = mvnrnd_prec_cond ( 
                          trans(aux_A.row(n)),
                          mean,
                          precision,
                          trans(sum(VA(n)))
    );
    
    aux_A.row(n)      = trans(draw);
  } // END n loop
  
  return aux_A;
} // END sample_A_heterosk1_mss





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube sample_A_heterosk1_mssa (
    arma::cube        aux_A,          // NxKxM
    const arma::cube& aux_B,          // NxNxM
    const arma::mat&  aux_xi,         // MxT
    arma::field<arma::mat>  prior_precision, // (N,M)(N,N)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VA  // restrictions on A
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  const int T         = Y.n_cols;
  const int M         = aux_xi.n_rows;
  const vec Tm        = sum(aux_xi, 1);
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  
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
      
      mat     precision = prior_precision(n, m) + trans(Wn_sigma) * Wn_sigma;
      precision         = 0.5 * (precision + precision.t());
      rowvec  location  = prior_A_mean.row(n) * prior_precision(n, m) + trans(zn_sigma) * Wn_sigma;
      vec     mean      = solve(precision, location.t(), solve_opts::likely_sympd);
      
      vec     draw      = mvnrnd_prec_cond ( 
        trans(aux_A.slice(m).row(n)),
        mean,
        precision,
        trans(sum(VA(n)))
      );
      
      aux_A.slice(m).row(n) = trans(draw);
    } // END n loop
  } // END m loop
  
  return aux_A;
} // END sample_A_heterosk1_mssa



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_hyperparameters_mss_s4_boost (
    Rcpp::List&             aux_hyper_list,      // (2 * N + 1) x 2
    const arma::cube&       aux_B,            // NxNxM
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
) {
  // the function changes the value of aux_hyper by reference
  mat       aux_hyper = as<mat>(aux_hyper_list["aux_hyper"]); 
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
  
  return List::create(
    _["aux_hyper"] = aux_hyper
  );
} // END sample_hyperparameters_mss_s4_boost



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_hyperparameters_mssa_s4_boost (
    Rcpp::List&             aux_hyper_list,      // (2 * N + 1) x 2
    const arma::cube&       aux_B,            // NxNxM
    const arma::cube&       aux_A,            // NxKxM
    const arma::field<arma::mat>& VB,
    const arma::imat&       aux_SL,         // NxM row-specific S4 indicators
    const Rcpp::List&       prior,
    const bool              hyper_boost = true
) {
  // the function changes the value of aux_hyper by reference
  mat       aux_hyper = as<mat>(aux_hyper_list["aux_hyper"]); 
  
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
  
  return List::create(
    _["aux_hyper"] = aux_hyper
  );
} // END sample_hyperparameters_mssa_s4_boost



arma::field<arma::mat> hyper2precisionB_mss_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
) {
  
  mat         aux_hyper     = as<mat>(aux_hyper_list["aux_hyper"]); 
  mat         prior_SS_inv  = as<mat>(prior["B_V_inv"]);
  mat         prior_PR_TR   = as<mat>(prior["PR_TR"]);
  
  int         N             = prior_SS_inv.n_rows; 
  int         M             = prior_PR_TR.n_rows;
  field<mat>  precisionB(N, M);
  
  for (int n=0; n<N; n++) {
    for (int m=0; m<M; m++) {
      precisionB(n, m)      = pow(aux_hyper(n, 0), -1) * prior_SS_inv;
    }
  }
  
  return precisionB;
} // END hyper2precisionB_boost


arma::field<arma::mat> hyper2precisionA_msa_boost (
    Rcpp::List              aux_hyper_list,
    const Rcpp::List&       prior
) {
  mat         aux_hyper     = as<mat>(aux_hyper_list["aux_hyper"]); 
  mat         prior_A_V_inv = as<mat>(prior["A_V_inv"]);
  mat         prior_A       = as<mat>(prior["A"]);
  mat         prior_PR_TR   = as<mat>(prior["PR_TR"]);
  
  int         N             = prior_A.n_rows; 
  int         M             = prior_PR_TR.n_rows;
  field<mat>  precisionA(N, M);
  
  for (int n=0; n<N; n++) {
    for (int m=0; m<M; m++) {
      precisionA(n, m)      = pow(aux_hyper(n, 1), -1) * prior_A_V_inv;
    }
  }  
  
  return precisionA;
} // END hyper2precisionA_msa_boost



arma::field<arma::mat> hyper2precisionA_boost (
    Rcpp::List              aux_hyper_list,      // (2 * N + 1) x 2
    const Rcpp::List&       prior
) {
  
  mat         aux_hyper     = as<mat>(aux_hyper_list["aux_hyper"]); 
  mat         prior_A_V_inv = as<mat>(prior["A_V_inv"]);
  mat         prior_A       = as<mat>(prior["A"]);
  int         N             = prior_A.n_rows; 
  field<mat>  precisionA(N);
  
  for (int n=0; n<N; n++) {
    precisionA(n)           = pow(aux_hyper(n, 1), -1) * prior_A_V_inv;
  }
  
  return precisionA;
} // END hyper2precisionA_boost



// [[Rcpp::export]]
Rcpp::List construct_LR (
    const arma::mat& zeroB_R
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
  List out = List::create(
    _["R_E"] = R_E,
    _["D_E"] = D_E
  );
  return out;
} // END construct_LR



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_posterior_kernel_Theta0 (
    const int n,
    const arma::mat& aux_Theta0,
    const arma::mat& aux_B,
    const arma::mat& shocks,
    const arma::mat& aux_sigma,
    const arma::mat& prior_B0,
    const arma::mat& prior_VB0
) {
  
  int     T = shocks.n_cols;
  double  ldet;
  double  ldet_sign;
  log_det( ldet, ldet_sign, aux_Theta0 );
  
  mat struc_shocks = solve(aux_Theta0, aux_B * shocks);
  struc_shocks.each_row() /= aux_sigma.row(n);
  
  double  lpk = - 0.5 * accu(square(struc_shocks));                   // likelihood
  lpk        -= T * ldet;                                                              // likelihood log determinant term
  lpk        -= 0.5 * accu(square(aux_Theta0 - prior_B0) / prior_VB0);  // prior
  
  return lpk;
} // END log_posterior_kernel_Theta0_n



// [[Rcpp::export]]
arma::mat sample_Theta0_Hou24_heterosk1_coln (
    const int         n,
    arma::mat&        aux_Theta0,     // NxN
    const arma::mat&  shocks,         // NxT B(Y-AX)
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::mat&  R_E_n,
    const arma::mat&  D_E_n
) {
  
  bool debug = false;
  
  if (debug) {Rcout << "START" << endl;}
  const int N           = aux_Theta0.n_rows;
  // const int TT          = Y.n_cols;
  // const double T        = TT;
  
  mat prior_VB0         = as<mat>(prior["VB0"]);
  mat prior_B0          = as<mat>(prior["B0"]);
  
  // mat shocks            = Y - aux_A * X;
  mat Y_XA;
  
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
  mat Vw0             = Ti * Vb10 * Ti.t();
  Vw0                 = 0.5 * (Vw0 + Vw0.t());
  mat invVw0          = inv_sympd(Vw0);
  mat invBtilde22     = inv(Btilde.submat(1, 1, N - 1, N - 1));
  
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
  mat Dw_11_inv       = Z_w11.t() * Z_w11 + invVw_110;
  Dw_11_inv           = 0.5 * (Dw_11_inv + Dw_11_inv.t());
  mat Dw_11           = inv_sympd(Dw_11_inv);
  Dw_11               = 0.5 * (Dw_11 + Dw_11.t());
  vec w_11hat         = Dw_11 * (Z_w11.t() * q_w11 + invVw_110 * w_110);
  vec w_11, neww_11;
  w_11                = w_11hat + chol(Dw_11, "lower") * randn(w_11hat.n_elem);
  if (N - nri - 2 >= 0) {
    w_1                 = join_vert(w_11, D_E_n);
  } else {
    w_1                 = D_E_n;
  }
  
  vec w1_vec(1); w1_vec(0) = w1;
  w                   = join_vert(w1_vec, w_1);
  vec aux_Theta0n     = solve(Ti, w);
  // aux_Theta0.col(n)   = sign(aux_Theta0n(n)) * aux_Theta0n;
  if (debug) {Rcout << " aux_Theta0n: " << endl << aux_Theta0n << endl;}
  aux_Theta0.col(n)   = aux_Theta0n / aux_Theta0n(n);
  // No MH here!
  return aux_Theta0;
} // END sample_Theta0_Hou24_heterosk1_coln




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_Theta0_Hou24_heterosk1_s4 (
    arma::mat                     aux_Theta0,     // NxN
    arma::mat                     aux_B,          // NxN
    arma::ivec                    aux_SL,         // Nx1 row-specific S4 indicators aux_SL.slice(1).col(m)
    arma::vec                     aux_SLlpr,       // N col os S4 indicators probs aux_SL.slice(1).col(m)
    const arma::mat&              shocks,         // NxT B(Y-AX)
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    const Rcpp::List&             restrictions    // output of construct_LR
) {
  // the function draws new values of aux_Theta0 and aux_SL
  
  bool debug = false;
  
  field<mat> R_E        = as<field<mat>>(restrictions["R_E"]);
  field<mat> D_E        = as<field<mat>>(restrictions["D_E"]);
  
  mat prior_VB0         = as<mat>(prior["VB0"]);
  mat prior_B0          = as<mat>(prior["B0"]);
  
  const int N           = aux_Theta0.n_rows;
  int       Ltmp        = R_E.n_elem - 1;
  vec       Lm          = R_E(Ltmp);
  double    L           = accu(Lm);
  field<mat> RE         = R_E.rows(0, L-1);
  mat aux_Theta0_tmp_inv;
  
  for (int n=0; n<N; n++) {
    if ( debug ) Rcout << " n: " << n << endl;
    
    mat     aux_Theta0_nL(N, Lm(n));
    vec     log_posterior_kernel_nL(Lm(n));
    mat     aux_Theta0_tmp    = aux_Theta0;
    double  lpk_old           = aux_SLlpr(n) + log_posterior_kernel_Theta0 ( n, aux_Theta0_tmp, aux_B, shocks, aux_sigma, prior_B0, prior_VB0 );
    
    for (int l=0; l<Lm(n); l++) {
      int ll = 0;
      if (n == 0) {
        ll                    = l;
      } else {
        vec Lm_cs             = cumsum(Lm);
        ll                    = Lm_cs(n-1) + l;
      }
      
      if ( debug ) Rcout << " before sample_Theta0_Hou24_heterosk1_coln" << endl;
      try {
        aux_Theta0_tmp        = sample_Theta0_Hou24_heterosk1_coln( n, aux_Theta0_tmp, shocks, aux_sigma, prior, RE(ll), D_E(ll) );
        aux_Theta0_nL.col(l)  = aux_Theta0_tmp.col(n);
      } catch (std::runtime_error &e) {}
      if ( debug ) Rcout << " after sample_Theta0_Hou24_heterosk1_coln" << endl;
      
      // posterior kernel
      log_posterior_kernel_nL(l)  = log_posterior_kernel_Theta0 ( n, aux_Theta0_tmp, aux_B, shocks, aux_sigma, prior_B0, prior_VB0 );
      
    } // END loop l
    
    // Sample S4 indicator
    int     index_s4          = 0;
    double  lmn = Lm(n);
    vec     pr_s4(Lm(n), fill::value(1/lmn));
    if (Lm(n) > 1) {
      // Compute S4 components probabilities
      log_posterior_kernel_nL -= log_posterior_kernel_nL.max();
      pr_s4                   = exp(log_posterior_kernel_nL)/accu(exp(log_posterior_kernel_nL));
      
      // Sample S4 indicator
      NumericVector seq_1S    = wrap(seq_len(Lm(n)) - 1);
      index_s4                = csample_num1(seq_1S, wrap(pr_s4));
    }
    
    double lprs               = log(pr_s4(index_s4));
    if ( exp((log_posterior_kernel_nL(index_s4) + lprs) - lpk_old) > randu() ) {
      aux_SL(n)                 = index_s4;
      aux_SLlpr(n)              = lprs;  
      aux_Theta0.col(n)         = aux_Theta0_nL.col(index_s4);
    }
    
  } // END n loop
  
  return List::create(
    _["aux_Theta0"] = aux_Theta0,
    _["aux_SL"]     = aux_SL,
    _["aux_SLlpr"]   = aux_SLlpr
  );
} // END sample_Theta0_Hou24_heterosk1_s4



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_BTheta0_tvi (
    arma::cube                    aux_B,          // NxNxM
    arma::cube                    aux_Theta0,     // NxNxM
    arma::icube                   aux_SL,         // NxMx2 row-specific S4 indicators
    arma::cube                    aux_SLlpr,       // NxMx2 row-specific S4 indicators for prior 
    const arma::mat&              shocks,         // NxT shocks = Y - aux_A * X;
    const arma::mat&              aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&              aux_xi,         // MxT
    const Rcpp::List&             prior,          // a list of priors - original dimensions
    arma::field<arma::mat>        prior_precision, // (N,M)(N,N)
    const arma::field<arma::mat>& VB, // restrictions on B0
    const Rcpp::List&             VTheta0
) {
  
  bool debug = false;
  
  const int T         = shocks.n_cols;
  const int N         = shocks.n_rows;
  const int M         = aux_Theta0.n_slices;
  const vec T_m       = sum(aux_xi,1);
  imat aux_SL_B       = aux_SL.slice(0);
  imat aux_SL_Theta0  = aux_SL.slice(1);
  cube aux_Theta0_inv(N, N, M);
  cube aux_struc(N, N, M);
  
  if ( debug ) Rcout << " start!" << endl;
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
    
    if ( debug ) Rcout << " aux_Theta0" << endl;
    List Theta0SL_m         = sample_Theta0_Hou24_heterosk1_s4( aux_Theta0.slice(m), aux_B.slice(m), aux_SL_Theta0.col(m), aux_SLlpr.slice(1).col(m), shocks_m, aux_sigma_m, prior, VTheta0 );
    aux_Theta0.slice(m)     = as<mat>(Theta0SL_m["aux_Theta0"]);
    aux_Theta0_inv.slice(m) = inv(aux_Theta0.slice(m));
    aux_SL_Theta0.col(m)    = as<ivec>(Theta0SL_m["aux_SL"]);
    aux_SLlpr.slice(1).col(m) = as<vec>(Theta0SL_m["aux_SLlpr"]);
    
    if ( debug ) Rcout << " aux_B" << endl;
    List BSL_m              = sample_B_heterosk1_s4( aux_B.slice(m), aux_SL_B.col(m), aux_Theta0.slice(m), aux_Theta0_inv.slice(m), aux_SLlpr.slice(0).col(m), shocks_m, aux_sigma_m, prior_precision.col(m), prior, VB );
    aux_B.slice(m)          = as<mat>(BSL_m["aux_B"]);
    aux_SL_B.col(m)         = as<ivec>(BSL_m["aux_SL"]);
    aux_SLlpr.slice(0).col(m) = as<vec>(BSL_m["aux_SLlpr"]);
    
    aux_struc.slice(m)      = aux_Theta0_inv.slice(m) * aux_B.slice(m);
    
  }
  aux_SL.slice(0)       = aux_SL_B;
  aux_SL.slice(1)       = aux_SL_Theta0;
  
  return List::create(
    _["aux_B"]          = aux_B,
    _["aux_Theta0"]     = aux_Theta0,
    _["aux_struc"]      = aux_struc,
    _["aux_SL"]         = aux_SL,
    _["aux_SLlpr"]      = aux_SLlpr
  );
} // END sample_BTheta0_tvi

