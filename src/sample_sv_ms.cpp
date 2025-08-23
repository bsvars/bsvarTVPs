
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include "RcppTN.h"
#include "sample_ABhyper.h"

using namespace Rcpp;
using namespace arma;



//---------------------------------------------------------------------------------------------------
// a transformed sample implementation taken from Rcpp Gallery:
// https://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/
// fixed to one draw, sampling without replacement, and changed output type to int
// IMPORTANT: always #include <RcppArmadilloExtensions/sample.h>
//---------------------------------------------------------------------------------------------------
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
int csample_num1 (
    Rcpp::NumericVector x,
    Rcpp::NumericVector prob = NumericVector::create()
) {
  bool replace = false;
  NumericVector ret = Rcpp::RcppArmadillo::sample(x, 1, replace, prob);
  int out           = ret(0);
  return out;
} // END csample_num1



/*______________________function find_mixture_indicator_cdf______________________*/
// utility function from file utils_latent_states.cc from the source code of package stochvol
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec find_mixture_indicator_cdf (
    const arma::vec& datanorm           // provide all that is conditionally normal
){
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  
  const int T = datanorm.n_elem;
  vec mixprob(10 * T);
  for (int j = 0; j < T; j++) {  // TODO slow (10*T calls to exp)!
    const int first_index = 10*j;
    mixprob(first_index) = std::exp(pr_s(0) - (datanorm(j) - alpha_s(0)) * (datanorm(j) - alpha_s(0)) / sigma_s(0) );
    for (int r = 1; r < 10; r++) {
      mixprob(first_index+r) = mixprob(first_index+r-1) + std::exp(pr_s(r) - (datanorm(j) - alpha_s(r)) * (datanorm(j) - alpha_s(r)) / sigma_s(r) );
    }
  }
  return mixprob;
} // END find_mixture_indicator_cdf




/*______________________function inverse_transform_sampling______________________*/
// utility function from file utils_latent_states.cc from the source code of package stochvol
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::uvec inverse_transform_sampling (
    const arma::vec&  mixprob,
    const int         T
) {
  uvec r(T);
  for (int j = 0; j < T; j++) {
    int index = (10-1)/2;  // start searching in the middle
    const double unnorm_cdf_value = R::unif_rand()*mixprob[9 + 10*j];  // current (non-normalized) value
    bool larger = false;  // indicates that we already went up
    bool smaller = false; // indicates that we already went down
    while(true) {
      if (unnorm_cdf_value > mixprob[index +  10*j]) {
        index++;
        if (smaller) {
          break;
        } else {
          larger = true;
        }
      } else if (larger || index == 0) {
        break;
      } else {
        index--;
        smaller = true;
      }
    }
    r[j] = index;
  }
  return r;
} // END inverse_transform_sampling




/*______________________function do_rgig1______________________*/
// utility function copied from package factorstochvol
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double do_rgig1(
    double lambda, 
    double chi, 
    double psi
) { 
  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
  return as<double>(fun(1, lambda, chi, psi));
} // END do_rgig1





/*______________________function cholesky_tridiagonal______________________*/
// utility function from file precision_sampler.cpp
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List cholesky_tridiagonal(
    const arma::vec&    omega_diag,
    const double&       omega_offdiag
) {
  const int T = omega_diag.n_elem - 1;
  vec chol_diag(T+1);
  vec chol_offdiag(T+1);
  chol_diag[0] = sqrt(omega_diag[0]);
  for (int j = 1; j < T+1; j++) {
    chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
    chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
  }
  return List::create(_["chol_diag"]=chol_diag, _["chol_offdiag"]=chol_offdiag);
} // END cholesky_tridiagonal



/*______________________function forward_algorithm______________________*/
// utility function from file precision_sampler.cpp
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector
) {
  const int T = chol_diag.n_elem - 1;
  vec htmp(T+1);
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < T+1; j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
  return htmp;
} // END forward_algorithm



/*______________________function backward_algorithm______________________*/
// utility function from file precision_sampler.cpp
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp
) {
  const int T = chol_diag.size() - 1;
  vec h(T+1);
  h[T] = htmp[T] / chol_diag[T];
  for (int j = T-1; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
  }
  return h;
} // END backward_algorithm

/*______________________function precision_sampler_ar1______________________*/
// utility function from file precision_sampler.cpp
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec precision_sampler_ar1(
    const arma::vec&     precision_diag,
    const double&        precision_offdiag,
    const arma::vec&     location
) {
  int T               = location.n_rows;
  vec  epsilon(T, fill::randn);
  List precision_chol = cholesky_tridiagonal(precision_diag, precision_offdiag);    // Cholesky decomposition using a dedicated technique
  vec  aa             = forward_algorithm(precision_chol["chol_diag"],              // this forward substitution can be performed outside of the loop
                                          precision_chol["chol_offdiag"],
                                                        location);
  vec draw_ssar1      = backward_algorithm(precision_chol["chol_diag"],
                                           precision_chol["chol_offdiag"],
                                                         aa + epsilon);     // this has to be done in the loop as function backward_algorithm requires covector to be a vector (not a matrix)
  return draw_ssar1;
} // END precision_sampler_ar1



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_nc1 (
    arma::rowvec    aux_h_n,            // 1xT
    double          aux_rho_n,
    double          aux_omega_n,
    double          aux_sigma2v_n,
    double          aux_sigma2_omega_n, // omega prior hyper-parameter 
    double          aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec   aux_S_n,            // 1xT
    const arma::rowvec&   u,                  // 1xT
    const Rcpp::List&     prior,
    bool            sample_s_ = true
) {
  // sampler for the non-centred parameterisation of the SV process
  
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_cols;
  const rowvec  U = log(pow(u + ccc, 2));
  
  const double  prior_sv_a_ = prior["sv_a_"];
  const double  prior_sv_s_ = prior["sv_s_"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho_n;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = find_mixture_indicator_cdf(trans(U - aux_omega_n*aux_h_n));
  aux_S_n               = trans(inverse_transform_sampling(mixprob, T));
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = alpha_s(aux_S_n(t));
    sigma_S_inv.col(t)  = 1/sigma_s(aux_S_n(t));
  }
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_n               = (prior_sv_s_ + 2 * aux_sigma2_omega_n)/chi2rnd(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega_n    = do_rgig1( prior_sv_a_-0.5, pow(aux_omega_n,2), 2/aux_s_n );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  double    aux_rho_var = as_scalar(pow(hm1*hm1.t(), -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  double    upper_bound = pow(1-aux_sigma2_omega_n, 0.5);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-upper_bound,upper_bound);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_omega
  double    V_omega_inv = 1/( as_scalar(aux_h_n * diagmat(sigma_S_inv) * aux_h_n.t()) + pow(aux_sigma2_omega_n, -1) );
  double    omega_bar   = as_scalar(aux_h_n * diagmat(sigma_S_inv) * (U - alpha_S).t());
  double    omega_aux   = randn( distr_param(V_omega_inv*omega_bar, sqrt(V_omega_inv) ));
  
  // sample aux_h
  mat       V_h         = pow(omega_aux, 2) * diagmat(sigma_S_inv) + HH_rho;
  vec       h_bar       = omega_aux * diagmat(sigma_S_inv) * (U - alpha_S).t();
  rowvec    h_aux       = trans(precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar));
  
  // ASIS
  rowvec    aux_h_tilde = omega_aux * h_aux;
  double    hHHh        = as_scalar( aux_h_tilde * HH_rho * aux_h_tilde.t() );
  aux_sigma2v_n         = do_rgig1( -0.5*(T-1), hHHh, 1/aux_sigma2_omega_n );
  int       ss=1;
  if (R::runif(0,1)<0.5) ss *= -1;
  aux_omega_n           = ss * sqrt(aux_sigma2v_n);
  aux_h_n               = aux_h_tilde / aux_omega_n;
  
  // ASIS: resample aux_rho
  hm1                   = aux_h_n.cols(0,T-2);
  aux_rho_var           = as_scalar(pow(hm1*hm1.t(), -1));
  aux_rho_mean          = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  upper_bound           = pow(1-aux_sigma2_omega_n, 0.5);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-upper_bound,upper_bound);
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2v_n"]        = aux_sigma2v_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END sv_nc1




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_nc1_mss (
    arma::rowvec&         aux_h_n,            // 1xT
    double&               aux_rho_n,
    arma::rowvec&         aux_omega_n,        // 1xM nth equation regime-dependent omegas
    arma::rowvec&         aux_sigma2v_n,      // 1xM nth equation regime-dependent omegas^2
    double&               aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&               aux_s_n,            // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&        aux_S_n,            // 1xT
    const arma::mat&      aux_xi,             // MxT
    const arma::rowvec&   u,                  // 1xT
    const Rcpp::List&     prior,
    bool                  sample_s_ = true
) {
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_cols;
  const int     M = aux_xi.n_rows;
  const rowvec  U = log(square(u + ccc));
  const vec    Tm= sum(aux_xi, 1);
  
  rowvec        omega_aux(M);
  
  const double  prior_sv_a_ = prior["sv_a_"];
  const double  prior_sv_s_ = prior["sv_s_"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho_n;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  rowvec        omega_T(T);
  for (int t=0; t<T; t++) {
    omega_T(t)          = aux_omega_n(aux_xi.col(t).index_max());
  }
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = find_mixture_indicator_cdf(trans(U - aux_h_n % omega_T));
  aux_S_n               = trans(inverse_transform_sampling(mixprob, T));
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = alpha_s(aux_S_n(t));
    sigma_S_inv.col(t)  = 1/sigma_s(aux_S_n(t));
  }
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_n               = (prior_sv_s_ + 2 * aux_sigma2_omega_n)/R::rchisq(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega_n    = do_rgig1( prior_sv_a_ - 0.5 * M, accu(square(aux_omega_n)), 2/aux_s_n );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  double    aux_rho_var = as_scalar(pow(hm1*hm1.t(), -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  double    upper_bound = sqrt(1-aux_sigma2_omega_n);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5), -upper_bound, upper_bound);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_omega
  for (int m=0; m<M; m++) {
    rowvec  aux_h_n_m(Tm(m));
    rowvec  U_m(Tm(m));
    rowvec  alpha_S_m(Tm(m));
    rowvec  sigma_S_inv_m(Tm(m));
    
    int ii = 0;
    for (int t=0; t<T; t++) {
      if (aux_xi(m,t)==1) {
        aux_h_n_m(ii)   = aux_h_n(t);
        U_m(ii)         = U(t);
        alpha_S_m(ii)   = alpha_S(t);
        sigma_S_inv_m(ii) = sigma_S_inv(t);
        ii++;
      }
    }
    double    V_omega_inv = 1/( as_scalar(aux_h_n_m * diagmat(sigma_S_inv_m) * aux_h_n_m.t()) + pow(aux_sigma2_omega_n, -1) );
    double    omega_bar   = as_scalar(aux_h_n_m * diagmat(sigma_S_inv_m) * (U_m - alpha_S_m).t());
    omega_aux(m)          = R::rnorm(V_omega_inv*omega_bar, sqrt(V_omega_inv) );
  } // END m loop
  
  for (int t=0; t<T; t++) {
    omega_T(t)          = omega_aux(aux_xi.col(t).index_max());
  }
  
  // sample aux_h
  mat       V_h         = diagmat(square(omega_T)) * diagmat(sigma_S_inv) + HH_rho;
  vec       h_bar       = diagmat(omega_T) * diagmat(sigma_S_inv) * (U - alpha_S).t();
  rowvec    h_aux       = trans(precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar));
  
  // ASIS
  rowvec    aux_h_tilde = h_aux % omega_T;
  rowvec    aux_h_tilde_ar = trans(H_rho * aux_h_tilde.t());
  
  for (int m=0; m<M; m++) {
    double  hHHh = 0;
    for (int t=0; t<T; t++) {
      if (aux_xi(m,t)==1) {
        hHHh           += pow(aux_h_tilde(t), 2);
      }
    }
    double  sigma2_aux  = do_rgig1( -0.5*(Tm(m)-1), hHHh, 1/aux_sigma2_omega_n );
    int     ss = 1;
    if (R::runif(0,1)<0.5) ss *= -1;
    aux_omega_n(m)      = ss * sqrt(sigma2_aux);
    aux_sigma2v_n(m)    = sigma2_aux;
  } // END m loop
  
  for (int t=0; t<T; t++) {
    omega_T(t)          = aux_omega_n(aux_xi.col(t).index_max());
  }
  aux_h_n               = aux_h_tilde % pow(omega_T, -1);
  
  // resample aux_rho
  hm1                   = aux_h_n.cols(0,T-2);
  aux_rho_var           = as_scalar(pow(hm1*hm1.t(), -1));
  aux_rho_mean          = as_scalar(aux_rho_var * hm1*aux_h_n.cols(1,T-1).t());
  upper_bound           = sqrt(1-aux_sigma2_omega_n);
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5), -upper_bound, upper_bound);
  
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2v_n"]        = aux_sigma2v_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END sv_nc1




/*______________________function svar_ce1______________________
 * This function has been copied from the R package bsvars on 2025-05-05 and modified subsequently.*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_ce1 (
    arma::rowvec&       aux_h_n,            // 1xT
    double&             aux_rho_n,
    double&             aux_omega_n,
    double&             aux_sigma2v_n,
    double&             aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&             aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&      aux_S_n,            // 1xT
    const arma::rowvec& u,                  // 1xT
    const Rcpp::List&   prior,
    bool                sample_s_ = true
) {
  // sampler for the centred parameterisation of the SV process
  
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_cols;
  const rowvec  U = log(pow(u + ccc, 2));
  
  const double  prior_sv_a_ = prior["sv_a_"];
  const double  prior_sv_s_ = prior["sv_s_"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho_n;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = find_mixture_indicator_cdf(trans(U - aux_h_n));
  aux_S_n               = trans(inverse_transform_sampling(mixprob, T));
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = alpha_s(aux_S_n(t));
    sigma_S_inv.col(t)  = 1/sigma_s(aux_S_n(t));
  }
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_n               = (1 + 2 * aux_sigma2_omega_n) / chi2rnd(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega_n    = randg( distr_param(1 + 0.5 * prior_sv_a_, pow(pow(prior_sv_s_,-1) + pow(2 * aux_sigma2v_n,-1), -1)  ) );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  double    aux_rho_var = as_scalar(pow( hm1 * hm1.t() / aux_sigma2v_n, -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * (hm1 * aux_h_n.cols(1,T-1).t() / aux_sigma2v_n) );
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-1,1);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_sigma2v
  aux_sigma2v_n         = (aux_sigma2_omega_n + as_scalar(aux_h_n * HH_rho * aux_h_n.t())) / chi2rnd( 3 + T );
  aux_omega_n           = pow(aux_sigma2v_n, 0.5);
  
  // sample aux_h
  mat       V_h         = diagmat(sigma_S_inv) + (HH_rho / aux_sigma2v_n);
  vec       h_bar       = diagmat(sigma_S_inv) * (U - alpha_S).t();
  aux_h_n               = trans(precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar));
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2v_n"]        = aux_sigma2v_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END svar_ce1



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_ce1_mss (
    arma::rowvec&       aux_h_n,            // 1xT
    double&             aux_rho_n,
    arma::rowvec&       aux_omega_n,        // 1xM nth equation regime-dependent omegas
    arma::rowvec&       aux_sigma2v_n,      // 1xM nth equation regime-dependent omegas^2
    double&             aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&             aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&      aux_S_n,            // 1xT
    const arma::mat&    aux_xi,             // MxT
    const arma::rowvec& u,                  // 1xT
    const Rcpp::List&   prior,
    bool                sample_s_ = true
) {
  // sampler for the centred parameterisation of the SV process
  
  // fixed values for auxiliary mixture
  const NumericVector alpha_s = NumericVector::create(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000);
  const NumericVector sigma_s = NumericVector::create(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342);
  const NumericVector pr_s    = NumericVector::create(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115);
  const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_cols;
  const rowvec  U = log(pow(u + ccc, 2));
  const int     M = aux_xi.n_rows;
  const vec    Tm = sum(aux_xi, 1);
  
  const double  prior_sv_a_ = prior["sv_a_"];
  const double  prior_sv_s_ = prior["sv_s_"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho_n;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  rowvec        sigma2v_T_inv(T);
  for (int t=0; t<T; t++) {
    sigma2v_T_inv(t)          = pow( aux_sigma2v_n(aux_xi.col(t).index_max()), -1);
  }
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = find_mixture_indicator_cdf(trans(U - aux_h_n));
  aux_S_n               = trans(inverse_transform_sampling(mixprob, T));
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = alpha_s(aux_S_n(t));
    sigma_S_inv.col(t)  = 1/sigma_s(aux_S_n(t));
  }
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_n               = (prior_sv_s_ + 2 * aux_sigma2_omega_n)/R::rchisq(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega_n    = randg( distr_param(1 + 0.5 * prior_sv_a_, pow(pow(prior_sv_s_,-1) + accu(pow(2 * aux_sigma2v_n,-1)), -1)  ) );
  
  // sample aux_rho
  rowvec    hm1         = aux_h_n.cols(0,T-2);
  mat       sigma_v2_inv_diag = diagmat(sigma2v_T_inv.cols(0,T-2));
  double    aux_rho_var = as_scalar(pow( hm1 * sigma_v2_inv_diag * hm1.t(), -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * (hm1 * sigma_v2_inv_diag * aux_h_n.cols(1,T-1).t() ) );
  aux_rho_n             = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-1,1);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho_n;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_sigma2v
  for (int m=0; m<M; m++) {
    rowvec  aux_h_n_m(Tm(m));
    
    int ii = 0;
    for (int t=0; t<T; t++) {
      if (aux_xi(m,t)==1) {
        aux_h_n_m(ii)   = aux_h_n(t);
        ii++;
      }
    }
    
    aux_sigma2v_n(m)      = (aux_sigma2_omega_n + as_scalar(aux_h_n_m * HH_rho.submat(0,0,Tm(m)-1,Tm(m)-1) * aux_h_n_m.t())) / chi2rnd( 3 + Tm(m) );
  } // END m loop
  
  aux_omega_n        = pow(aux_sigma2v_n, 0.5);
  for (int t=0; t<T; t++) {
    sigma2v_T_inv(t)          = pow( aux_sigma2v_n(aux_xi.col(t).index_max()), -1);
  }
  
  // sample aux_h
  mat       V_h         = diagmat(sigma_S_inv) + (H_rho_new.t() * diagmat(sigma2v_T_inv) * H_rho_new);
  vec       h_bar       = diagmat(sigma_S_inv) * (U - alpha_S).t();
  vec aux_h_n_tmp       = precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar);
  
  if (aux_h_n_tmp.has_nan()) throw std::runtime_error("Error: aux_h_n_tmp contains missing observations, nan.");
  aux_h_n               = aux_h_n_tmp.t();
  
  return List::create(
    _["aux_h_n"]              = aux_h_n,
    _["aux_rho_n"]            = aux_rho_n,
    _["aux_omega_n"]          = aux_omega_n,
    _["aux_sigma2v_n"]        = aux_sigma2v_n,
    _["aux_sigma2_omega_n"]   = aux_sigma2_omega_n,
    _["aux_s_n"]              = aux_s_n,
    _["aux_S_n"]              = aux_S_n
  );
} // END svar_ce1_mss



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::mat count_regime_transitions (
    const arma::mat& xi
) {
  const int M = xi.n_rows;
  const int T = xi.n_cols;
  
  mat count(M, M);
  urowvec s   = index_max( xi, 0 );
  for (int t=1; t<T; t++) {
    count( s(t-1), s(t))++;
  }
  return count;
} // END count_regime_transitions



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::rowvec rDirichlet1 (
    const arma::rowvec&   alpha     // Kx1
) {
  const int K   = alpha.size();
  rowvec    draw(K);
  for (int k=0; k<K; k++) {
    draw(k)     = randg(distr_param(alpha(k), 1.0));
  }
  return draw/sum(draw);
} // END rDirichlet1




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat filtering (
    const arma::cube& Z,                  // NxTxM state-specific standardised residuals
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  pi_0                // Mx1
) {
  
  // filtered probabilities for a model with MS structural matrix of SVAR-MSS-SV model
  
  const int   T = Z.n_cols;
  const int   N = Z.n_rows;
  const int   M = aux_PR_TR.n_rows;
  
  mat         eta_t(M, T);
  mat         xi_t_t(M, T);
  
  // This loop evaluates mvnormal pdf at Z being multivariate standard normal distribution
  for (int m=0; m<M; m++) {
    rowvec log_d    = -0.5 * sum(square(Z.slice(m)), 0);
    log_d          += -0.5 * N * log(2*M_PI);
    NumericVector   exp_log_d   = wrap(exp(log_d));
    exp_log_d[exp_log_d==0]     = 1e-300;
    eta_t.row(m)    = as<rowvec>(exp_log_d);
  } // END m loop
  
  vec xi_tm1_tm1    = pi_0;
  
  for (int t=0; t<T; t++) {
    vec     num     = eta_t.col(t) % (aux_PR_TR.t() * xi_tm1_tm1);
    double  den     = sum(num);
    xi_t_t.col(t)   = num/den;
    xi_tm1_tm1      = xi_t_t.col(t);
  } // END t loop
  
  return xi_t_t;
} // END filtering



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat smoothing (
    const arma::mat&  filtered,           // MxT
    const arma::mat&  aux_PR_TR           // MxM
) {
  // NOT the same as for msh (but could be the same if you get rid of arg U in the other one)
  const int   T = filtered.n_cols;
  const int   M = aux_PR_TR.n_rows;
  
  mat   smoothed(M, T);
  smoothed.col(T-1)   = filtered.col(T-1);
  
  for (int t=T-2; t>=0; --t) {
    smoothed.col(t)   = (aux_PR_TR * (smoothed.col(t+1)/(aux_PR_TR.t() * filtered.col(t)) )) % filtered.col(t);
  } // END t loop
  
  return smoothed;
} // smoothing










// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Markov_process_mss (
    arma::mat         aux_xi,             // MxT
    const arma::mat&  E,                  // NxT
    const arma::cube& aux_B,              // NxNxM
    const arma::mat&  aux_sigma,          // NxM
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  aux_pi_0,           // Mx1
    const bool        finiteM = true
) {
  
  int minimum_regime_occurrences = 10;
  int max_iterations = 50;
  if ( finiteM ) {
    minimum_regime_occurrences = 10;
    max_iterations = 50;
  }
  
  const int   T   = E.n_cols;
  const int   N   = E.n_rows;
  const int   M   = aux_PR_TR.n_rows;
  mat aux_xi_tmp = aux_xi;
  mat aux_xi_out = aux_xi;
  
  cube  Z(N, T, M);
  for (int m=0; m<M; m++) {
    Z.slice(m)    = pow(aux_sigma, -1) % (aux_B.slice(m) * E);
  }
  
  mat filtered    = filtering(Z, aux_PR_TR, aux_pi_0);
  mat smoothed    = smoothing(filtered, aux_PR_TR);
  mat    aj       = eye(M, M);
  
  // Rcout << "aux_xi: " << aux_xi << endl;
  // Rcout << "aux_B: " << aux_B << endl;
  // Rcout << " smoothed.col(T-1): " << smoothed.col(T-1) << endl;
  
  mat xi(M, T);
  int draw        = csample_num1(wrap(seq_len(M)), wrap(smoothed.col(T-1)));
  aux_xi_tmp.col(T-1)     = aj.col(draw-1);
  
  if ( minimum_regime_occurrences==0 ) {
    for (int t=T-2; t>=0; --t) {
      vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
      draw          = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
      aux_xi_tmp.col(t)   = aj.col(draw-1);
    }
    aux_xi_out = aux_xi_tmp;
  } else {
    int regime_occurrences  = 1;
    int iterations  = 1;
    while ( (regime_occurrences<minimum_regime_occurrences) & (iterations<max_iterations) ) {
      for (int t=T-2; t>=0; --t) {
        vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
        draw          = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
        aux_xi_tmp.col(t)   = aj.col(draw-1);
      }
      mat transitions       = count_regime_transitions(aux_xi_tmp);
      regime_occurrences    = min(transitions.diag());
      iterations++;
    } // END while
    if ( iterations<max_iterations ) aux_xi_out = aux_xi_tmp;
  }
  
  return aux_xi_out;
} // END sample_Markov_process_mss





// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Markov_process_mssa (
    arma::mat         aux_xi,             // MxT
    const arma::cube& aux_B,              // NxNxM
    const arma::cube& aux_A,              // NxKxM
    const arma::mat&  Y,
    const arma::mat&  X,
    const arma::mat&  aux_sigma,          // NxM
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  aux_pi_0,           // Mx1
    const bool        finiteM = true
) {
  
  int minimum_regime_occurrences = 10;
  int max_iterations = 50;
  if ( finiteM ) {
    minimum_regime_occurrences = 10;
    max_iterations = 50;
  }
  
  const int   T   = Y.n_cols;
  const int   N   = Y.n_rows;
  const int   M   = aux_PR_TR.n_rows;
  mat aux_xi_tmp = aux_xi;
  mat aux_xi_out = aux_xi;
  
  cube  Z(N, T, M);
  for (int m=0; m<M; m++) {
    Z.slice(m)    = pow(aux_sigma, -1) % (aux_B.slice(m) * (Y - aux_A.slice(m) * X));
  }
  
  mat filtered    = filtering(Z, aux_PR_TR, aux_pi_0);
  mat smoothed    = smoothing(filtered, aux_PR_TR);
  mat    aj       = eye(M, M);
  
  mat xi(M, T);
  int draw        = csample_num1(wrap(seq_len(M)), wrap(smoothed.col(T-1)));
  aux_xi_tmp.col(T-1)     = aj.col(draw-1);
  
  if ( minimum_regime_occurrences==0 ) {
    for (int t=T-2; t>=0; --t) {
      vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
      draw          = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
      aux_xi_tmp.col(t)   = aj.col(draw-1);
    }
    aux_xi_out = aux_xi_tmp;
  } else {
    int regime_occurrences  = 1;
    int iterations  = 1;
    while ( (regime_occurrences<minimum_regime_occurrences) & (iterations<max_iterations) ) {
      for (int t=T-2; t>=0; --t) {
        vec xi_Tmj    = (aux_PR_TR * (aux_xi.col(t+1)/(aux_PR_TR.t() * filtered.col(t)))) % filtered.col(t);
        draw          = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
        aux_xi_tmp.col(t)   = aj.col(draw-1);
      }
      mat transitions       = count_regime_transitions(aux_xi_tmp);
      regime_occurrences    = min(transitions.diag());
      iterations++;
    } // END while
    if ( iterations<max_iterations ) aux_xi_out = aux_xi_tmp;
  }
  
  return aux_xi_out;
} // END sample_Markov_process_mssa





















// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_transition_probabilities (
    arma::mat           aux_PR_TR,    // MxM 
    arma::vec           aux_pi_0,     // Mx1
    const arma::mat&    aux_xi,       // MxT
    const Rcpp::List&   prior,         // a list of priors - original dimensions
    const bool          MSnotMIX = true
) {
  // the function changes the value of aux_PR_TR and aux_pi_0 by reference (filling it with a new draw)
  const int   M           = aux_PR_TR.n_rows;
  const mat   prior_PR_TR = as<mat>(prior["PR_TR"]);
  
  if ( MSnotMIX ) {
    mat transitions       = count_regime_transitions(aux_xi);
    mat posterior_alpha   = transitions + prior_PR_TR;
    
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = rDirichlet1(posterior_alpha.row(m));
    }
    vec prob_xi1          = aux_PR_TR *aux_xi.col(0);
    prob_xi1             /= sum(prob_xi1);
    int S0_draw           = csample_num1(wrap(seq_len(M)), wrap(prob_xi1));
    rowvec posterior_alpha_0(M, fill::value((1.0)));
    posterior_alpha_0(S0_draw-1)++;
    aux_pi_0              = trans(rDirichlet1(posterior_alpha_0));
  } else {
    rowvec occurrences    = trans(sum(aux_xi, 1));
    rowvec posterior_alpha= occurrences + prior_PR_TR.row(0);
    aux_pi_0              = trans(rDirichlet1(posterior_alpha));
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = aux_pi_0.t();
    }
  }
  
  return List::create(
    _["aux_PR_TR"]        = aux_PR_TR,
    _["aux_pi_0"]         = aux_pi_0
    );
} // END sample_transition_probabilities



