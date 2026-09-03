
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include "RcppTN.h"
#include "sample_ABTheta0hyper.h"

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
    vec logprob(10);
    for (int r = 0; r < 10; r++) {
      logprob(r) = std::log(pr_s(r)) - 0.5 * std::log(sigma_s(r))
        - 0.5 * std::pow(datanorm(j) - alpha_s(r), 2) / sigma_s(r);
    }
    // sigma_s contains variances; rescale before exponentiating to avoid underflow.
    logprob -= logprob.max();
    mixprob(first_index) = std::exp(logprob(0));
    for (int r = 1; r < 10; r++) {
      mixprob(first_index+r) = mixprob(first_index+r-1) + std::exp(logprob(r));
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
  
  if ( 
      !(std::isfinite(lambda) && std::isfinite(chi) && std::isfinite(psi)) ||
       (chi <  0. || psi < 0)      || 
       (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.) 
    ) {
    throw std::logic_error("do_rgig1 invalid input.");
  }
  
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



// All four SV samplers use g_t = omega_{s_t} h_t, h_0 = 0, and
// h_t = rho*h_{t-1} + v_t, v_t ~ N(0,1). The external h is always non-centred.
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double sample_sv_rho(
    const arma::rowvec& h
) {
  const int T = h.n_elem;
  if (T < 2) return R::runif(0, 1);
  const rowvec lag = h.cols(0, T-2);
  const double precision = dot(lag, lag);
  
  if (precision == 0) return R::runif(0, 1);
  return RcppTN::rtn1(dot(lag, h.cols(1, T-1))/precision,
                     1/std::sqrt(precision), 0, 1);
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::rowvec sample_sv_omega_nc(
    arma::rowvec omega,
    const arma::rowvec& h,
    const arma::rowvec& demeaned, 
    const arma::rowvec& precision,
    const arma::urowvec& state, 
    const double sigma2_omega
) {
  for (uword m=0; m<omega.n_elem; m++) {
    double vo = 1/sigma2_omega;
    double mo = 0;
    for (uword t=0; t<h.n_elem; t++) {
      if (state(t) == m) {
        vo += h(t)*h(t)*precision(t);
        mo += h(t)*precision(t)*demeaned(t);
      }
    }
    omega(m) = R::rnorm(mo/vo, 1/std::sqrt(vo));
  }
  return omega;
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::rowvec sample_sv_omega_ce(
    arma::rowvec omega,
    const arma::rowvec& g,
    const arma::urowvec& state, 
    const double rho, 
    const double sigma2_omega
) {
  const int T = g.n_elem;
  for (uword m=0; m<omega.n_elem; m++) {
    // z contains this regime's g; rest contains h from all other regimes.
    // Thus H*h = (H*z)/omega_m + H*rest, including regime boundaries.
    rowvec z(T, fill::zeros), rest(T, fill::zeros);
    int Tm = 0;
    for (int t=0; t<T; t++) {
      if (state(t) == m) {
        z(t) = g(t);
        Tm++;
      } else {
        rest(t) = g(t)/omega(state(t));
      }
    }
    for (int t=T-1; t>0; t--) {
      z(t) -= rho*z(t-1);
      rest(t) -= rho*rest(t-1);
    }
    const double a = dot(z, z);
    const double b = -dot(z, rest);
    // p(omega_m | g, ...) is proportional to
    // |omega_m|^{-Tm} exp[-(a/omega_m^2 + omega_m^2/sigma2_omega)/2 + b/omega_m].
    // A signed sqrt(GIG) proposal accounts for every term except b/omega_m.
    // For M=1 (or rho=0), b=0 and this is an exact Gibbs draw.
    const double variance = do_rgig1(0.5*(1-Tm), a, 1/sigma2_omega);
    const double proposal = (R::runif(0, 1) < 0.5 ? -1 : 1)*std::sqrt(variance);
    if (b == 0 || std::log(R::runif(0, 1)) < b*(1/proposal - 1/omega(m))) {
      omega(m) = proposal;
    }
  }
  return omega;
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::rowvec sample_sv_latent(
    const arma::rowvec& omega_T, 
    const arma::rowvec& demeaned,
    const arma::rowvec& precision, 
    const double rho, 
    const bool centred
) {
  const int T = demeaned.n_elem;
  vec diagonal(T, fill::value(1 + rho*rho));
  diagonal(T-1) = 1;
  if (!centred) {
    diagonal += trans(square(omega_T) % precision);
    return trans(precision_sampler_ar1(diagonal, -rho,
                 trans(omega_T % precision % demeaned)));
  }

  // Q_g = D^{-1} H'H D^{-1} + diag(precision), D = diag(omega_T).
  // At a regime boundary the off-diagonal is -rho/(omega_t*omega_{t-1}).
  diagonal = diagonal/trans(square(omega_T)) + precision.t();
  vec chol_diag(T), chol_offdiag(T, fill::zeros);
  chol_diag(0) = std::sqrt(diagonal(0));
  
  for (int t=1; t<T; t++) {
    chol_offdiag(t-1) = -rho/(omega_T(t)*omega_T(t-1)*chol_diag(t-1));
    chol_diag(t) = std::sqrt(diagonal(t) - std::pow(chol_offdiag(t-1), 2));
  }
  
  const vec location = trans(precision % demeaned);
  const vec forward = forward_algorithm(chol_diag, chol_offdiag, location);
  const vec epsilon(T, fill::randn);
  return trans(backward_algorithm(chol_diag, chol_offdiag, forward + epsilon));
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_asis(
    arma::rowvec& h, 
    double& rho, 
    arma::rowvec& omega, 
    arma::rowvec& sigma2v,
    double& sigma2_omega, 
    double& s, 
    arma::urowvec& S,
    const arma::mat& xi, 
    const arma::rowvec& u, 
    const Rcpp::List& prior,
    const bool centred
) {
  const int T = u.n_elem;
  const int M = omega.n_elem;
  
  const urowvec state = index_max(xi, 0);
  const rowvec U = log(square(u + 0.000000001));
  
  const rowvec alpha = {1.92677,1.34744,0.73504,0.02266,-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000};
  const rowvec variance = {0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342};
  
  const double a = prior["sv_a_"];
  const double prior_s = prior["sv_s_"];
  
  rowvec omega_T(T), demeaned(T), precision(T);
  
  for (int t=0; t<T; t++) omega_T(t) = omega(state(t));
  
  rowvec g = h % omega_T;
  S = trans(inverse_transform_sampling(find_mixture_indicator_cdf(trans(U-g)), T));
  
  for (int t=0; t<T; t++) {
    demeaned(t) = U(t) - alpha(S(t));
    precision(t) = 1/variance(S(t));
  }

  // Same normal-gamma hierarchy in both parameterisations:
  // omega_m | sigma2_omega ~ N(0,sigma2_omega), sigma2_omega | s ~ Gamma(a,scale=s).
  s = (prior_s + 2*sigma2_omega)/R::rchisq(3 + 2*a);
  sigma2_omega = do_rgig1(a - 0.5*M, accu(square(omega)), 2/s);
  rho = sample_sv_rho(h);

  if (centred) {
    
    // C -> NC: update omega and g in C, then hold h fixed for the NC update.
    omega = sample_sv_omega_ce(omega, g, state, rho, sigma2_omega);
    for (int t=0; t<T; t++) omega_T(t) = omega(state(t));
    g = sample_sv_latent(omega_T, demeaned, precision, rho, true);
    h = g/omega_T;
    omega = sample_sv_omega_nc(omega, h, demeaned, precision, state, sigma2_omega);
    
  } else {
    
    // NC -> C: update omega and h in NC, then hold g fixed for the C update.
    omega = sample_sv_omega_nc(omega, h, demeaned, precision, state, sigma2_omega);
    for (int t=0; t<T; t++) omega_T(t) = omega(state(t));
    h = sample_sv_latent(omega_T, demeaned, precision, rho, false);
    g = h % omega_T;
    omega = sample_sv_omega_ce(omega, g, state, rho, sigma2_omega);
    for (int t=0; t<T; t++) omega_T(t) = omega(state(t));
    h = g/omega_T;
    
  }
  
  sigma2v = square(omega);
  rho = sample_sv_rho(h);
  
  return List::create(
    _["aux_h_n"] = h, 
    _["aux_rho_n"] = rho, 
    _["aux_omega_n"] = omega,
    _["aux_sigma2v_n"] = sigma2v, 
    _["aux_sigma2_omega_n"] = sigma2_omega,
    _["aux_s_n"] = s, 
    _["aux_S_n"] = S
  );
}


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
    const Rcpp::List&     prior
) {
  rowvec omega(1, fill::value(aux_omega_n));
  rowvec sigma2v(1, fill::value(aux_sigma2v_n));
  const mat xi(1, u.n_elem, fill::ones);
  
  List out = svar_asis(aux_h_n, aux_rho_n, omega, sigma2v,
    aux_sigma2_omega_n, aux_s_n, aux_S_n, xi, u, prior, false);
  
  aux_omega_n = omega(0);
  aux_sigma2v_n = sigma2v(0);
  
  out["aux_omega_n"] = aux_omega_n;
  out["aux_sigma2v_n"] = aux_sigma2v_n;
  return out;
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
    const Rcpp::List&     prior
) {
  return svar_asis(aux_h_n, aux_rho_n, aux_omega_n, aux_sigma2v_n,
    aux_sigma2_omega_n, aux_s_n, aux_S_n, aux_xi, u, prior, false);
} // END sv_nc1




/*______________________function svar_ce1______________________
 * This function has been copied from the R package bsvars on 2025-05-05 and modified subsequently.*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_ce1 (
    arma::rowvec&       aux_h_tilde,            // 1xT
    double&             aux_rho_n,
    double&             aux_omega_n,
    double&             aux_sigma2v_n,
    double&             aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&             aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&      aux_S_n,            // 1xT
    const arma::rowvec& u,                  // 1xT
    const Rcpp::List&   prior
) {
  rowvec omega(1, fill::value(aux_omega_n));
  rowvec sigma2v(1, fill::value(aux_sigma2v_n));
  const mat xi(1, u.n_elem, fill::ones);
  
  List out = svar_asis(aux_h_tilde, aux_rho_n, omega, sigma2v,
    aux_sigma2_omega_n, aux_s_n, aux_S_n, xi, u, prior, true);
  
  aux_omega_n = omega(0);
  aux_sigma2v_n = sigma2v(0);
  
  out["aux_omega_n"] = aux_omega_n;
  out["aux_sigma2v_n"] = aux_sigma2v_n;
  return out;
} // END svar_ce1



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_ce1_mss (
    arma::rowvec&       aux_h_tilde,            // 1xT
    double&             aux_rho_n,
    arma::rowvec&       aux_omega_n,        // 1xM nth equation regime-dependent omegas
    arma::rowvec&       aux_sigma2v_n,      // 1xM nth equation regime-dependent omegas^2
    double&             aux_sigma2_omega_n, // omega prior hyper-parameter 
    double&             aux_s_n,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::urowvec&      aux_S_n,            // 1xT
    const arma::mat&    aux_xi,             // MxT
    const arma::rowvec& u,                  // 1xT
    const Rcpp::List&   prior
) {
  return svar_asis(aux_h_tilde, aux_rho_n, aux_omega_n, aux_sigma2v_n,
    aux_sigma2_omega_n, aux_s_n, aux_S_n, aux_xi, u, prior, true);
} // END svar_ce1_mss



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::mat count_regime_transitions (
    const arma::mat& xi
) {
  const int M = xi.n_rows;
  const int T = xi.n_cols;
  
  mat count(M, M, fill::zeros);
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




// Filter in log space so rare observations do not flatten regime likelihoods.
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat filtering_log_density(
    const arma::mat& log_density,
    const arma::mat& transition,
    const arma::vec& pi_0
) {
  mat filtered(log_density.n_rows, log_density.n_cols);
  vec previous = pi_0;
  for (uword t=0; t<log_density.n_cols; t++) {
    vec log_prob = log_density.col(t) + log(transition.t() * previous);
    vec prob = exp(log_prob - log_prob.max());
    filtered.col(t) = prob / accu(prob);
    previous = filtered.col(t);
  }
  return filtered;
}


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
  
  mat log_density(M, T);
  for (int m=0; m<M; m++) {
    log_density.row(m) = -0.5 * sum(square(Z.slice(m)), 0)
      - 0.5 * N * log(2*M_PI);
  }
  return filtering_log_density(log_density, aux_PR_TR, pi_0);
} // END filtering



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat filtering_studentt (
    const arma::cube& Z,                  // NxTxM state-specific standardised residuals
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  pi_0,               // Mx1
    const arma::mat&  aux_df              // NxM
) {
  
  // filtered probabilities for a model with MS structural matrix of SVAR-MSS-SV model
  
  const int   T = Z.n_cols;
  const int   N = Z.n_rows;
  const int   M = aux_PR_TR.n_rows;
  
  mat log_density(M, T, fill::zeros);
  for (int m=0; m<M; m++) {
    for (int n=0; n<N; n++) {
      const double df = aux_df(n,m);
      // Unit-variance Student-t: lambda has IG shape df/2 and scale (df-2)/2.
      log_density.row(m) += lgamma(0.5*(df+1)) - lgamma(0.5*df)
        - 0.5*log(M_PI*(df-2))
        - 0.5*(df+1)*log(1 + square(Z.slice(m).row(n))/(df-2));
    }
  }
  return filtering_log_density(log_density, aux_PR_TR, pi_0);
} // END filtering_studentt





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
arma::mat sample_markov_filtered(
    const arma::mat& filtered,
    const arma::mat& aux_xi,
    const arma::mat& transition,
    const bool finiteM
) {
  const int T = filtered.n_cols;
  const int M = filtered.n_rows;
  const mat states = eye(M, M);
  mat proposal(M, T);
  
  for (int attempt=0; attempt<(finiteM ? 10 : 1); attempt++) {
    
    int draw = csample_num1(wrap(seq_len(M)), wrap(filtered.col(T-1))) - 1;
    proposal.col(T-1) = states.col(draw);
    
    for (int t=T-2; t>=0; t--) {
      // Condition on the newly drawn next state, not the previous MCMC path.
      vec prob = filtered.col(t) % transition.col(draw);
      draw = csample_num1(wrap(seq_len(M)), wrap(prob)) - 1;
      proposal.col(t) = states.col(draw);
    }
    
    if (!finiteM) return proposal;
    
    mat transitions = count_regime_transitions(proposal);
    if (min(transitions.diag()) >= 2) return proposal;
  }
  return aux_xi;
}


// The SV callers supply log|det(B_m)| - sum_n log(sigma_{n,t,m}).
// Condition on the current scale-mixture variables, as do the other Gibbs blocks.
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Markov_process_sv(
    const arma::cube& Z,
    const arma::mat& log_jacobian,
    const arma::mat& aux_lambda,
    const arma::mat& aux_df,
    const arma::uvec& studentt,
    const arma::mat& aux_xi,
    const arma::mat& aux_PR_TR,
    const arma::vec& aux_pi_0,
    const bool finiteM
) {
  mat log_density = log_jacobian;
  for (uword m=0; m<Z.n_slices; m++) {
    for (uword n=0; n<Z.n_rows; n++) {
      if (studentt(n)) {
        const double shape = 0.5*aux_df(n,m);
        const double scale = 0.5*(aux_df(n,m)-2);
        log_density.row(m) += -0.5*square(Z.slice(m).row(n))/aux_lambda.row(n)
          + shape*log(scale) - lgamma(shape)
          - (shape+1)*log(aux_lambda.row(n)) - scale/aux_lambda.row(n);
      } else {
        log_density.row(m) -= 0.5*square(Z.slice(m).row(n));
      }
    }
  }
  // Normal density constants and log(lambda)/2 cancel across regimes.
  return sample_markov_filtered(
    filtering_log_density(log_density, aux_PR_TR, aux_pi_0), aux_xi, aux_PR_TR, finiteM);
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Markov_process (
    const arma::cube& Z,                  // NxTxM
    arma::mat         aux_xi,             // MxT
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  aux_pi_0,           // Mx1
    const bool        finiteM = true
) {
  return sample_markov_filtered(filtering(Z, aux_PR_TR, aux_pi_0), aux_xi, aux_PR_TR, finiteM);
} // END sample_Markov_process






// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Markov_process_studentt (
    const arma::cube& Z,                  // NxTxM
    arma::mat         aux_xi,             // MxT
    const arma::mat&  aux_PR_TR,          // MxM
    const arma::vec&  aux_pi_0,           // Mx1
    const arma::mat&  aux_df,             // NxM
    const bool        finiteM = true
) {
  return sample_markov_filtered(filtering_studentt(Z, aux_PR_TR, aux_pi_0, aux_df), aux_xi, aux_PR_TR, finiteM);
} // END sample_Markov_process_studentt







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
    vec prob_xi1          = aux_pi_0 % (aux_PR_TR * aux_xi.col(0));
    prob_xi1             /= sum(prob_xi1);
    int S0_draw           = csample_num1(wrap(seq_len(M)), wrap(prob_xi1));
    mat transitions       = count_regime_transitions(aux_xi);
    transitions(S0_draw-1, aux_xi.col(0).index_max())++;
    mat posterior_alpha   = transitions + prior_PR_TR;
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = rDirichlet1(posterior_alpha.row(m));
    }
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
