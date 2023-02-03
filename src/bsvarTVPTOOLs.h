
#ifndef _BSVARTVPTOOLS_H_
#define _BSVARTVPTOOLS_H_

#include <RcppArmadillo.h>


arma::field<arma::cube> bsvarTVPs_ir (
    arma::field<arma::cube>&  posterior_B,        // (S)(N, N, M)
    arma::cube&               posterior_A,        // (N, K, S)
    const int                 horizon,
    const int                 p
);


#endif  // _BSVARTVPTOOLS_H_