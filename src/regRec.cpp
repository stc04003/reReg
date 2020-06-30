#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export]]
arma::vec scRate1(const arma::vec& T,
		  const arma::vec& Y,
		  const arma::vec& W,
		  const arma::vec& T0) {
  int n = Y.n_elem;
  int m = T0.n_elem;
  arma::vec out(m);
  out.zeros();
  arma::vec de(n);
  de.zeros();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ((T[i] <= Y[j]) && (T[i] >= T[j])) {
	de(i) += W[j];
      }
    }
  }
  for (int k = 0; k < m; k++) {
    for (int i = 0; i < n; i++) {
      if (T[i] >= T0[k] && de(i) > 0) {
        out[k] += W(i) / de(i);
      }
    }
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
arma::rowvec sc1Log3(const arma::vec& a,
		     const arma::mat& X,
		     const arma::vec& T,
		     const arma::vec& Y) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::vec texa = log(T) + X * a;
  arma::vec yexa = log(Y) + X * a;
  arma::rowvec out(p);
  out.zeros();
  arma::rowvec nu(p);
  for (int i = 0; i < n; i++) {
    nu.zeros();
    double de = 0;
    for (int j = 0; j < n; j++) {
      if ((texa[i] <= yexa[j]) & (texa[i] >= texa[j])) {
        nu += X.row(j);
        de += 1;
      }
    }
    out += X.row(i) - nu / de; 
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
arma::rowvec sc22(const arma::vec& b,
		  const arma::vec& R,
		  const arma::mat& X,
		  const arma::vec& W) {
  return (W % (R - exp(X * b))).t() * X / W.n_elem;
}

//' @noRd
// [[Rcpp::export]]
arma::rowvec sc1Gehan2(const arma::vec& a,
		       const arma::mat& X,
		       const arma::vec& T,
		       const arma::vec& Y) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::rowvec out(p);
  out.zeros();
  arma::vec texa = log(T) + X * a;
  arma::vec yexa = log(Y) + X * a;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ((texa[i] <= yexa[j]) & (texa[i] >= texa[j])) {
        out += X.row(i) - X.row(j);
      }
    }
  }
  return out;
}
