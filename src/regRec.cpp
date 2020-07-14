#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Used for recurrent event processes
//' @noRd
// [[Rcpp::export]]
arma::vec reRate(const arma::vec& T,
		 const arma::vec& Y,
		 const arma::vec& W,
		 const arma::vec& T0) {
  int n = Y.n_elem;
  int m = T0.n_elem;
  arma::vec out(m, arma::fill::zeros);
  arma::vec de(n, arma::fill::zeros);
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
arma::rowvec reLog(const arma::vec& a,
		   const arma::mat& X,
		   const arma::vec& T,
		   const arma::vec& Y,
		   const arma::vec& W) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::vec texa = log(T) + X * a;
  arma::vec yexa = log(Y) + X * a;
  arma::rowvec out(p, arma::fill::zeros);
  arma::rowvec nu(p);
  for (int i = 0; i < n; i++) {
    nu.zeros();
    double de = 0;
    for (int j = 0; j < n; j++) {
      if ((texa[i] <= yexa[j]) & (texa[i] >= texa[j])) {
        nu += W(j) * X.row(j);
        de += W(j);
      }
    }
    out += W(i) * (X.row(i) - nu / de); 
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
arma::rowvec re2(const arma::vec& b,
		 const arma::vec& R,
		 const arma::mat& X,
		 const arma::vec& W) {
  return (W % (R - exp(X * b))).t() * X / W.n_elem;
}

//' @noRd
// [[Rcpp::export]]
arma::vec reGehan(const arma::vec& a,
		  const arma::mat& X,
		  const arma::vec& T,
		  const arma::vec& Y,
		  const arma::vec& W) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::vec out(p, arma::fill::zeros);
  arma::vec texa = log(T) + X * a;
  arma::vec yexa = log(Y) + X * a;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ((texa[i] <= yexa[j]) & (texa[i] >= texa[j])) {
        out += W(i) * W(j) * (X.row(i) - X.row(j));
      }
    }
  }
  return out;
}

//' @noRd
// [[Rcpp::export]]
arma::rowvec am1(const arma::vec& a,
		 const arma::vec& T,
		 const arma::vec& Y,
		 const arma::vec& W,
		 const arma::mat& X,
		 const arma::vec& m) {
  int nm = accu(m);
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec m2 = cumsum(m); 
  arma::mat Xi(nm, p, arma::fill::zeros);
  arma::vec Yi(nm, arma::fill::zeros);
  arma::vec Wi(nm, arma::fill::zeros);
  arma::vec T0 = log(Y) + X * a;
  int mn = m.n_elem;
  for (int i = 0; i < mn; i ++) {
    if (i == 0 && m(i) > 0) {
      Wi.subvec(0, m2(i) - 1).fill(W(i));
      Yi.subvec(0, m2(i) - 1).fill(Y(i));
      Xi.submat(0, 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
    }
    if (i > 0 && m(i) > 0) {
      Wi.subvec(m2(i - 1), m2(i) - 1).fill(W(i));
      Yi.subvec(m2(i - 1), m2(i) - 1).fill(Y(i));
      Xi.submat(m2(i - 1), 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
    }
  }
  arma::vec texa = log(T) + Xi * a;
  arma::vec yexa = log(Yi) + Xi * a;  
  arma::vec Lam(n, arma::fill::zeros);
  arma::vec de(nm, arma::fill::zeros);
  for (int i = 0; i < nm; i++) {
    for (int j = 0; j < nm; j++) {
      if ((texa[i] <= yexa[j]) && (texa[i] >= texa[j])) {
	de(i) += Wi[j];
      }
    }
  }
  for (int k = 0; k < n; k++) {
    for (int i = 0; i < nm; i++) {
      if (texa[i] >= T0[k] && de(i) > 0) {
        Lam[k] += Wi(i) / de(i);
      }
    }
  }
  Lam = exp(-Lam); 
  arma::vec R = m / Lam;
  int Rn = R.n_elem;
  for (int i = 0; i < Rn; i ++) {
    if (R[i] > 100000) R[i] = (m[i] + 0.01) / (Lam[i] + 0.01);
  }
  // return (W % (R - mean(W % R))).t() * X / n;
  return (W % (mean(W) * R - mean(W % R))).t() * X / n;
  // return (W % (R - mean(R))).t() * X / n;
}

// Used for terminal events

//' @noRd
// [[Rcpp::export]]
arma::vec temHaz(const arma::vec& a,
		 const arma::vec& b,
		 const arma::mat& X,
		 const arma::vec& Y,
		 const arma::vec& Z,
		 const arma::vec& D,
		 const arma::vec& W,
		 const arma::vec& T0) {
  int n = Y.n_elem;
  int m = T0.n_elem;
  arma::vec out(m, arma::fill::zeros);
  arma::vec de(n, arma::fill::zeros);
  arma::vec ebax = exp(X * (b - a));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if ((Y[i] <= Y[j])) {
	de(i) += Z(j) * W(j) * ebax(j);
      }
    }
  }
  for (int k = 0; k < m; k++) {
    for (int i = 0; i < n; i++) {
      if (Y[i] <= T0[k] && de(i) > 0) {
        out[k] += D(i) * W(i) / de(i);
      }
    }
  }
  return out;
}


//' @noRd
// [[Rcpp::export]]
arma::rowvec temScLog(const arma::vec& a,
		      const arma::vec& b,
		      const arma::mat& X,
		      const arma::vec& Y,
		      const arma::vec& Z,
		      const arma::vec& D,
		      const arma::vec& W) {
  int n = Y.n_elem;
  arma::vec yexa = Y % exp(X * a);
  arma::vec ebax = exp(X * (b - a)); 
  arma::mat Iij = arma::conv_to<arma::mat>::from(repmat(yexa, 1, n) <= repmat(yexa, 1, n).t());
  arma::mat nu = Iij * (X % repmat(ebax % Z % W, 1, X.n_cols));
  arma::mat de = Iij * (ebax % Z % W);
  arma::mat tmp = nu / repmat(de, 1, nu.n_cols);
  tmp.replace(arma::datum::nan, 0);
  arma::mat D2 = repmat(D % W, 1, X.n_cols);
  arma::mat D3 = repmat(yexa % D % W, 1, X.n_cols);
  return join_rows(sum(X % D2, 0) - sum(tmp % D2, 0),
                   sum(X % D3, 0) - sum(tmp % D3, 0)) / n;
}

// [[Rcpp::export]]
Rcpp::NumericVector temScGehan(const arma::vec& a,
			       const arma::vec& b,
			       const arma::mat& X,
			       const arma::vec& Y,
			       const arma::vec& Z,
			       const arma::vec& D,
			       const arma::vec& W) {
  int n = Y.n_elem;
  int p = a.n_elem;
  Rcpp::NumericVector out(2 * p);
  arma::mat ebax = repmat(D, 1, n) % repmat(exp(X * (b - a)) % Z, 1, n).t();
  arma::vec yexa = Y % exp(X * a);
  arma::mat Iij = arma::conv_to<arma::mat>::from(repmat(yexa, 1, n) <= repmat(yexa, 1, n).t());
  for (int i = 0; i < p; i++) {
    arma::mat xdif=arma::conv_to<arma::mat>::from(repmat(W % X.col(i), 1, n) -
						  repmat(W % X.col(i), 1, n).t()); 
    out[i] = accu(xdif % ebax % Iij);
    out[i + p] = accu(repmat(yexa, 1, n) % xdif % ebax % Iij);
  }
  return out / n / n;
}

// [[Rcpp::export]]
arma::rowvec temLog(const arma::vec& a,
		    const arma::vec& b,
		    const arma::mat& X,
		    const arma::vec& Y,
		    const arma::vec& Z,
		    const arma::vec& D,
		    const arma::vec& W) {
  int n = Y.n_elem;
  arma::vec yexa = Y % exp(X * a);
  arma::vec ebax = exp(X * (b - a)); 
  arma::mat Iij = arma::conv_to<arma::mat>::from(repmat(yexa, 1, n) <= repmat(yexa, 1, n).t());
  // arma::mat nu = Iij * (X % repmat(ebax % Z % W, 1, X.n_cols));
  // arma::mat de = Iij * (ebax % Z % W);
  arma::mat nu = Iij * (X % repmat(ebax % Z, 1, X.n_cols));
  arma::mat de = Iij * (ebax % Z);
  arma::mat tmp = nu / repmat(de, 1, nu.n_cols);
  tmp.replace(arma::datum::nan, 0);
  arma::mat D2 = repmat(D % W, 1, X.n_cols);
  return (sum(X % D2, 0) - sum(tmp % D2, 0)) / n;
}

// [[Rcpp::export]]
Rcpp::NumericVector temGehan(const arma::vec& a,
			     const arma::vec& b,
			     const arma::mat& X,
			     const arma::vec& Y,
			     const arma::vec& Z,
			     const arma::vec& D,
			     const arma::vec& W) {
  int n = Y.n_elem;
  int p = a.n_elem;
  Rcpp::NumericVector out(p);
  arma::mat ebax = repmat(D, 1, n) % repmat(exp(X * (b - a)) % Z, 1, n).t();
  arma::vec yexa = Y % exp(X * a);
  arma::mat Iij = arma::conv_to<arma::mat>::from(repmat(yexa, 1, n) <= repmat(yexa, 1, n).t());
  for (int i = 0; i < p; i++) {
    arma::mat xdif=arma::conv_to<arma::mat>::from(repmat(W % X.col(i), 1, n) -
						  repmat(W % X.col(i), 1, n).t()); 
    out[i] = accu(xdif % ebax % Iij);
  }
  return out / n / n;
}

// Resampling 
