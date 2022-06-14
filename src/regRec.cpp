#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace Rcpp;

using cmp_par = std::pair<double, arma::uword>;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Used for recurrent event processes
//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec reRate(const arma::vec& T,
		 const arma::vec& Y,
		 const arma::vec& W,
		 const arma::vec& T0) {
  arma::uword const n = Y.n_elem;
  arma::uword const m = T0.n_elem;
  arma::vec out(m, arma::fill::zeros);
  arma::vec de(n, arma::fill::zeros); 
  arma::uvec const idx = arma::sort_index(T);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
	       return x.first <= y.first;
	     };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{}; 
  {
    auto const idx_i = idx[0];
    indices.emplace(Y[idx_i], idx_i);
    w_sum += W(idx_i);
    de(idx_i) = w_sum;
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(Y[idx_i], idx_i);
    if(Y[idx_i] > indices_head->first)
      while(indices_head->first < T[idx_i]){
        w_sum -= W(indices_head->second);
        ++indices_head;
      }
    else
      --indices_head;
    w_sum += W(idx_i);
    de(idx_i) = w_sum;
  }
  for (arma::uword k = 0; k < m; k++) {
    for (arma::uword i = 0; i < n; i++) {
      if (T[i] >= T0[k] && de(i) > 0) {
        out[k] += W(i) / de(i);
      }
    }
  }
  return out;
}


// arma::vec reRate(const arma::vec& T,
// 		 const arma::vec& Y,
// 		 const arma::vec& W,
// 		 const arma::vec& T0) {
//   int n = Y.n_elem;
//   int m = T0.n_elem;
//   arma::vec out(m, arma::fill::zeros);
//   arma::vec de(n, arma::fill::zeros);
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       if ((T[i] <= Y[j]) && (T[i] >= T[j])) {
// 	de(i) += W[j];
//       }
//     }
//   }
//   for (int k = 0; k < m; k++) {
//     for (int i = 0; i < n; i++) {
//       if (T[i] >= T0[k] && de(i) > 0) {
//         out[k] += W(i) / de(i);
//       }
//     }
//   }
//   return out;
// }

arma::mat matvec(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_cols; i++) {
    out.col(i) = x.col(i) % y;
  }
  return out;
}

// Log-rank estimating equation 
//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec reLog(const arma::vec& a,
		const arma::mat& X,
		const arma::vec& T,
		const arma::vec& Y,
		const arma::vec& W) {
  arma::uword const n = Y.n_elem;
  arma::uword const p = a.n_elem;
  arma::vec out(p, arma::fill::zeros);
  if(n < 1) return out;
  arma::vec texa = log(T) + X.t() * a;
  arma::vec yexa = log(Y) + X.t() * a;
  arma::uvec const idx = arma::sort_index(texa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
	       return x.first <= y.first;
	     };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};
  arma::vec x_col_sum(p, arma::fill::zeros);
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    x_col_sum += W(idx_i) * X.col(idx_i);
    w_sum += W(idx_i);
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first)
      while(indices_head->first < texa[idx_i]){
        x_col_sum -= W(indices_head->second) * X.col(indices_head->second);
        w_sum -= W(indices_head->second);
        ++indices_head;
      }
    else
      --indices_head;
    x_col_sum += W(idx_i) * X.col(idx_i);
    w_sum += W(idx_i);
    if (w_sum > 0) 
      out += W(idx_i) * (X.col(idx_i) - x_col_sum / w_sum);
  }
  return out;
}

// arma::rowvec reLog(const arma::vec& a,
// 		   const arma::mat& X,
// 		   const arma::vec& T,
// 		   const arma::vec& Y,
// 		   const arma::vec& W) {
//   int n = Y.n_elem;
//   int p = a.n_elem;
//   arma::vec texa = log(T) + X * a;
//   arma::vec yexa = log(Y) + X * a;
//   texa.replace(-arma::datum::inf, 0);
//   yexa.replace(-arma::datum::inf, 0);  
//   arma::rowvec out(p, arma::fill::zeros);
//   arma::mat XW = matvec(X, W);
//   for (int i = 0; i < n; i++) {
//     arma::uvec w = find((texa - texa[i]) % (yexa - texa[i]) <= 0);
//     out += W(i) * (X.row(i) - sum(XW.rows(w), 0) / sum(W(w)));
//     // out += W(i) * (X.row(i) - sum(X.rows(w), 0) / w.n_elem);
//   }
//   return out;
// }

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::rowvec re2(const arma::vec& b,
		 const arma::vec& R,
		 const arma::mat& X,
		 const arma::vec& W) {
  return (W % (R - exp(X * b))).t() * X / W.n_elem;
}

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec reGehan(const arma::vec& a,
		  const arma::mat& X,
		  const arma::vec& T,
		  const arma::vec& Y,
		  const arma::vec& W) {
  arma::uword const n = Y.n_elem;
  arma::uword const p = a.n_elem;
  arma::vec out(p, arma::fill::zeros);
  if(n < 1) return out;
  arma::vec texa = log(T) + X.t() * a;
  arma::vec yexa = log(Y) + X.t() * a;
  arma::uvec const idx = arma::sort_index(texa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};
  arma::vec x_col_sum(p, arma::fill::zeros);
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    x_col_sum += W(idx_i) * X.col(idx_i);
    w_sum += W(idx_i);
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first)
      while(indices_head->first < texa[idx_i]){
        x_col_sum -= W(indices_head->second) * X.col(indices_head->second);
        w_sum -= W(indices_head->second);
        ++indices_head;
      }
    else
      --indices_head;
    x_col_sum += W(idx_i) * X.col(idx_i);
    w_sum += W(idx_i);
    out += W(idx_i) * (w_sum * X.col(idx_i) - x_col_sum);
  }
  return out;
}

// arma::rowvec reGehan(const arma::vec& a,
// 		  const arma::mat& X,
// 		  const arma::vec& T,
// 		  const arma::vec& Y,
// 		  const arma::vec& W) {
//   int n = Y.n_elem;
//   int p = a.n_elem;
//   arma::rowvec out(p, arma::fill::zeros);
//   arma::vec texa = log(T) + X * a;
//   arma::vec yexa = log(Y) + X * a;
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       if ((texa[i] <= yexa[j]) & (texa[i] >= texa[j])) {
//         out += W(i) * W(j) * (X.row(i) - X.row(j));
//       }
//     }
//   }
//   return out;
// }

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::rowvec reGehan_s(const arma::vec& a,
		       const arma::mat& X,
		       const arma::vec& T,
		       const arma::vec& Y,
		       const arma::vec& W,
		       double nc) {
  int n = Y.n_elem;
  int p = a.n_elem;
  arma::vec texa = log(T) + X * a;
  arma::vec yexa = log(Y) + X * a;
  arma::rowvec out(p, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      arma::rowvec xdif = (X.row(i) - X.row(j));
      arma::vec rij = xdif * xdif.t();
      rij[0] = sqrt(rij[0] / nc);
      double H = 0;
      if (rij[0] != 0) {
	H = arma::normcdf((yexa[j] - texa[i]) / rij(0)) -
	  arma::normcdf((texa[j] - texa[i]) / rij(0));
      }
      out += W(i) * W(j) * xdif * H; 
    }
  }
  return out;
}

//' @noRd
// [[Rcpp::export(rng = false)]]

arma::rowvec am1(const arma::vec& a,
		 const arma::vec& T,
		 const arma::vec& Y,
		 const arma::vec& W,
		 const arma::mat& X,
		 const arma::vec& m) {
  arma::uword const nm = accu(m);
  arma::uword const n = X.n_rows;
  arma::uword const p = X.n_cols;
  arma::vec m2 = cumsum(m); 
  arma::mat Xi(nm, p, arma::fill::zeros);
  arma::vec Yi(nm, arma::fill::zeros);
  // arma::vec Wi(nm, arma::fill::zeros);
  arma::vec T0 = log(Y) + X * a;
  arma::uword const mn = m.n_elem;
  for (arma::uword i = 0; i < mn; i++) {
    if (i == 0 && m(i) > 0) {
      Yi.subvec(0, m2(i) - 1).fill(Y(i));
      Xi.submat(0, 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
    }
    if (i > 0 && m(i) > 0) {
      Yi.subvec(m2(i - 1), m2(i) - 1).fill(Y(i));
      Xi.submat(m2(i - 1), 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
    }
  }
  arma::vec texa = log(T) + Xi * a;
  arma::vec yexa = log(Yi) + Xi * a;
  arma::vec Lam(n, arma::fill::zeros);
  arma::vec de(nm, arma::fill::zeros);
  arma::uvec const idx = arma::sort_index(texa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
	       return x.first <= y.first;
	     };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{}; 
  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    w_sum += 1; // W(idx_i);
    de(idx_i) = w_sum;
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < nm; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first)
      while(indices_head->first < texa[idx_i]){
        w_sum -= 1; // W(indices_head->second);
        ++indices_head;
      }
    else
      --indices_head;
    w_sum += 1; //W(idx_i);
    de(idx_i) = w_sum;
  }
  for (arma::uword k = 0; k < n; k++) {
    for (arma::uword i = 0; i < nm; i++) {
      if (texa[i] >= T0[k] && de(i) > 0) {
        Lam[k] += 1 / de(i); // W(i) / de(i);
      }
    }
  }
  Lam = exp(-Lam); 
  arma::vec R = m / Lam;  
  return ((W % R - mean(W % R))).t() * X / n;
}

// arma::rowvec am1(const arma::vec& a,
// 		 const arma::vec& T,
// 		 const arma::vec& Y,
// 		 const arma::vec& W,
// 		 const arma::mat& X,
// 		 const arma::vec& m) {
//   int nm = accu(m);
//   int n = X.n_rows;
//   int p = X.n_cols;
//   arma::vec m2 = cumsum(m); 
//   arma::mat Xi(nm, p, arma::fill::zeros);
//   arma::vec Yi(nm, arma::fill::zeros);
//   // arma::vec Wi(nm, arma::fill::zeros);
//   arma::vec T0 = log(Y) + X * a;
//   int mn = m.n_elem;
//   for (int i = 0; i < mn; i ++) {
//     if (i == 0 && m(i) > 0) {
//       Yi.subvec(0, m2(i) - 1).fill(Y(i));
//       Xi.submat(0, 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
//     }
//     if (i > 0 && m(i) > 0) {
//       Yi.subvec(m2(i - 1), m2(i) - 1).fill(Y(i));
//       Xi.submat(m2(i - 1), 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
//     }
//   }
//   arma::vec texa = log(T) + Xi * a;
//   arma::vec yexa = log(Yi) + Xi * a;  
//   arma::vec Lam(n, arma::fill::zeros);
//   arma::vec de(nm, arma::fill::zeros);
//   for (int i = 0; i < nm; i++) {
//     for (int j = 0; j < nm; j++) {
//       if ((texa[i] <= yexa[j]) && (texa[i] >= texa[j])) {
// 	// de(i) += Wi[j];
// 	de(i) += 1;
//       }
//     }
//   }
//   for (int k = 0; k < n; k++) {
//     for (int i = 0; i < nm; i++) {
//       if (texa[i] >= T0[k] && de(i) > 0) {
//         // Lam[k] += Wi(i) / de(i);
// 	Lam[k] += 1 / de(i);
//       }
//     }
//   }
//   Lam = exp(-Lam); 
//   arma::vec R = m / Lam;
//   return ((W % R - mean(W % R))).t() * X / n;
// }

// Used for terminal events
//' @noRd
// [[Rcpp::export(rng = false)]]
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
// [[Rcpp::export(rng = false)]]
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
  // arma::mat nu = Iij * (X % repmat(ebax % Z % W, 1, X.n_cols));
  arma::mat nu = Iij * matvec(X, ebax % Z % W);
  arma::mat de = Iij * (ebax % Z % W);
  arma::mat tmp = nu / repmat(de, 1, nu.n_cols);
  tmp.replace(arma::datum::nan, 0);
  arma::mat D2 = repmat(D % W, 1, X.n_cols);
  arma::mat D3 = repmat(yexa % D % W, 1, X.n_cols);
  return join_rows(sum(X % D2, 0) - sum(tmp % D2, 0),
                   sum(X % D3, 0) - sum(tmp % D3, 0)) / n;
}

// [[Rcpp::export(rng = false)]]
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

// [[Rcpp::export(rng = false)]]
arma::rowvec temLog(const arma::vec& a,
		    const arma::vec& b,
		    const arma::mat& X,
		    const arma::vec& Y,
		    const arma::vec& Z,
		    const arma::vec& D,
		    const arma::vec& W) {
  int n = Y.n_elem;
  int p = X.n_cols;
  arma::vec yexa = Y % exp(X * a);
  arma::vec ebaxZ = Z % exp(X * (b - a));
  arma::uvec ind = stable_sort_index(yexa, "descend");
  arma::vec ordD = D(ind);
  arma::vec ordW = W(ind);
  arma::mat xz = X % repmat(ebaxZ, 1, p);
  xz = cumsum(xz.rows(ind), 0);
  // Rcpp::Rcout << ind;
  arma::mat c1 = X.rows(ind);
  arma::vec tmp = cumsum(ebaxZ(ind));
  arma::mat r = c1 - xz / repmat(tmp, 1, p);
  r.replace(arma::datum::nan, 0);
  return sum(repmat(ordW % ordD, 1, p) % r, 0) / n;
}

// [[Rcpp::export(rng = false)]]
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
