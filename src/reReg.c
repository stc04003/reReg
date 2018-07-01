#include <R.h>
#include <Rmath.h>
#include <math.h>

void plLambda(double *sl, double *tij, double *yi, double *weights,
	      int *n, int *N, // n is the length of sl
	      // output
	      double *res) {
  int i, j;
  double dl = 0, Rl = 0;
  for (j = 0; j < *n; j++) {
    res[j] = 1;
    dl = 0;
    Rl = 0;
    for (i = 0; i < *N; i++) {
      if (sl[j] >= tij[i] && sl[j] <= yi[i]) {
	Rl = Rl + weights[i];
      }
      if (sl[j] == tij[i] && tij[i] != yi[i]) {
	dl = dl + weights[i];
      }
    }
    if (Rl > 0) {
      if (j == 0)
	res[j] = (dl / Rl);
      if (j > 0)
	res[j] = res[j - 1] + (dl / Rl);
    }
  }
}

void sarm1(double *X, double *weights, double *xr, 
	   int *ratio, int *n, int *p, double *res) {
  int i, r; 
  for (i = 0; i < *n; i++) {
    for (r = 0; r < *p; r++) {
      res[r] += X[i + r * *n] * weights[i] * (ratio[i] - xr[i]);
    } 
  } 
}

void alphaEqC(double *X, double *ratio, int *n, int *p, double *weights, double *res) {
  int i, j, r;
  double snd;
  double wgtMean = 0;
  for (i = 0; i < *n; i++) {
    wgtMean += weights[i];
  }
  wgtMean = wgtMean / *n;
  for (i = 0; i < *n; i++) {
    snd = 0;
    for (j = 0; j < *n; j++) {
      snd += weights[j] * ratio[j];
    }
    for (r = 0; r < *p; r++) {
      res[r] += weights[i] * X[i + r * *n] * (wgtMean * ratio[i] - snd / *n);
    }
  }
}

void betaEst(double *Y, double *X, double *delta, double *z, double *weights,
	     int *n, int *p, double *res) {
  int i, j, r;
  double *nu = Calloc(*p, double); 
  double de;
  for (i = 0; i < *n; i++) {
    if (delta[i] > 0) {
      de = 0.0;
      for (j = 0; j < *n; j++) {
	if (Y[i] <= Y[j]) {
	  for (r = 0; r < *p; r++) {
	    nu[r] += weights[j] * z[j] * X[j + r * *n];
	  }
	  de += weights[j] * z[j];
	}
      }
      for (r = 0; r < *p; r++) {
	if (de == 0) {
	  res[r] += weights[i] *  X[i + r * *n];
	}
	if (de != 0) {
	  res[r] += weights[i] * (X[i + r * *n] - (nu[r] / de));
	}
	nu[r] = 0;
      }
    } // end delta
  }
  Free(nu);
}

void HWb(double *Y, double *X, double *delta, double *z, double *xb, double *weights, 
	 int *n, int *p, double *res) {
  int i, j, r;
  double *nu = Calloc(*p, double); 
  double de;
  for (i = 0; i < *n; i++) {
    if (delta[i] != 0) {
      for ( r = 0; r < *p; r++) {
	nu[r] = 0.0;
      }
      de = 0.0;
      for (j = 0; j < *n; j++) {
	if (Y[i] <= Y[j]) {
	  for (r = 0; r < *p; r++) {
	    nu[r] += weights[j] * exp(xb[j]) * z[j] * X[j + r * *n];
	  }
	  de += weights[j] * exp(xb[j]) * z[j];
	}
      }
      for (r = 0; r < *p; r++) {
	if (de == 0) {
	  res[r] += weights[i] * X[i + r * *n];
	}
	if (de != 0) {
	  res[r] += weights[i] * (X[i + r * *n] - (nu[r] / de));
	}
      }
    } // end delta
  }  
  Free(nu);
}

// \code{sc1Log} gives the first estimating equation in the scale-change model Xu et al. (2018).
// Equation 3.4:
// S_n(a) = (1/n) * sum_i sum_k [X_i - \frac{sum_j sum_l X_j I(...)}{sum_j sum_l X_j I(...)}
// This is in log-rank form.
// Notations are similar to that in \code{scRate}.
void sc1Log(int *n, int *p, int *start, int *M,
	    double *y, double *tij, double *X, double *W, double *result) {
  int i, j, k, l, r;
  double de;
  double *nu = Calloc(*p, double);
  for (i = 0; i < n[0]; i++) {
    for (j = 0; j < M[i]; j++) {
      for (k = 0; k < n[0]; k++) {
	for (l = 0; l < M[k]; l++) {
	  if (tij[start[k] + l] <= tij[start[i] + j] && tij[start[i] + j] <= y[k]) {
	    de = de + W[k];
	    for (r = 0; r < p[0]; r++) {
	      nu[r] += W[k] * X[k + r * n[0]]; 
	    }
	  } // end if
	}
      }
      for (r = 0; r < p[0]; r++) {
	if (de > 0) result[r] += W[i] * (X[i + r * n[0]] - nu[r] / de);
	nu[r] = 0;
      }
      de = 0;
    }
  }
  Free(nu);
}

// \code{sc1Gehan} gives the first estimating equation in the scale-change model Xu et al. (2018).
// Equation 3.4:
// S_n(a) = (1/n) * sum_i sum_k sum_j sum_l (X_i - X_j) * I(...)
// This is in Gehan form.
// Notations are similar to that in \code{scRate}.
void sc1Gehan(int*n, int *p, int *start, int *M, double *y,
		      double *tij, double *X, double *W, double *result) {
  int i, j, k, l, r;
  for (i = 0; i < n[0]; i++) {
    for (j = 0; j < M[i]; j++) {
      for (k = 0; k < n[0]; k++) {
	for (l = 0; l < M[k]; l++) {
	  if (tij[start[k] + l] <= tij[start[i] + j] && tij[start[i] + j] <= y[k]) {
	    for (r = 0; r < p[0]; r++) {
	      result[r] += W[i] * W[k] * (X[i + r * n[0]] - X[k + r * n[0]]);
	    }	    
	  } // end if
	}
      }
    }
  }
}

// \code{sc2} gives the 2nd estimating equation in the scale-change model Xu et al. (2018).
// Equation 3.5
// (1/n) sum_i \bar{X} [\frac{m_i}{\Lambda} - exp(\bar{X} \theta)
// Notations similar to that in scRate
void sc2(int *n, int *p, double *X, double *xb, 
	 double *ratio, double *W, double *result) {
  int i, r;
  for (i = 0; i < *n; i++) {
    for (r = 0; r < *p; r++) {
      result[r] += W[i] * X[i + r * *n] * (ratio[i] - xb[i]);
    }
  }
}

// \code{scRate} gives rates for the scale-change model (method = sc.XCYH)
// From the paper, this is
// \hat{H}_n(t; \hat\alpha) = \int_t^\tau \frac{\sum_{i=1}^n dN*....}{\sum_{i=1}^nR_i^*...}
//
// notations are defined similarly as that in glRate
void scRate(int *n, int *start, int *M, int *nt0, double *W,
	    double *yi, double *tij, double *t0, double *result) {
  int i, j, k, l, r;
  int nM = 0;
  for (i = 0; i < *n; i++) nM += M[i];
  double *risk = Calloc(nM, double);
  for (i = 0; i < nM; i++) {
    for (j = 0; j < *n; j++) {
      if (tij[i] <= yi[j]) {
	for (l = 0; l < M[j]; l++) {
	  risk[i] += W[j] * (tij[start[j] + l] <= tij[i]);
	}
      }
    }
  }
  for (r = 0; r < *nt0; r++) {
    for (i = 0; i < *n; i++) {
      for (k = 0; k < M[i]; k++) {
	if (tij[start[i] + k] >= t0[r] && risk[start[i] + k] > 0) result[r] += W[i] / risk[start[i] + k];
      }
    }
  }
  Free(risk);
}

/* void scRate(int *n, int *start, int *M, int *nt0, double *W, */
/* 	    double *yi, double *tij, double *t0, double *result) { */
/*   int i, j, k, l, r; */
/*   double de = 0; */
/*   for (r = 0; r < *nt0; r++) { */
/*     for (i = 0; i < *n; i++) { */
/*       for (k = 0; k < M[i]; k++) { */
/* 	if (tij[start[i] + k] >= t0[r] && tij[start[i] + k] <= yi[i]) { */
/* 	  for (j = 0; j < *n; j++) { */
/* 	    for (l = 0; l < M[j]; l++) { */
/* 	      if (tij[start[j] + l] <= tij[start[i] + k] && tij[start[i] + k] <= yi[j]) de += W[j]; */
/* 	    } */
/* 	  } */
/* 	  if (de > 0) result[r] += W[i] / de; */
/* 	  de = 0; */
/* 	} */
/*       } */
/*     } */
/*   } */
/* } */


// cumulative baseline function in cox.hw
// zhat = Z_j * exp(X %*% beta)
void hwHaz(double *t0, double *yi, double *zi, double *delta, double *wgt,
	   int *n, int *nt0, double *result) {
  int i, j, r;
  double de = 0;
  for (r = 0; r < *nt0; r++) {
    for (i = 0; i < *n; i++) {
      if (delta[i] > 0 && yi[i] <= t0[r]) {
	for (j = 0; j < *n; j++) {
	  de += zi[j] * wgt[j] * (yi[i] <= yi[j]);
	}
      }
      if (de > 0) result[r] += wgt[i] / de;
      de = 0;
    }
  }
}


// C version of outer
// This function minics R's outer(t1, t2, "<=") but returns a vecter.
void outerC1(double *t1, double *t2, int *n1, int*n2, double *result) {
  int i, j, ind;
  ind = 0;
  for (j = 0; j < *n2; j++) {
    for (i = 0; i < *n1; i++) {
      result[ind] = (t1[i] <= t2[j]);
      ind += 1;
    }
  }  
}

// C version of outer
// This function minics R's outer(t1, t2, ">=") but returns a vecter.
void outerC2(double *t1, double *t2, int *n1, int*n2, double *result) {
  int i, j, ind;
  ind = 0;
  for (j = 0; j < *n2; j++) {
    for (i = 0; i < *n1; i++) {
      result[ind] = (t1[i] >= t2[j]);
      ind += 1;
    }
  }  
}
