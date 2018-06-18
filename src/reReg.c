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

void sarm1(double *X, double *Lambda, double *weights, double *gamma, 
	   int *mt, int *n, int *p, int *B, 
	   double *res) {
  int i, r, b, k; 
  double xr; 
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      xr = 0;
      for (k = 0; k < *p; k++) {
	xr += X[i + k * *n] * gamma[k];
      }
      for (r = 0; r < *p; r++) {
	if (Lambda[i] > 0) 
	  res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * (mt[i]/Lambda[i] - exp(xr));
	if (Lambda[i] == 0)
	  res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * (-1 * exp(xr));
      } // end r
    } // end i
  } // end b
}

void sarm2(double *X, double *T, double *Y, double *weights, 
	   int *n, int *p, int *B,
	   double *res) {
  double *nu = Calloc(*p, double);
  int i, j, r, b;
  double de;
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      for (r = 0; r < *p; r++) {
	nu[r] = 0;
      }
      de = 0.0;
      for (j = 0; j < *n; j++) {
	if (T[i] <= Y[j] && T[i] >= T[j]) {
	  for (r = 0; r < *p; r++) {
	    nu[r] += X[j + r * *n];
	  }
	  de += 1;
	}
      }
      for (r = 0; r < *p; r ++) {
      	res[r] += X[i + r * *n] - nu[r] / de;
      	// res[r] += (X[i + r * *n] - X[j + r * *n]) * de[r];
      }
    }
    }
  Free(nu);
}

void alphaEqC(double *X, double *Lambda, int *mt, int *n, int *p, double *res) {
  int i, j, r; 
  for (i = 0; i < *n; i++) {
    for (j = 0; j < *n; j++) {
      for (r = 0; r < *p; r++) {
	if (Lambda[i] != 0 && Lambda[j] != 0) {
	  res[r] += X[i + r * *n] * (mt[i] / Lambda[i] -  mt[j] / Lambda[j]);
	}
	if (Lambda[i] != 0 && Lambda[j] == 0) {
	  res[r] += X[i + r * *n] * (mt[i] / Lambda[i]);
	}
	if (Lambda[i] == 0 && Lambda[j] != 0) {
	  res[r] += X[i + r * *n] * (0 - mt[j] / Lambda[j]);
	}
      }
    }
  }
}		 

void betaEst(double *Y, double *X, double *delta, double *z, double *weights,
	     int *n, int *p, int *B, 
	     // output
	     double *res) {
  int i, j, r, b;
  double *nu = Calloc(*p, double); 
  double de;
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      if (delta[i] != 0) {
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	for (j = 0; j < *n; j++) {
	  if (Y[i] <= Y[j]) {
	    for (r = 0; r < *p; r++) {
	      nu[r] += weights[j + b * *n] * z[j] * X[j + r * *n];
	    }
	    de += weights[j + b * *n] * z[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de == 0) {
	    res[r + b * *p] += weights[i + b * *n] *  X[i + r * *n];
	  }
	  if (de != 0) {
	    res[r + b * *p] += weights[i + b * *n] * (X[i + r * *n] - (nu[r] / de));
	  }
	}
      } // end delta
    }
  }
  Free(nu);
}

void HWb(double *Y, double *X, double *delta, double *z, double *xb, double *weights, 
	     int *n, int *p, int *B, 
	     // output
	 double *res) {
  int i, j, r, b;
  double *nu = Calloc(*p, double); 
  double de;
  for (b = 0; b < *B; b++) {
    for (i = 0; i < *n; i++) {
      if (delta[i] != 0) {
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	for (j = 0; j < *n; j++) {
	  if (Y[i] <= Y[j]) {
	    for (r = 0; r < *p; r++) {
	      nu[r] += weights[j] * exp(xb[j + b * *n]) * z[j] * X[j + r * *n];
	    }
	    de += weights[j] * exp(xb[j + b * *n]) * z[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de == 0) {
	    res[r + b * *p] += weights[i] * X[i + r * *n];
	  }
	  if (de != 0) {
	    res[r + b * *p] += weights[i] * (X[i + r * *n] - (nu[r] / de));
	  }
	}
      } // end delta
    }
  }
  Free(nu);
}

void scaleChangeLog(int *n, int *p, int *start, int *M,
		    double *y, double *tij, double *X, double *W, double *result) {
  int i, j, k, l, r;
  double de;
  double *nu = Calloc(*p, double);
  for (i = 0; i < n[0]; i++) {
    for (j = 0; j < M[i]; j++) {
      for (k = 0; k < n[0]; k++) {
	for (l = 0; l < M[k]; l++) {
	  if (tij[start[k] + l] <= tij[start[i] + j] && tij[start[i] + j] <= y[start[k]]) {
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

void scaleChangeGehan(int*n, int *p, int *start, int *M, double *y,
		      double *tij, double *X, double *W, double *result) {
  int i, j, k, l, r;
  for (i = 0; i < n[0]; i++) {
    for (j = 0; j < M[i]; j++) {
      for (k = 0; k < n[0]; k++) {
	for (l = 0; l < M[k]; l++) {
	  if (tij[start[k] + l] <= tij[start[i] + j] && tij[start[i] + j] <= y[start[k]]) {
	    for (r = 0; r < p[0]; r++) {
	      result[r] += W[i] * W[k] * (X[i + r * n[0]] - X[k + r * n[0]]);
	    }	    
	  } // end if
	}
      }
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
