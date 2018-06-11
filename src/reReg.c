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

// SARM 3
/* void sarm2(double *X, double *T, double *Y, double *weights, double *lambda,  */
/* 	   double *mt, int *n, int *p, int *B, */
/* 	   double *res) { */
/*   double *nu = Calloc(*p, double); */
/*   int i, j, r, b; */
/*   double de; */
/*   for (b = 0; b < *B; b++) { */
/*     for (i = 0; i < *n; i++) { */
/*       for (r = 0; r < *p; r++) { */
/* 	nu[r] = 0; */
/*       } */
/*       de = 0.0; */
/*       for (j = 0; j < *n; j++) { */
/* 	if (T[i] <= Y[j] & lambda[j] > 0) { */
/* 	  for (r = 0; r < *p; r++) { */
/* 	    // nu[r] += mt[j] * X[j + r * *n] / lambda[j]; */
/* 	    nu[r] += X[j + r * *n] / lambda[j]; */
/* 	  } */
/* 	  // de += mt[j] / lambda[j]; */
/* 	  de += 1 / lambda[j]; */
/* 	} */
/*       } */
/*       for (r = 0; r < *p; r ++) { */
/*       	res[r] += X[i + r * *n] - nu[r] / de; */
/*       	// res[r] += (X[i + r * *n] - X[j + r * *n]) * de[r]; */
/*       } */
/*     } */
/*     } */
/*   Free(nu); */
/*   res; */
/* } */

// SARM 2

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

/* void alphaEq1(double *X, double *Lambda, double *weights,  */
/* 	      int *mt, int *n, int *p, int *B,  */
/* 	      double *res) { */
/*   int i, j, r, b, iId = 0, jId = 0;  */
/*   for (b = 0; b < *B; b++) { */
/*     for (i = 0; i < *n; i++) { */
/*       for (j = 0; j < *n; j++) { */
/* 	for (r = 0; r < *p; r++) { */
/* 	  if (Lambda[i] != 0 && Lambda[j] != 0) */
/* 	    res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * weights[j + b * *n] * (mt[i] / Lambda[i] -  mt[j] / Lambda[j]); */
/* 	  if (Lambda[i] != 0 && Lambda[j] == 0) */
/* 	    res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * weights[j + b * *n] * (mt[i] / Lambda[i]); */
/* 	  if (Lambda[i] == 0 && Lambda[j] != 0) */
/* 	    res[r + b * *p] += X[i + r * *n] * weights[i + b * *n] * weights[j + b * *n] * (0 - mt[j] / Lambda[j]); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   res; */
/* }		  */


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

void HWb(double *Y, double *X, double *delta, double *z, double *weights,
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
	      nu[r] += exp(weights[j + b * *n]) * z[j] * X[j + r * *n];
	    }
	    de += exp(weights[j + b * *n]) * z[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de == 0) {
	    res[r + b * *p] += X[i + r * *n];
	  }
	  if (de != 0) {
	    res[r + b * *p] += X[i + r * *n] - (nu[r] / de);
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

// from aftsrr

// log-rank type estimating function (non-smooth); old name = ulognsfun
void log_ns_est(double *beta, double *Y, double *X, double *delta, int *clsize,
		int *n, int *p, int *N, double *weights, double *gw, double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, r;
  double *e = Calloc(*N, double), *nu = Calloc(*p, double);
  double de = 0.0;
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset nu */
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    if (e[ik_idx] - e[jl_idx] <= 0) {
	      for ( r = 0; r < *p; r++) {
		nu[r] += X[jl_idx + r * *N] * weights[jl_idx];
	      }
	      de += weights[jl_idx];
	    }
	    jl_idx++;
	  }
	}  // end jl
	for (r =  0; r < *p; r++) {
	  sn[r] += weights[ik_idx] * gw[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	}
      }
      ik_idx++;
    }
  }
  Free(nu);
  Free(e);
}
