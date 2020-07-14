#include <R.h>
#include <Rmath.h>
#include <math.h>

// \code{glU2} gives equation U2 in GL (2003) [eq. 4].
//
// n: the number of id
// p: the dimension of X
// start: the index pointing to the starting of tij
// M: the # of recurrent event times (tij) in the ith subject (cluster size)
// yi: the transformed censoring (terminal) times of size n || Y_i * exp(-\eta * X - d)
// tij: the transformed recurrent event times || t_ij * exp(-\theta * X)
// X: the covariate matrix of dimension p and n rows
// results: what to return
void glU2(int *n, int *p, int *start, int *M,
	  double *yi, double *tij, double *X,
	  double *weight, double *result) {
  int i, j, k, r;
  double *nu = Calloc(*p, double); 
  double de;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < M[i]; k++) {
      if(tij[start[i] + k] <= yi[i]) {
	for (j = 0; j < *n; j++) {
	  if (tij[start[i] + k] <= yi[j]) {
	    de += weight[j];
	    for (r = 0; r < *p; r++) {
	      nu[r] += weight[j] * X[j + r * *n];
	    }
	  }
	} // end j
	for (r = 0; r < *p; r++) {
	  if (de > 0) result[r] += weight[i] * (X[i + r * *n] - nu[r] / de);
	  if (de <= 0) result[r] += weight[i] * X[i + r * *n];
	  nu[r] = 0;
	}
	de = 0;
      }
    }
  }
  Free(nu);
}

// \code{glRate} gives rates in GL (2003) [the integral term in \hat R_0].
// From the paper, this is
// \sum_{k = 1}^{m_i}\frac{I(...)}{\sum_{j = 1}^nI(...)}
//
// notations are defined similarly as that in glU2,
// except that *result returns a vector with length equals to all possible jumps, length(unique(tij))
void glRate(int *n, int *start, int *M, int *nt0, 
	    double *yi, double *tij, double *t0, double *result) {
  int i, j, k, r;
  double de = 0;
  for (r = 0; r < *nt0; r++) {
    for (i = 0; i < *n; i++) {
      for (k = 0; k < M[i]; k++) {
	if(tij[start[i] + k] <= yi[i] && tij[start[i] + k] <= t0[r]) {
	  for (j = 0; j < *n; j++) {
	    if (tij[start[i] + k] <= yi[j]) de += 1;
	  }
	  if (de > 0) result[r] += 1 / de;
	  de = 0;
	}
      } // end k
    }
  }
}

// \code{glHaz} gives hazard in GL (2003) [the integral term in \hat \Lambda]
// From the paper, this is
// \frac{\Delta_i}{\sum_{j = 1}^nI(...)}
//
// notations similar to that in glRate
void glHaz(int *n, int *status, int *ny0, double *yi, double *y0, double *result) {
  int i, j, r;
  double de = 0;
  for (r = 0; r < *ny0; r++) {
    for (i = 0; i < *n; i++) {
      if (status[i] == 1 && yi[i] <= y0[r]) {
	for (j = 0; j < *n; j++) {
	  if (yi[j] >= yi[i]) de += 1;
	}
	if (de > 0) result[r] += 1 / de;
	de = 0;
      }
    } // end i
  }
}


// from aftsrr
// log-rank type estimating function (non-smooth)
// used for regression with method = am.GL
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
		nu[r] += X[jl_idx + r * *N] * weights[j];
	      }
	      de += weights[j];
	    }
	    jl_idx++;
	  }
	}  // end jl
	for (r =  0; r < *p; r++) {
	  sn[r] += weights[i] * gw[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	}
      }
      ik_idx++;
    }
  }
  Free(nu);
  Free(e);
}

/* void lwyy(double *Tik, double *Y, double *X, double *wgt, int *cl, int *clsz, */
/* 	  int *n, int *p, double *res) { */
/*   int i, j, k, r; */
/*   double *nu = Calloc(*p, double); */
/*   double de; */
/*   for (i = 0; i < *n; i++) { */
/*     for (k = 0; k < cl[i]; k++) { */
/*       if (Y[i] >= Tik[clsz[i] + k]) { */
/* 	de = 0.0; */
/* 	for (r = 0; r < *p; r++) { */
/* 	  nu[r] = 0; */
/* 	} */
/* 	for (j = 0; j < *n; j++) { */
/* 	  if (Y[j] >= Tik[clsz[i] + k]) { */
/* 	    for (r = 0; r < *p; r++) { */
/* 	      nu[r] += X[j + r * *n] * wgt[j]; */
/* 	    } */
/* 	    de += wgt[j]; */
/* 	  } */
/* 	} */
/* 	for (r = 0; r < *p; r++) { */
/* 	  if (de > 0) { */
/* 	    res[r] += X[i + r * *n] - (nu[r] / de); */
/* 	  } */
/* 	  if (de <= 0) { */
/* 	    res[r] += X[i + r * *n]; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   Free(nu); */
/* } */

// Equation 8 of Ghosh & Lin (2002); Marginal regression models for recurrent and terminal events.
// IPSW estimator (ARF in Luo et al (2015)
// The weight 'matrix', w_i(t_ij), is a n by length(m) matrix.
// The ith column gives w_i and the jth row evaluates w_i at t_ij
void coxGL(double *Tik, double *Y, double *X, double *xb, double *wgt,
	   int *len_Tik, int *cl, int *clsz, int *n, int *p, double *res) {
  int i, j, k, r;
  double *nu = Calloc(*p, double);
  double de;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < cl[i]; k++) {
      if (Y[i] >= Tik[clsz[i] + k]) {
	de = 0.0;
	for (r = 0; r < *p; r++) {
	  nu[r] = 0;
	}
	for (j = 0; j < *n; j++) {
	  if (Y[j] >= Tik[clsz[i] + k]) {
	    for (r = 0; r < *p; r++) {
	      nu[r] += X[j + r * *n] * wgt[j * *len_Tik + clsz[i] + k] * xb[j];
	    }
	    de += wgt[j * *len_Tik + clsz[i] + k] * xb[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de > 0) {
	    res[r] += wgt[i * *len_Tik + clsz[i] + k] * (X[i + r * *n] - (nu[r] / de));
	  }
	  if (de <= 0) {
	    res[r] += wgt[i * *len_Tik + clsz[i] + k] * X[i + r * *n];
	  }
	}
      }
    }
  }
  Free(nu);
}


// Equation 9 of Ghosh & Lin (2002); Marginal regression models for recurrent and terminal events.
// Baseline rate function
// arguments follow similar structure in `coxGL`
void glCoxRate(double *Tik, double *Y, double *xb, double *wgt, double *T0, int *len_T0,
	       int *len_Tik, int *cl, int *clsz, int *n, double *res) {
  int i, j, k, r;
  double de;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < cl[i]; k++) {
      if (Y[i] >= Tik[clsz[i] + k]) {
	de = 0.0;
	for (j = 0; j < *n; j++) {
	  if (Y[j] >= Tik[clsz[i] + k]) {
	    de += wgt[j * *len_Tik + clsz[i] + k] * xb[j];
	  }
	}
	for (r = 0; r < *len_T0; r++) {
	  if (T0[r] >= Tik[clsz[i] + k]) {
	    res[r] += wgt[i * *len_Tik + clsz[i] + k] / de;
	  }
	}
      }
    }
  }
}

