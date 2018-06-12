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
	  double *yi, double *tij, double *X, double *result) {
  int i, j, k, r;
  double *nu = Calloc(*p, double); 
  double de;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < M[i]; k++) {
      if(tij[start[i] + k] <= yi[i]) {
	for (j = 0; j < *n; j++) {
	  if (tij[start[i] + k] <= yi[j]) {
	    de = de + 1;
	    for (r = 0; r < *p; r++) {
	      nu[r] += X[j + r * *n];
	    }
	  }
	} // end j
	for (r = 0; r < *p; r++) {
	  if (de > 0) result[r] += X[i + r * *n] - nu[r] / de;
	  if (de <= 0) result[r] += X[i + r * *n];
	  nu[r] = 0;
	}
	de = 0;
      }
    }
  }
  Free(nu);
}

// \code{glRate} gives individual rates in GL (2003) [the integral term in \hat R_0].
// From the paper, this is
// \sum_{k = 1}^{m_i}\frac{I(...)}{\sum_{j = 1}^nI(...)}
//
// notations are defined similarly as in glU2,
// except that *result returns a vector with length equals to all possible jumps, length(unique(tij))
void glRate(int *n, int *p, int *start, int *M, int *nt0, 
	    double *yi, double *tij, double *X, double *t0, 
	    double *result) {
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

// \code{glHaz} gives individual hazard in GL (2003) [the integral term in \hat \Lambda]
// From the paper, this is
// \frac{\Delta_i}{\sum_{j = 1}^nI(...)}
//
// notations similar to glRate
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

void lwyy(double *Tik, double *Y, double *X, double *wgt, int *cl, int *clsz,
	  int *n, int *p, double *res) {
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
	      nu[r] += X[j + r * *n] * wgt[j];
	    }
	    de += wgt[j];
	  }
	}
	for (r = 0; r < *p; r++) {
	  if (de > 0) {
	    res[r] += X[i + r * *n] - (nu[r] / de);
	  }
	  if (de <= 0) {
	    res[r] += X[i + r * *n];
	  }
	  // nu[r] = 0;
	}
      }
    }
  }
  Free(nu);
}


/* void lwyy(double *Tik, double *Y, double *X, double *wgt, int *cl, int *clsz, int *n, int *p, double *res) { */
/*   int i, j, k, l, r;  */
/*   double *nu = Calloc(*p, double);  */
/*   double de; */
/*   for (i = 0; i < *n; i++) { */
/*     for (k = 0; k < cl[i]; k++) { */
/*       //if (Y[i] >= Tik[clsz[i] + k]) { */
/* 	de = 0.0; */
/* 	for (r = 0; r < *p; r++) { */
/* 	  nu[r] = 0; */
/* 	}   */
/* 	for (j = 0; j < *n; j++) { */
/* 	  for (l = 0; l < cl[j]; l++) { */
/* 	    if (Tik[clsz[j] + l] >= Tik[clsz[i] + k]) { */
/* 	      for (r = 0; r < *p; r++) { */
/* 		nu[r] += X[j + r * *n] * wgt[j]; */
/* 	      } */
/* 	      de += wgt[j]; */
/* 	    } */
/* 	  } */
/* 	} */
/* 	for (r = 0; r < *p; r++) { */
/* 	  if (de > 0) { */
/* 	    res[r] += X[i + r * *n] - (nu[r] / de); */
/* 	  } */
/* 	  if (de <= 0) { */
/* 	    res[r] += X[i + r * *n]; */
/* 	  }	   */
/* 	  // nu[r] = 0; */
/* 	} */
/* 	//} */
/*     } */
/*   } */
/*   Free(nu); */
/*   res; */
/* } */


