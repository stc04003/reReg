#include <R.h>
#include <Rmath.h>
#include <math.h>

void ghosh(double *Tik, double *Y, double *X, int *cl, int *clsz, int *n, int *p, double *res) {
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
	      nu[r] += X[j + r * *n];
	    }
	    de += 1;
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

void lwyy(double *Tik, double *Y, double *X, double *wgt, int *cl, int *clsz, int *n, int *p, double *res) {
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


