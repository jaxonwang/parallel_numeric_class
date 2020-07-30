#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* matrix-vector multiply : y = A * x, where 
   A is symmetric and only lower half are stored */
void symMatVec(int n, double *a, double *x, double *y) {
  int i, j;
  
  for (i=0; i< n; i++) {
    double t = 0.0;
    for (j=0; j<= i; j++)
      t += a[i*n+j] * x[j];
    
    for (j= i+1; j< n; j++)
      t += a[j*n+i] * x[j];
    
    y[i] = t;
  }
}

/* solve Ax = b */
void solveSym(int n, double *a, double *x, double *b) {
  int i, j, k;
  
  /* LDLT decomposition: A = L * D * L^t */
  for (i=0; i< n; i++) {
    double invp = 1.0 / a[i*n+i];
    
    for (j= i+1; j< n; j++) {
      double aji = a[j*n+i];
      a[j*n+i] *= invp;
      
      for (k= i+1; k<= j; k++)
	a[j*n+k] -= aji * a[k*n+i];
    }
  }

  /* forward solve L y = b: but y is stored in x
     can be merged to the previous loop */
  for (i=0; i< n; i++) {
    double t = b[i];

    for (j=0; j< i; j++)
      t -= a[i*n+j] * x[j];

    x[i] = t;
  }

  /* backward solve D L^t x = y */
  for (i= n-1; i>= 0; i--) {
    double t = x[i] / a[i*n+i];

    for (j= i+1; j< n; j++)
      t -= a[j*n+i] * x[j];

    x[i] = t;
  }
}

int main(void) {

  int n, i, j;
  n = 30;

  /* matrix */
  double *a = malloc(sizeof(double) * n * n);
  assert(a != NULL);

  /* fill lower triangular elements */
  for (i=0; i< n; i++)
    for (j=0; j< i; j++)
      a[i*n+j] = rand()/(RAND_MAX + 1.0);

  /* fill diagonal elements */
  for (i=0; i< n; i++) {
    double s = 0.0;
    for (j=0; j< i; j++)
      s += a[i*n+j];
    
    for (j= i+1; j< n; j++)
      s += a[j*n+i];		/* upper triangular */
    
    a[i*n+i] = s + 1.0;		/* diagonal dominant */
  }

  /* first make the solution */
  double *xx = malloc(sizeof(double) * n);
  assert(xx != NULL);

  for (i=0; i< n; i++)
    xx[i] = 1.0;			/* or anything you like */

  /* make right hand side b = Ax */
  double *b = malloc(sizeof(double) * n);
  assert(b != NULL);

  symMatVec(n, a, xx, b);

  /* solution vector, pretend to be unknown */
  double *x = malloc(sizeof(double) * n);
  assert(x != NULL);

  /* solve: the main computation */
  solveSym(n, a, x, b);

  /* check error norm */
  double e = 0;
  for (i=0; i< n; i++)
    e += (x[i] - xx[i]) * (x[i] - xx[i]);
  e = sqrt(e);

  printf("error norm = %e\n", e);
  printf("--- good if error is around n * 1e-16 or less\n");
  
  return 0;
}
