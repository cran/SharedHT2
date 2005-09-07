#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

void ludcmp(double *, int, int *, double *);
void lubksb(double *, int, int *, double *);
/*                                                        */
/*  pm : pointer to int, m.                              */
/*   a : pointer to double of length m^2                  */
/* yvv : pointer to a double of length m^2                */
/*       all space must be allocated in calling routine   */
/* NOTE: the routine writes over the contents of 'a'      */
/*                                                        */
void matinv(double *a, double *yvv, int *pm)
{
	double *yvcol, *dvv;
	int i,j, m, *indx;

	m = *pm;

	indx = (int *) Calloc(m, int);
	 dvv = (double *) Calloc(1, double);
 

	for (i=0;i<m;i++){
	  for (j=0;j<m;j++){
	    *(yvv+m*j+i) = 0.0;
	  }
	  *(yvv+m*i+i) = 1.0;
	}

	ludcmp(a, *pm, indx, dvv);

	for (j=0;j<m;j++){
	   yvcol = (yvv+m*j);
	   lubksb(a, *pm, indx, yvcol);
	}
	Free(indx);
	Free(dvv);
}

void ludcmp(double *a, int n, int *indx, double *pd)
{
	double tiny, *vv, aamax, sum, dum, d;
	int i, j, k, imax=0;

	   d = *pd;

	  vv = (double *) Calloc(n, double);

	tiny = pow(10,-20);
	d = 1.0;
	for (i=0;i<n;i++){
	   aamax = 0.0;
	   for (j=0;j<n;j++)
	      if (fabs(*(a+n*j+i)) > aamax) aamax = fabs(*(a+n*j+i));
	   *(vv+i) = 1.0/aamax;
	}

	for (j=0;j<n;j++){
	   if (j>0)
	   for (i=0;i<j;i++){
	      sum = *(a+j*n+i);
	      if (i>0)
	      for (k=0;k<i;k++)
		 sum = sum - (*(a+n*k+i))*(*(a+n*j+k));
	      *(a+n*j+i) = sum;
	   }
	   aamax = 0.0;

	   for (i=j;i<n;i++){
	      sum = *(a+n*j+i);
	      for (k=0;k<j;k++)
		 sum = sum - (*(a+n*k+i))*(*(a+n*j+k));
	      *(a+n*j+i) = sum;
	      dum = fabs(sum)*(*(vv+i));
	      if (dum>=aamax) {
		 imax = i;
		 aamax = dum;
	      }
	   }
	   if (j!=imax) {
	      for (k=0;k<n;k++){
		 dum = *(a+n*k+imax);
		 *(a+n*k+imax) = *(a+n*k+j);
		 *(a+n*k+j) = dum;
	      }
	      d = -1.0*d;
	      *(vv+imax) = *(vv+j);
	   }
	   *(indx+j) = imax;
	   if (*(a+n*j+j)==0.0) *(a+n*j+j) = tiny;
	   if (j!=(n-1)) {
	      dum = 1.0/(*(a+n*j+j));
	      for (i=(j+1);i<n;i++)
		*(a+n*j+i) = dum*(*(a+n*j+i));
	   }
	}
	Free(vv);
}

void lubksb(double *a, int n, int *indx, double *b)
{
	int ii, ll, i, j;
	double sum;

	ii = -1;
	for (i=0;i<n;i++){
	   ll = *(indx+i);
	   sum = *(b+ll);
	   *(b+ll) = *(b+i);
	   if (ii!=-1) {
	      for (j=ii;j<i;j++)
		 sum = sum - (*(a+n*j+i))*(*(b+j));
	   }
	   else if (sum != 0.0) ii = i;
	   *(b+i) = sum;
	}
	for (i=(n-1);i>=0;i--){
	   sum = *(b+i);
	   for (j=(i+1);j<n;j++)
	      sum = sum - (*(a+n*j+i))*(*(b+j));
	      *(b+i) = sum/(*(a+n*i+i));
	}
}
