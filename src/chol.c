#include <R.h>
#include <Rmath.h>
#define MDXMIN 2.470328e-323

void chol(double *s, double *t, int *pd)
{
	int d, d2, i, j, k;
	double sum, ansij;

	d = *pd;
	d2 = d*d;

	for (i=0;i<d2;i++) *(t+i) = 0.0;

	for (i=0;i<d;i++){
	   sum = 0.0;

	   for (k=0;k<i;k++)
	      sum = sum + pow(*(t+d*i+k),2);
	   *(t+d*i+i) = pow(*(s+d*i+i) - sum, 0.5);
	   for (j=(i+1);j<d;j++){
	      sum = 0.0;
	      for (k=0;k<i;k++)
	         sum = sum + (*(t+d*i+k))*(*(t+d*j+k));
              ansij = (*(s+d*j+i) - sum)/(*(t+d*i+i));
	      *(t+d*j+i) = (R_FINITE(ansij) ? ansij : MDXMIN);
	   }
	}
}
