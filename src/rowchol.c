void chol(double *s, double *t, long *pd);

void rowchol(double *x, double *z, long *pn, long *pnx)
{

        long n, nx, l;

        n = *pn;
        nx = *pnx;

	if (nx==1){
	   for (l=0;l<n;l++)
	      *(z+l) = 1/(*(x+l));
	}
	if (nx>1){
	   for (l=0;l<n;l++)
	      chol((x+nx*nx*l), (z+nx*nx*l), pnx);
	}
}
