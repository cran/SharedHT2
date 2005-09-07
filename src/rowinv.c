void matinv(double *a, double *yvv, int *pm);

void rowinv(double *x, double *z, int *pn, int *pnx)
{

        int n, nx, l;

        n = *pn;
        nx = *pnx;

	if (nx==1){
	   for (l=0;l<n;l++)
	      *(z+l) = 1/(*(x+l));
	}
	if (nx>1){
	   for (l=0;l<n;l++)
	      matinv((x+nx*nx*l), (z+nx*nx*l), pnx);
	}
}
