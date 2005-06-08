void rowprod(double *x, double *y, double *z, long *pn, long *pna1, long *pnb1, long *pnb2){

        long n, na1, nb1, nb2, l, i, j, k;
        double s;
	n = *pn;
        na1 = *pna1;
	nb1 = *pnb1;
	nb2 = *pnb2;
	for (l=0;l<n;l++)
	for (i=0;i<na1;i++)
	for (j=0;j<nb2;j++){
	s = 0.0;
	for (k=0;k<nb1;k++)
	s = s + *(x+n*na1*k+n*i+l)* *(y+n*nb1*j+n*k+l);
	*(z+n*na1*j+n*i+l) = s;
	}
}

