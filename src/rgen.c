#include<R.h>
#include<Rmath.h>
void rnormn(int *pn, double *ans);
void rnormns(int *pn, double *ans);
void rgamman(int *pn, double *shape, double *scale, double *ans);
void rgammans(int *pn, double *shape, double *scale, double *ans);
void trwish1(double *pdf, int *pd, double *pSqrtSigma, double *pW);
void rwishart1(double *pdf, int *pd, double *pSqrtSigma, double *pW);
void rwishart1s(double *pdf, int *pd, double *pSqrtSigma, double *pW);
void rwishartn(int *pn, double *pdf, int *pd, double *pSqrtSigma, int *prows,
	       double *pW);

void rwishartns(int *pn, double *pdf, int *pd, double *pSqrtSigma, int *prows, 
		double *pW)
{
  int n, d, d2, i, j, k, h, irows;
  double df, sm, *Z, *Zsig;
  n = *pn;
  df = *pdf;
  d = *pd;
  d2 = d*d;
  irows = *prows;

  Z = (double *)Calloc(d2, double);
  Zsig = (double *)Calloc(d2, double);

  GetRNGstate();
  for (i=0;i<n;i++){
    for(j=0;j<d2;j++){
      *(Z+j) = 0.0;
      *(Zsig+j) = 0.0;
    }

    for (j=0;j<d;j++){
      *(Z + d*j + j) = pow(rgamma(0.5*(df- (double) j),2),0.5);
      for(k=(j+1);k<d;k++) 
	*(Z + d*k + j) = rnorm(0.0, 1.0);
    }

    for(j=0;j<d;j++) 
      for(k=0;k<d;k++){
	sm = 0.0;
        for(h=0;h<d;h++)
	  sm += *(Z +d*h +j) * (*(pSqrtSigma +d*k +h));
	*(Zsig +d*k +j) = sm;
      }

    if(irows==1)
      for(j=0;j<d;j++) 
        for(k=0;k<d;k++){
	  sm = 0.0;
          for(h=0;h<d;h++)
	    sm += *(Zsig +d*j +h +i) *(*(Zsig +d*k +h));
	  *(pW +n*d*k +n*j +i) = sm;
	}

    if(irows==0)
      for(j=0;j<d;j++) 
        for(k=0;k<d;k++){
	  sm = 0.0;
          for(h=0;h<d;h++)
	    sm += *(Zsig +d*j +h) *(*(Zsig +d*k +h));
	  *(pW +d*d*i +d*k +j) = sm;
	}
  }
  PutRNGstate();
  Free(Z);
  Free(Zsig);
}

void rwishartn(int *pn, double *pdf, int *pd, double *pSqrtSigma, int *prows,
                double *pW)
{
  int n, d, d2, i, j, k, h, irows;
  double df, sm, *Z, *Zsig;
  n = *pn;
  df = *pdf;
  d = *pd;
  d2 = d*d;
  irows = *prows;

  Z = (double *)Calloc(d2, double);
  Zsig = (double *)Calloc(d2, double);

  for (i=0;i<n;i++){

    for(j=0;j<d2;j++){
      *(Z+j) = 0.0;
      *(Zsig+j) = 0.0;
    }

    for (j=0;j<d;j++){
      *(Z + d*j + j) = pow(rgamma(0.5*(df- (double) j),2),0.5);
      for(k=(j+1);k<d;k++)
        *(Z + d*k + j) = rnorm(0.0, 1.0);
    }
    for(j=0;j<d;j++)
      for(k=0;k<d;k++){
	sm = 0.0;
        for(h=0;h<d;h++)
          sm += *(Z +d*h +j) * (*(pSqrtSigma +d*k +h));
	*(Zsig +d*k +j) = sm;
      }

    if(irows==1)
      for(j=0;j<d;j++)
        for(k=0;k<d;k++){
	  sm = 0.0;
          for(h=0;h<d;h++)
            sm += *(Zsig +d*j +h) *(*(Zsig +d*k +h));
	  *(pW +n*d*k +n*j +i) = sm;
	}

    if(irows==0)
      for(j=0;j<d;j++)
        for(k=0;k<d;k++){
	  sm = 0.0;
          for(h=0;h<d;h++)
            sm += *(Zsig +d*j +h) *(*(Zsig +d*k +h));
	  *(pW +d*d*i +d*k +j) = sm;
	}
  }
  Free(Z);
  Free(Zsig);
}

void rwishart1s(double *pdf, int *pd, double *pSqrtSigma, double *pW)
{
  int d, d2, i, j, k, h;
  double df, sm, *Z, *Zsig;
  df = *pdf;
  d = *pd;
  d2 = d*d;

  Z = (double *)Calloc(d2, double);
  Zsig = (double *)Calloc(d2, double);

  for(j=0;j<d2;j++){
    *(Z+j) = 0.0;
    *(Zsig+j) = 0.0;
  }

  GetRNGstate();
  for (j=0;j<d;j++){
    *(Z + d*j + j) = pow(rgamma(0.5*(df- (double) j),2),0.5);
    for(k=(j+1);k<d;k++) 
      *(Z + d*k + j) = rnorm(0.0, 1.0);
  }

  for(i=0;i<d2;i++) {
    *(Zsig + i) = 0.0;
    *(pW + i) = 0.0;
  }

  for(j=0;j<d;j++) 
    for(k=0;k<d;k++){
      sm = 0.0;
      for(h=0;h<d;h++)
        sm += *(Z + d*h + j) * (*(pSqrtSigma + d*k + h));
      *(Zsig + d*k + j) = sm;
    }

  for(j=0;j<d;j++) 
    for(k=0;k<d;k++){
      sm = 0.0;
      for(h=0;h<d;h++)
	sm += *(Zsig + d*j + h) * (*(Zsig + d*k + h));
      *(pW + d*k + j) = sm;
    }

  PutRNGstate();
  Free(Z);
  Free(Zsig);
}

void rwishart1(double *pdf, int *pd, double *pSqrtSigma, double *pW)
{
  int d, d2, j, k, h;
  double df, sm, *Z, *Zsig;
  df = *pdf;
  d = *pd;
  d2 = d*d;

  Z = (double *)Calloc(d2, double);
  Zsig = (double *)Calloc(d2, double);

  for(j=0;j<d2;j++){
    *(Z+j) = 0.0;
    *(Zsig+j) = 0.0;
  }

  for (j=0;j<d;j++){
    *(Z + d*j + j) = pow(rgamma(0.5*(df- (double) j),2),0.5);
    for(k=(j+1);k<d;k++){
      *(Z + d*k + j) = rnorm(0.0, 1.0);
    }
  }

  for(j=0;j<d;j++)
    for(k=0;k<d;k++){
      sm = 0.0;
      for(h=0;h<d;h++)
        sm += *(Z + d*h + j) * (*(pSqrtSigma + d*k + h));
      *(Zsig + d*k + j) = sm;
    }

  for(j=0;j<d;j++)
    for(k=0;k<d;k++){
      sm = 0.0;
      for(h=0;h<d;h++)
        sm += *(Zsig + d*j + h) * (*(Zsig + d*k + h));
      *(pW + d*k + j) = sm;
    }

  Free(Z);
  Free(Zsig);
}

void trwish1(double *pdf, int *pd, double *pSqrtSigma, double *pW)
{
  GetRNGstate();
  rwishart1(pdf, pd, pSqrtSigma, pW);
  PutRNGstate();
}

void rgammans(int *pn, double *shape, double *scale, double *ans)
{
  int n, i;
  n = *pn;
  GetRNGstate();
  for(i=0;i<n;i++) *(ans+i) = rgamma(*shape, *scale);
  PutRNGstate();
}

void rgamman(int *pn, double *shape, double *scale, double *ans)
{
  int n, i;
  n = *pn;
  for(i=0;i<n;i++) *(ans+i) = rgamma(*shape, *scale);
}

void rnormns(int *pn, double *ans)
{
	int i,n;
	n = *pn;

        GetRNGstate();
	for (i=0;i<n;i++)
	{
		*(ans+i) = rnorm(0.0, 1.0);
	}
	PutRNGstate();
}

void rnormn(int *pn, double *ans)
{
  int i,n;
  n = *pn;

  for (i=0;i<n;i++)
    *(ans+i) = rnorm(0.0, 1.0);
}
