#include<R.h>
#include<Rmath.h>
#define EPS 1.0e-7
#define MDXMIN 2.470328e-323
typedef struct{
  double *MVM;
  int *pN;
  int *pd;
  int *nreps;
} Data;

typedef double optimfn(int n, double *par, void *ex);
typedef void optimgr(int n, double *par, double *gr, void *ex);

void tloglik(double *ptheta, double *MVM, int *pN, int *pd, 
             int *pnreps, double *pans);

void tGloglik(double *ptheta, double *MVM, int *pN, int *pd, 
              int *pnreps, double *pG);

void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
           int *fail, double abstol, double intol, void *ex,
           double alpha, double beta, double gamma, int trace,
           int *fncount, int maxit);

void vmmin(int n, double *x, double *Fmin, optimfn fn, optimgr gr, 
           int maxit, int trace, int *mask, double abstol, 
           double reltol, int nREPORT, void *ex, int *fncount, 
           int *grcount, int *fail);

void fHESS(double *x, Data *y, double *G, double *H, optimgr *grad);

double lmgamma(double a, int m);
double dimgamma(double a, int m);
double trimgamma(double a, int m);
double det(double *x, double *xd2buff, int *pd);
void chol(double *s, double *t, int *pd);
void matinv(double *a, double *ainv, int *pd);
void printmat(double *pA, int nr, int nc, char *name);

optimfn loglik;
optimgr Gloglik, *grad;

void Fit_MVF(double *ptheta0, double *MVM, int *pN, int *pd,
		int *pnreps, int *pverbose, double *objval, 
                double *estimate, int *fail, int *fncnt, int *grcnt, 
                int *mask, int *usegr, double *G, double *H)
{
  int npar, d, d2, N, i;
  int inpar, *ifail, *ifncnt, *igrcnt, *imask, verb;
  double xd, xN;
  Data *y;

  verb = (int) *pverbose;
  d = *pd;
  d2 = d*d;
  xd = (double) d;
  N = *pN;
  xN = (double) N;
  npar = d*(d+1)/2 + 1;
  inpar = (int) npar;

  ifail = (int *) fail;
  ifncnt = (int *) fncnt;
  igrcnt = (int *) grcnt;
  imask = (int *) mask;

  y = (Data *)Calloc(1, Data);

  y->MVM = MVM;
  y->pN = pN;
  y->pd = pd;
  y->nreps = pnreps;

  *estimate = 0.0;

  for(i=1;i<npar;i++) *(estimate+i) = 0.0;

  if(*usegr==0){
    nmmin(inpar, ptheta0, estimate, objval, loglik, ifail, -1e200, 1e-8, y,
          1.0, 0.5, 2.0, verb, ifncnt, 10000);
  }
  else{
    vmmin(inpar, ptheta0, objval, loglik, Gloglik, 10000, verb, imask, -1e200,
          1e-8, 10, y, ifncnt, igrcnt, ifail);

    for(i=0;i<npar;i++)
      *(estimate+i) = *(ptheta0+i);
  }

  Gloglik(inpar, estimate, G, y);
  grad = &Gloglik;
  fHESS(estimate, y, G, H, grad);

  Free(y);
}

double loglik(int p, double *theta, void *yy)
{
  int N, d, d2, h, i, j, k, l;
  double xd, xnreps, nu, detL, detS, detSplL, u, v, logl, logK;

  int *nreps, *pd;
  double *MVM, *Lambdahlf, *Lambda, *xd2buff, *SplL;
  Data *y;

  y = (Data *)yy;
  MVM = y->MVM;
  N = *(y->pN);
  pd = (y->pd);
  d = *pd;
  nreps = y->nreps;
  d2 = d*d;
  xd = (double) d;

  Lambdahlf = (double *)Calloc(d2, double);
  Lambda    = (double *)Calloc(d2, double);
  xd2buff   = (double *)Calloc(d2, double);
  SplL      = (double *)Calloc(d2, double);

  for(i=0;i<d2;i++) {
    *(Lambdahlf +i) = 0.0;
    *(Lambda + i) = 0.0;
  }

  nu = (2.0*xd + 2.0)*(exp(*theta)+1.0);

  l=0;
  for(i=0;i<d;i++){
    *(Lambdahlf + d*i + i) = exp(*(theta + i + 1));
    for(j=(i+1);j<d;j++) {
      *(Lambdahlf + d*j + i) = *(theta + d + 1 + l);
      l++;
    }
  }

  for(i=0;i<d;i++)
    for(j=0;j<d;j++)
      for(k=0;k<d;k++) *(Lambda + d*j + i) += *(Lambdahlf + d*i + k) *(
                                              *(Lambdahlf + d*j + k));

  detL = det(Lambda, xd2buff, pd);
  logl = 0.0;
  for(h=0;h<N;h++){
    xnreps = (double) (*(nreps +h));
    logK = lmgamma((nu + xnreps - xd -2.0)/2.0, d) - lmgamma((xnreps - 1.0)/2.0, d) -
           lmgamma((nu - xd - 1.0)/2.0, d);
    for (i=0;i<d2;i++) *(xd2buff+i) = 0.0;
    detS = pow(xnreps - 1.0, xd) * det((MVM+d2*h), xd2buff, pd);
    for (i=0;i<d2;i++) {
      *(xd2buff+i) = 0.0;
      *(SplL+i) = 0.0;
    }
    for(i=0;i<d2;i++) *(SplL +i) = (xnreps-1.0)*(*(MVM +d2*h +i)) +(*(Lambda +i));

    detSplL = det(SplL, xd2buff, pd);
    u = detS/detSplL;
    v = detL/detSplL;
    logl += logK + (xnreps-xd-2.0)/2.0 * log(u) + (nu-xd-1.0)/2.0 * log(v) - (xd+1.0)/2.0 * log(detSplL);
  }

  Free(Lambdahlf);
  Free(Lambda);
  Free(xd2buff);
  Free(SplL);

  return logl*(-1.0);
}

void Gloglik(int inpar, double *theta, double *G, void *yy)
{
  int   N, d, d2, d4, npar, npar2, nparm1;
  int   j, k, l, h, i1, i2, j2;
  double xd, xnreps, nu, detL, detS, detSplL, v, eth0;

  int   *nreps, *pd;
  double *MVM, *Lambdahlf, *Lambda, *Lambdainv, *xd2buff;
  double *SplL, *SplLinv, *dvecBdtheta, *dvecBprdtheta, *FactorII;
  double *Prod, *IxBpr, *BprxI, *IxLi, *IxSplLi, *FactorI;
  Data *y;

  y = (Data *)yy;
  MVM = y->MVM;
  N = *(y->pN);
  pd = (y->pd);
  d = *pd;
  nreps = y->nreps;
  npar = d*(d+1)/2 + 1;
  npar2 = npar*npar;
  nparm1 = npar - 1;
  d2 = d*d;
  d4 = d2*d2;
  xd = (double) d;

  Lambdahlf     = (double *)Calloc(       d2, double); /* d x d          */
  Lambda        = (double *)Calloc(       d2, double); /* d x d          */
  Lambdainv     = (double *)Calloc(       d2, double); /* d x d          */
  xd2buff       = (double *)Calloc(       d2, double); /* d x d          */
  SplL          = (double *)Calloc(       d2, double); /* d x d          */
  SplLinv       = (double *)Calloc(       d2, double); /* d x d          */
  dvecBdtheta   = (double *)Calloc(d2*nparm1, double); /* d^2 x (npar-1) */
  dvecBprdtheta = (double *)Calloc(d2*nparm1, double); /* d^2 x (npar-1) */
  FactorII      = (double *)Calloc(d2*nparm1, double); /* d^2 x (npar-1) */
  Prod          = (double *)Calloc(d2*nparm1, double); /* d^2 x (npar-1) */
  IxBpr         = (double *)Calloc(       d4, double); /* d^2 x d^2      */
  BprxI         = (double *)Calloc(       d4, double); /* d^2 x d^2      */
  IxLi          = (double *)Calloc(       d4, double); /* d^2 x d^2      */
  IxSplLi       = (double *)Calloc(       d4, double); /* d^2 x d^2      */
  FactorI       = (double *)Calloc(       d4, double); /* d^2 x d^2      */

  /* initialize some matrices to zero  */
  for(j=0;j<d2;j++) {
    *(Lambdahlf +j) = 0.0;
    *(Lambda + j) = 0.0;
  }

  /* calculate nu, the d.f. parameter */
  eth0 = exp(*theta);
  nu = (2.0*xd+2.0)*(eth0+1.0);

  /* calculate Lambdahlf, the chol. sq. root of Lambda, the matrix parameter */
  l=0;
  for(j=0;j<d;j++){
    *(Lambdahlf + d*j + j) = exp(*(theta + j + 1));
    for(k=(j+1);k<d;k++) {
      *(Lambdahlf + d*k + j) = *(theta + d + 1 + l);
      l++;
    }
  }

  /* printmat(Lambdahlf, d, d, "Lambdahlf"); */

  /* calculate Lambda = t(Lambdahlf) %*% Lambdahlf */
  for(j=0;j<d;j++)
    for(k=0;k<d;k++)
      for(h=0;h<d;h++) *(Lambda + d*k + j) += *(Lambdahlf + d*j + h) *(*(Lambdahlf + d*k + h));

  /* printmat(Lambda, d, d, "Lambda"); */

  detL = det(Lambda, xd2buff, pd);

  /* Rprintf("detL = %g\n",detL); */

  /* zero the gradient vector */
  for(j=0;j<npar;j++) *(G+j) = 0.0;

  /* preparing to take derivatives w.r.t. theta_j, j=2,...,(d*(d+1)/2 + 1)              */
  /* which are the parameters of the Lambda matrix.  Lambda = Lambdahlf' Lambdahlf,     */
  /* where Lambdahlf is the (upper triangular) cholesky square-root of Lambda.          */
  /* Lambdahlf is parameterized with diagonal elements exp(theta_j), j=2,...,(d+1) and  */
  /* off diagonal elements theta_j, j=(d+2),...,(d*(d+1)/2 + 1)                         */
  /* The names of the variables below reflect the fact that in my notes, Lambdahlf      */
  /* is a.k.a.  B.                                                                      */
  /*                                                                                    */
  /* First, compute derivatives of B (the cholesky square root of Lambda),              */
  /* as this independent of gene, and thus can be done before the loop over             */
  /* genes.                                                                             */
  for(j=0;j<d2;j++) *(xd2buff +j) = *(Lambda +j);
  matinv(xd2buff, Lambdainv, pd);
  for(j=0;j<d2;j++) {
    for(k=0;k<nparm1;k++) {
      *(dvecBdtheta + d2*k + j) = 0.0;
      *(dvecBprdtheta + d2*k + j) = 0.0;
      *(FactorII + d2*k + j) = 0.0;
      *(Prod + d2*k + j) = 0.0;
    }
  }
  for(j=0;j<d4;j++) {
    *(IxBpr + j)= 0.0;
    *(BprxI +j) = 0.0;
    *(IxLi+j) = 0.0;
  }

  for(i1=0;i1<d;i1++)
    for(i2=0;i2<d;i2++)
      for(j2=0;j2<d;j2++){
        *(IxBpr + d2*(d*i1+j2) + d*i1+i2) = *(Lambdahlf + d*i2+j2);
                                             /* I x Lambda^(-1) times the coefficients */
	*(IxLi  + d2*(d*i1+j2) + d*i1+i2) = *(Lambdainv + d*i2+j2) * (nu -xd - 1.0)/2.0;
	*(BprxI   + d2*(d*j2+i1) + d*i2+i1) = *(Lambdahlf + d*i2+j2);
      }

  /* printmat(IxBpr, d2, d2, "IxBpr");  */
  /* printmat(IxLi,  d2, d2, "IxLi");   */
  /* printmat(BprxI,   d2, d2, "BxI");  */

  l = 0;
  for(j=0;j<d;j++){
       /* derivatives of   B w.r.t. theta_1, ..., theta_d */
    *(dvecBdtheta + d2*j + d*j + j) =  *(Lambdahlf + d*j + j);
      /* derivatives of Bpr w.r.t. theta_1, ..., theta_d */
    *(dvecBprdtheta + d2*j + d*j + j) =  *(Lambdahlf + d*j + j);

    for(k=(j+1);k<d;k++){
        /* derivatives of   B w.r.t. theta_(d+1), ..., theta_p */
      *(dvecBdtheta + d2*(d + l) + d*k + j) = 1.0;
        /* derivatives of Bpr w.r.t. theta_(d+1), ..., theta_p */
      *(dvecBprdtheta + d2*(d + l) + d*j + k) = 1.0;
      l++;
    }
  }

  /* printmat(dvecBdtheta, d2, pm1, "dvecBdtheta");     */
  /* printmat(dvecBprdtheta, d2, pm1, "dvecBprdtheta"); */

  for(j=0;j<d2;j++)
    for(k=0;k<nparm1;k++)
      for(i1=0;i1<d2;i1++)
        *(FactorII + d2*k + j) += *(IxBpr + d2*i1 + j) * (*(dvecBdtheta   + d2*k + i1)) +
                                  *(BprxI   + d2*i1 + j) * (*(dvecBprdtheta + d2*k + i1));

  /* printmat(FactorII, d2, pm1, "FactorII"); */

  /* Loop over genes */
  for(h=0;h<N;h++){
    /* derivative w.r.t. theta_1, where nu = exp(theta_1) */
    xnreps = (double) (*(nreps +h));
    for (j=0;j<d2;j++) *(xd2buff+j) = 0.0;
    detS = pow(xnreps - 1.0, xd) * det((MVM+d2*h), xd2buff, pd);
    for (j=0;j<d2;j++) {
      *(xd2buff+j) = 0.0;
      *(SplL+j) = 0.0;
    }
    for(j=0;j<d2;j++) *(SplL +j) = (xnreps-1.0)*(*(MVM +d2*h +j)) +(*(Lambda +j));
    detSplL = det(SplL, xd2buff, pd);
    v = detL/detSplL;
    *G += -1.0 * (dimgamma((nu + xnreps - xd -2.0)/2.0, d) - dimgamma((nu - xd - 1.0)/2.0, d) +
           log(v))*eth0*(xd + 1.0);

    for(j=0;j<d2;j++) *(xd2buff +j) = *(SplL +j);
    matinv(xd2buff, SplLinv, pd);

    for(j=0;j<d4;j++) {
      *(IxSplLi + j) = 0.0;
      *(FactorI + j) = 0.0;
    }

    for(j=0;j<(d2*nparm1);j++)
      *(Prod + j) = 0.0;

    for(i1=0;i1<d;i1++)
      for(i2=0;i2<d;i2++)
        for(j2=0;j2<d;j2++)
                 /* I x ((n-1) S + Lambda)^(-1) times the coefficients */
	  *(IxSplLi  + d2*(d*i1+j2) + d*i1+i2) = *(SplLinv + d*j2+i2) * (nu + xnreps - xd - 2.0)/2.0 ;

    /* if (h<=5) printmat(IxSplLi, d2, d2, "IxSplLi"); */

    for (j=0;j<d4;j++)
      *(FactorI + j) = *(IxLi + j) - *(IxSplLi + j);

    /* if (h<=5) printmat(FactorI, d2, d2, "FactorI"); */

     for(j=0;j<d2;j++)
       for(k=0;k<nparm1;k++)
	for(i1=0;i1<d2;i1++)
	  *(Prod + d2*k + j) += *(FactorI + d2*i1 + j) * (*(FactorII   + d2*k + i1));

     /* if (h<=5) printmat(Prod, d2, nparm1, "Prod"); */

     for(k=0;k<nparm1;k++)
       for(j=0;j<d;j++)
         *(G + k + 1) += *(Prod + d2*k + (d*j + j)) * (-1.0);
  }

  Free(Lambdahlf);
  Free(Lambda);
  Free(Lambdainv);
  Free(xd2buff);
  Free(SplL);
  Free(SplLinv);
  Free(dvecBdtheta);
  Free(dvecBprdtheta);
  Free(FactorII);
  Free(Prod);
  Free(IxBpr);
  Free(BprxI);
  Free(IxLi);
  Free(IxSplLi);
  Free(FactorI);

}

double lmgamma(double a, int m)
{
  int i;
  double pi=3.141592653589793, s=0.0, xm, xi;

  xm = (double) m;

  for(i=0;i<m;i++){
    xi = (double) i;
    s += lgamma(a - 0.5 * xi);
  }
  s += xm*(xm-1.0)*log(pi)/4.0;
  return s;
}

double dimgamma(double a, int m)
{
  int i;
  double s=0.0, xi;

  for(i=0;i<m;i++){
    xi = (double) i;
    s += digamma(a - 0.5 * xi);
  }
  return s;
}

double trimgamma(double a, int m)
{
  int i;
  double s=0.0, xi;

  for(i=0;i<m;i++){
    xi = (double) i;
    s += trigamma(a - 0.5 * xi);
  }
  return s;
}

double det(double *x, double *xd2buff, int *pd)
{
        int d, i;
	double tmp, ans;

        d = *pd;
        ans = 0.0;
	if (d==1) ans = *x;
	if (d>1){
          tmp = 1.0;
          chol(x, xd2buff, pd);
          for (i=0;i<d;i++) tmp = tmp*(*(xd2buff+d*i+i));
            ans = tmp*tmp;
        }
        ans = (R_FINITE(ans) ? ans : MDXMIN);
	return ans;
}

void tdet(double *x, double *xd2buff, int *pd, double *pans)
{
	*pans = det(x, xd2buff, pd);
}

void tloglik(double *ptheta, double *MVM, int *pN, int *pd, 
             int *pnreps, double *pans)
{
  Data *y;
  int inpar;
  int d, npar;
  y = (Data *)Calloc(1, Data);

  y->MVM = MVM;
  y->pN = pN;
  y->pd = pd;
  y->nreps = pnreps;

  d = *pd;
  npar = d*(d+1)/2 + 1;
  inpar = (int) npar;

  *pans = loglik(inpar, ptheta, y);

  Free(y);

}

void tGloglik(double *ptheta, double *MVM, int *pN, int *pd, 
              int *pnreps, double *pG)
{
  Data *y;
  int d, npar;
  int inpar;

  y = (Data *)Calloc(1, Data);

  y->MVM = MVM;
  y->pN = pN;
  y->pd = pd;
  y->nreps = pnreps;
  d = *pd;
  npar = d*(d+1)/2 + 1;
  inpar = (int) npar;

  Gloglik(inpar, ptheta, pG, y);

  Free(y);

}

void fHESS(double *x, Data *y, double *G, double *H, optimgr *grad)
{
  int i,j,d,npar;
  double h,temp,*G1;
  int inpar;

  d = *(y->pd);
  npar = d*(d+1)/2 + 1;
  inpar = (int) npar;
  G1 = (double *)Calloc(npar, double);
  for (j=0;j<npar;j++) {
    temp=*(x+j);
    h=EPS*fabs(temp);
    if (h == 0.0) h=EPS;
    /* Trick to reduce finite precision error. */
    *(x+j)=temp+h;
    h=*(x+j)-temp;
    (*grad)(inpar, x, G1, y);
    *(x+j)=temp;
    /* Forward difference formula. */
    for (i=0;i<(j+1);i++) {
      *(H + npar*j + i) = (*(G1+i)-(*(G+i)))/h;
      *(H + npar*i + j) = *(H + npar*j + i);
    }
  }

  Free(G1);

}

void printmat(double *pA, int nr, int nc, char *name)
{
  int j, k;

  Rprintf("%s = \n", name);
  for (j=0;j<nr;j++) {
    for (k=0;k<nc;k++) Rprintf("%g ",*(pA + nr*k + j));
    Rprintf("\n");
  }
}

