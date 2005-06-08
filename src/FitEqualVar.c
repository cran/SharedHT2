#include<R.h>
#include<Rmath.h>
#define EPS 1.0e-7

typedef struct{
  double *S;
  long *pN;
  long *pd;
  long *nreps;
} DataEV;

typedef double optimfn(int n, double *par, void *ex);
typedef void optimgr(int n, double *par, double *gr, void *ex);

void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
           int *fail, double abstol, double intol, void *ex,
           double alpha, double beta, double gamma, int trace,
           int *fncount, int maxit);

void vmmin(int n, double *x, double *Fmin,
           optimfn fn, optimgr gr, int maxit, int trace,
           int *mask, double abstol, double reltol, int nREPORT,
           void *ex, int *fncount, int *grcount, int *fail);

void fHESSEV(double *x, DataEV *y, double *G, double *H, optimgr *grad);

optimfn loglikEV;
optimgr GloglikEV, *gradEV;

void FitEqualVar(double *ptheta0, double *S, long *pN, long *pd, long *pnreps,
                 long *pverbose, double *objval, double *estimate, long *fail,
                 long *fncount, long *grcount, long *mask, long *usegr,
                 double *G, double *H)
{
  int inpar=2, *ifail, *ifncount, *imask, verb;
  DataEV *y;

  verb = (int) *pverbose;
  ifail = (int *) fail;
  ifncount = (int *) fncount;
  imask = (int *) mask;

  y = (DataEV *) Calloc(1, DataEV);

  y->S = S;
  y->pN = pN;
  y->pd = pd;
  y->nreps = pnreps;

  *estimate = 0.0;
  *(estimate+1) = 0.0;

  if(*usegr==0){
    nmmin(inpar, ptheta0, estimate, objval, loglikEV, ifail, -1e200, 1e-8, y,
          1.0, 0.5, 2.0, verb, ifncount, 10000);
  }
  else{
    vmmin(inpar, ptheta0, objval, loglikEV, GloglikEV, 10000, verb, imask, -1e200,
          1e-8, 10, y, ifncount, ifncount, ifail);
    *estimate = *ptheta0;
    *(estimate+1) = *(ptheta0+1);
  }

  GloglikEV(inpar, estimate, G, y);
  gradEV = &GloglikEV;
  fHESSEV(estimate, y, G, H, gradEV);
  Free(y);
}

double loglikEV(int inpar, double *theta, void *yy)
{
  long N, d, h;
  long *nreps;
  double xd, loglik, logC, s, r, xn, xn1_o2, xn2_o2, q, q_o2r;
  double *S;
  DataEV *y;

  y = (DataEV *) yy;
  S = y->S;
  N = *(y->pN);
  d = *(y->pd);
  nreps = y->nreps;
  xd = (double) d;

  s = exp(*theta);
  r = exp(*(theta+1));
  xn2_o2 = s;
  loglik = 0.0;
  for(h=0;h<N;h++)
  {
    xn = (double)(*(nreps+h));
    xn1_o2 = xd*(xn-1.0)/2.0;
    logC = lgamma(xn1_o2+xn2_o2) - lgamma(xn1_o2) -
           lgamma(xn2_o2) - log(2.0*r);
    q = *(S+h);
    q_o2r = q/(2.0*r);
    loglik += logC + (xn1_o2-1.0)*log(q_o2r) - (xn1_o2 + xn2_o2)*log(1.0+q_o2r);
  }
  return(loglik*(-1.0));
}

void GloglikEV(int inpar, double *theta, double *G, void *yy)
{
  long N, d, h;
  long *nreps;
  double xd, g1, g2, s, r, xn, xn1_o2, xn2_o2, q, q_o2r;
  double *S;
  DataEV *y;

  y = (DataEV *) yy;
  S = y->S;
  N = *(y->pN);
  d = *(y->pd);
  nreps = y->nreps;
  xd = (double) d;

  s = exp(*theta);
  r = exp(*(theta+1));
  xn2_o2 = s;
  g1 = 0.0;
  g2 = 0.0;
  for(h=0;h<N;h++)
  {
    xn = (double)(*(nreps+h));
    xn1_o2 = xd*(xn-1.0)/2.0;
    q = *(S+h);
    q_o2r = q/(2.0*r);
    g1 += digamma(xn1_o2+xn2_o2) - digamma(xn2_o2) - log(1.0+q_o2r);
    g2 += (xn1_o2+xn2_o2)*q/(q*r+2.0*r*r) - xn1_o2/r;
  }
  *G = g1*s*(-1.0);
  *(G+1) = g2*r*(-1.0);
}

void fHESSEV(double *x, DataEV *y, double *G, double *H, optimgr *grad)
{
  long i,j,d,npar;
  double h,temp,*G1;
  int inpar=2;
  npar = (long)inpar;

  d = *y->pd;
  G1 = (double *)Calloc(npar, double);
  for (j=0;j<npar;j++) {
    temp = (*(x+j));
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
