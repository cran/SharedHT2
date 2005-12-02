#include<R.h>
#include<Rmath.h>

void chol(double *s, double *t, int *pd);

void rwishart1(double *pdf, int *pd, double *pSqrtSigma, double *pW);
void matinv(double *a, double *yvv, int *pm);
void rnormn(int *pn, double *ans);

void SimOneMVN_mxIW(double *nu, double *Lbdinvhlf, double *f1f2, int *pd, 
                    int *pnreps, int *pN, double *es, double *YY)
{
  int i, j, k, l, d, d2, N, nreps, mxnreps, Jrand;
  int *lbuff;

  double xd, sm, tstnu, nu_1, nu_2, zz, lambda, f_1, f_2;

  double *df, *pW, *SgmHlf, *Y, *xbuff, *Sigma, *sig, *SigInv;

  N = *pN;
  d = *pd;
  xd = (double) d;
  d2 = d*d;

  mxnreps=0;
  for(l=0;l<N;l++) if(mxnreps < *(pnreps+l)) mxnreps = *(pnreps+l);

  lbuff         = (int   *)S_alloc(        1,sizeof(int));

  df            = (double *)S_alloc(        1, sizeof(double));
  pW            = (double *)S_alloc(       d2, sizeof(double));
  xbuff         = (double *)S_alloc(        d, sizeof(double));
  SgmHlf        = (double *)S_alloc(       d2, sizeof(double));
  Y             = (double *)S_alloc(mxnreps*d, sizeof(double));
  Sigma         = (double *)S_alloc(       d2, sizeof(double));
  SigInv        = (double *)S_alloc(       d2, sizeof(double));
  sig           = (double *)S_alloc(        d, sizeof(double));

  f_1 = *f1f2;
  f_2 = *(f1f2+1);
  lambda = *nu/(2 * xd + 2.0);

  nu_1 = (2.0 * xd + 2.0)*(1.0 + (lambda - 1.0) *(1.0 - f_1)/(1-f_1/f_2));
  nu_2 = (2.0 * xd + 2.0)*(1.0 + (lambda - 1.0) *          f_2          );

  tstnu = f_1/(nu_2 - (2.0*xd+2.0)) + (1.0-f_1)/(nu_1 - (2.0*xd+2.0));
  tstnu = 1.0/tstnu + 2.0*xd + 2.0;
  Rprintf("nu_1=%g, nu_2=%g, nu ?= %g", nu_1, nu_2, tstnu);

  GetRNGstate();

  /* NOTE:                                                                              */
  /* this block computes the average std dev over genes from the model                  */
  /* its diagonal elements, passed to the pointer, sig (of size 3)                      */
  /* are used for the purposes of assigning mean value to Y's under the alternative     */
  /*                                                                                    */
  for(i=0;i<d;i++)                                                   
    for(j=0;j<d;j++){
      sm = 0.0;
      for(k=0;k<d;k++) 
        sm += *(Lbdinvhlf + d*i + k) * (*(Lbdinvhlf + d*j + k));
      *(SigInv + d*j + i) = sm;
    }
  matinv(SigInv, Sigma, pd);
  for(i=0;i<d;i++) *(sig + i) = pow((*(Sigma + d*i + i))/(*nu - 2.0*xd - 2.0), 0.5);

  for(l=0;l<N;l++){  

    /* Pick J = 1 w.p. 1-f_1, J= 2 w.p. f_1.                                            */
    /* Then draw an InvWish_d(nu_J, Lambda) matrix.  This is done                       */
    /* using the result:  if Sigma^(-1) ~ Wish_d(nu-d-1, Lambda^(-1)) then              */
    /* Sigma ~ InvWish_d(nu, Lambda).  I simulate N i.i.d. Wish_d(nu-d-1,Lambda^(-1))   */
    /* matrices and then invert to get Sigma.  One more catch, my rwishart routine      */
    /* uses the cholesky square root of the parameter matrix, Lambda instead of Lambda  */
    /* itself.  Since I want the parameter matrix in the Wishart to be Lambda^(-1) then */
    /* I should pass its cholesky square root which is Lbdinvhlf, e.g. the cholesky     */
    /* square root of Lambda inverse.  That is computed in the calling R script and     */
    /* passed in.  Notice the need to check that Lambda is nonsingular and that         */
    /* nu > 2*d + 2 (required so that the expected value of the inverse wishart         */
    /* is finite.)                                                                      */
    /*                                                                                  */
    zz = unif_rand();
    Jrand = 1*(zz < f_1);
    *df = (1-Jrand)*nu_1 + Jrand*nu_2 - xd - 1.0;
    rwishart1(df, pd, Lbdinvhlf, pW);
    matinv(pW, Sigma, pd);
    /*                                                                                  */
    /* Sigma ~ (1-f_1)*InvWish_d(nu_1, Lambda) + f_1*InvWish_d(nu_2, Lambda)            */
    /*                                                                                  */
    /* Next, use Sigma to simulate i.i.d. N(0_d, Sigma)'s                               */
    nreps = *(pnreps + l);
    *lbuff = nreps*d;
    rnormn(lbuff, Y); 
    chol(Sigma, SgmHlf, pd);

    for(i=0;i<nreps;i++){
      for(j=0;j<d;j++){
        sm = 0.0;
        for(k=0;k<d;k++) sm += (*(SgmHlf +d*j +k))*(*(Y +d*i +k));
        *(xbuff+j) = sm;
      }
      for(j=0;j<d;j++) *(Y + d*i + j) = *(xbuff + j) + *(es + l)*(*(sig + j));
    }
    for(i=0;i<(nreps*d);i++) *(YY + mxnreps*d*l + i) = *(Y+i);
  }
  PutRNGstate();

}
