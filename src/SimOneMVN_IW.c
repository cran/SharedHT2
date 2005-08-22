#include<R.h>
#include<Rmath.h>

void rwishart1(double *pdf, long *pd, double *pSqrtSigma, double *pW);
void matinv(double *a, double *yvv, long *pm);
void rnormn(long *pn, double *ans);

void printmat(double *pA, long nr, long nc, char *name);
double chol(double *, double *, long *);

void SimOneMVN_IW(double *nu, double *Lbdinvhlf, long *pd, long *pnreps,
               long *pN, double *es, double *YY)
{
  long i, j, k, l, d, npar, npar2, d2, N, nreps, mxnreps;
  long *lbuff;

  double mu, xd, sm;

  double *df, *pW, *SgmHlf, *xbuff, *Y;
  double *Sigma, *LbdHlf, *sig, *SigInv;

  N = *pN;
  d = *pd;
  xd = (double) d;
  d2 = d*d;

  mxnreps=0;
  for(l=0;l<N;l++) if(mxnreps < *(pnreps+l)) mxnreps = *(pnreps+l);

  lbuff       = (long   *)S_alloc(      1,sizeof(long));

  df          = (double *)S_alloc(        1, sizeof(double));
  pW          = (double *)S_alloc(       d2, sizeof(double));
  xbuff       = (double *)S_alloc(        d, sizeof(double));
  SgmHlf      = (double *)S_alloc(       d2, sizeof(double));
  Y           = (double *)S_alloc(mxnreps*d, sizeof(double));
  Sigma       = (double *)S_alloc(       d2, sizeof(double));
  SigInv      = (double *)S_alloc(       d2, sizeof(double));
  LbdHlf      = (double *)S_alloc(       d2, sizeof(double));
  sig         = (double *)S_alloc(        d, sizeof(double));

  /* NOTE:                                                                           */
  /* this block computes the average std dev over genes from the model               */
  /* its diagonal elements, passed to the pointer, sig (of size 3)                   */
  /* are used for the purposes of assigning mean value to Y's under the alternative  */
  /*                                                                                 */
  for(i=0;i<d;i++)
    for(j=0;j<d;j++){
      sm = 0.0;
      for(k=0;k<d;k++)
        sm += *(Lbdinvhlf + d*i + k) * (*(Lbdinvhlf + d*j + k));
      *(SigInv + d*j + i) = sm;
    }
  matinv(SigInv, Sigma, pd);
  for(i=0;i<d;i++) *(sig + i) = pow((*(Sigma + d*i + i))/(*nu - 2.0*xd - 2.0),0.5);

  GetRNGstate();

  for(l=0;l<N;l++){
    nreps = *(pnreps +l);
    *lbuff = nreps * d;
    *df = *nu - xd - 1.0;

    /*
    /* First an InvWish_d(nu, Lambda) matrix.  This is done                             */
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
    rwishart1(df, pd, Lbdinvhlf, pW);
    matinv(pW, Sigma, pd);
    /*                                                                                  */
    /* Sigma ~ InvWish_d(nu, Lambda)                                                    */
    /*                                                                                  */

    /* Next, use Sigma to simulate Y ~ i.i.d. N(0_d, Sigma)                             */
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
