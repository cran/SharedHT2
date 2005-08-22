/* SimOneNorm_IG simulates data directly off the Normal/Inverse Gamma                  */
/* model (see equation 1 below).  It was designed to compare the ShHT2           */
/* statistic's performance with several other statistics, the HT2, the           */
/* ShUT2, and the UT2 (see manuscript).                                          */
/*                                                                               */
/*   [ Y | sigma2 ] = N(0, sigma2)                                               */
/*                                                                          (1)  */
/*   [ sigma2 ]     = InvGamma(shape, rate)                                      */
/*                                                                               */
/*  except that the argument 'es' is used to pass an effect size giving          */
/*  the corresponding row a mean value as to reject the simple 'UT2' f-test      */
/*  at alternative given by ncp = es                                             */
/*                                                                               */
#include<stdio.h>
#include<string.h>
#include<R.h>
#include<Rmath.h>

#define EPS 1.0e-7
#define vabsmax(v, n, l, a, sgn) a=0.0;                             \
                                 for(l=0;l<n;l++)                   \
                                   if(fabs(*(v+l))>a){              \
                                     a=fabs(*(v+l));                \
                                     sgn = (*(v+l)>0 ? 1.0 : -1.0); \
                                   }                                \
                                 0

void rnormn(long *pn, double *ans);

void printmat(double *pA, long nr, long nc, char *name);

void SimOneNorm_IG(double *shape, double *rate, long *pd, long *pnreps,
                   long *pN, double *es, double *YY)
{
  long i, j, l, d, N, nreps, mxnreps;
  long *lbuff;

  double sig, sigma2, sigma, xd;
  double Top, stat1, stat2, stat3, stat4, pval1, pval2, pval3, pval4;

  double *xbuff, *Y;

  N = *pN;
  d = *pd;

  mxnreps=0;
  for(l=0;l<N;l++) if(mxnreps < *(pnreps+l)) mxnreps = *(pnreps+l);

  lbuff       = (long   *)S_alloc(        1, sizeof(long));
  xbuff       = (double *)S_alloc(        d, sizeof(double));
  Y           = (double *)S_alloc(mxnreps*d, sizeof(double));

  GetRNGstate();

  /* NOTE:                                                                             */
  /* this block computes the average std dev over genes from the model                 */
  /* it is used for the purposes of assigning mean value to Y's under the alternative  */
  /*                                                                                   */

  sig = pow(*rate/(*shape-1.0), 0.5);

  for(l=0;l<N;l++){  

    /*                                                                                  */
    /* First, simulate sigma2 ~ InvGamma(shape, rate).  This is done                    */
    /* using the result:  if sigma2^(-1) ~ Gamma(shape, rate) then                      */
    /* sigma2 ~ InvGamma(shape, rate).                                                  */

    sigma2 = 1.0/rgamma(*shape, 1.0/(*rate));  

    /*                                                                                   */
    /* sigma2 ~ InvGamma(shape, rate)                                                    */
    /*                                                                                   */
    /* Next, use sigma2 to simulate Y ~ i.i.d. N(0_d, sigma2)                            */

    nreps = *(pnreps+l);
    *lbuff = nreps*d;
    rnormn(lbuff, Y); 

    sigma = pow(sigma2, 0.5);

    for(i=0;i<nreps;i++){
	for(j=0;j<d;j++) *(Y + d*i + j) = *(Y + d*i + j)*(sigma) + *(es + l)*(sig);
    }
    for(i=0;i<(nreps*d);i++) *(YY + mxnreps*d*l + i) = *(Y+i);
  }
  PutRNGstate();

}
