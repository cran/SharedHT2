/* Simnu_mix, an attempt to explore the robustness of the ShHT2                        */
/* statistic to model assumptions. It is (luckily) a minor modification of             */
/* SimW_IW.                                                                            */
/*                                                                                     */
/* Recall SimW_IW simulates data directly off the Wishart/Inverse Wishart              */
/* model:                                                                              */
/*                                                                                     */
/*   [ Y_i | Sigma ] = i.i.d MVN(0, Sigma)                                             */
/*                                                                               (1)   */
/*   [ Sigma ]     = InvWish(nu, Lambda)                                               */
/*                                                                                     */
/* except that the argument 'es' is used to pass an effect size giving                 */
/* the corresponding row a mean value as to reject the simple 'UT2' f-test             */
/* at alternative given by ncp = es                                                    */
/*                                                                                     */
/* Simnu_mix simulates data under the following model, (2), which is designed          */
/* to have mean zero and va-/cova-riance matrix identical to that of the               */
/* original model (1) in order that a given effect size may be directly                */
/* compared in the two models.  This is done as follows                                */
/*                                                                                     */
/*                                                                                     */
/*   [ Y_i | Sigma ] = i.i.d. MVN(0_d, Sigma)                                          */
/*                                                                               (2)   */
/*     [ Sigma ]     = (1-f_1) * InvWish(nu_1, Lambda) + f_1 * InvWish(nu_2, Lambda)   */
/*                                                                                     */
/*                                                                                     */
/* where f_1, nu_1, and nu_2 are chosen so that var[Y] is equal under the two models   */
/* (1) and (2).  This is done, for given nu and Lambda specified in (1),               */
/* with a particular choice of 0 < f_1 < 1 and 0 < f_2 < 1, with f_1/f_2 < 1 by setting*/
/*                                                                                     */
/*  nu_1 = (2d+2)*(1+ (lambda-1)*(1-f_1)/(1-f_1/f_2))                                  */
/*  nu_2 = (2d+2)*(1+ (lambda-1)*f_2)                                                  */
/*                                                                                     */
/*  where lambda = nu/(2d+2).  Notice that nu = (2d+2)*(exp(theta) + 1) > 2d+2         */
/*  so that lambda > 1.  Thus 'f_2'  specifies the portion of the fractional part      */
/*  of lambda by which nu_2 exceeds 1.  Notice also in model (2) that the first        */
/*  term results in smaller variances while the second term results in larger          */
/*  variances.                                                                         */
/*                                                                                     */
/*                                                                                     */
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

#define MIN(x,y) (x<=y ? x : y)

typedef struct{
  double *MVM;
  long *pN;
  long *pd;
  long *nreps;
} Data;

typedef struct{
  double *S;
  long *pN;
  long *pd;
  long *nreps;
} DataEV;

typedef struct{
  long id;
  double mean;
  double ShHT2;
  double ShHT2pval;
  double HT2;
  double HT2pval;
  double ShUT2;
  double ShUT2pval;
  double UT2;
  double UT2pval;
} gene;

typedef double optimfn(int n, double *par, void *ex);
typedef void optimgr(int n, double *par, double *gr, void *ex);
typedef int CmprFun(const void *x, const void *y);

CmprFun cmprShHT2, cmprHT2, cmprShUT2, cmprUT2, *cmpr;

void tloglik(long *pnpar, double *ptheta, double *MVM, long *pN, long *pd, 
	     long *pnreps, double *pans);

void tGloglik(double *ptheta, double *MVM, long *pN, long *pd, long *pnreps, 
	      double *pG);

void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
           int *fail, double abstol, double intol, void *ex,
           double alpha, double beta, double gamma, int trace,
           int *fncnt, int maxit);

void vmmin(int n, double *x, double *Fmin,
           optimfn fn, optimgr gr, int maxit, int trace,
           int *mask, double abstol, double reltol, int nREPORT,
           void *ex, int *fncnt, int *grcnt, int *fail);

void fHESS(double *x, Data *y, double *G, double *H, 
	   void (*grad)(int p, double *x, double *G, Data *y));

void rwishart1(double *pdf, long *pd, double *pSqrtSigma, double *pW);
void matinv(double *a, double *yvv, long *pm);
void rnormn(long *pn, double *ans);

void printmat(double *pA, long nr, long nc, char *name);

void FitInvWish1(double *ptheta0, long *pverbose, Data *y, double *objval, 
		 double *estimate, double *estimater, long *fail, long *fncnt, 
		 long *grcnt, long *mask, long *usegr, double *G, double *H);

void FitEqualVar1(double *ptheta0, long *pverbose, DataEV *y, double *objval, 
                  double *estimate, long *fail, long *fncnt, long *grcnt, 
                  long *mask, long *usegr, double *G, double *H);

void printglist(gene *x, long N, char *strng);

void Simnu_mix(long *verb, long *fail, long *fncnt, long *grcnt, long *mask, 
               long *usegr, long *pnsim, double *nu, double *Lbdinvhlf, long *pd, 
               long *pnreps, long *pN, double *es, double *f1f2, double *coef, 
               double *coefEV, double *FDRlist, long *pnFDRlist, 
               double *fdrtbl, double *roctbl)
{
  long i, j, k, l, d, npar, npar2, d2, d4, N, nreps, mxnreps, nsim, isim, indx;
  long ntruepos, nFDRlist, flagsig, Nsig, nTP, nFP, Jrand;
  long *lbuff, *pnpar;
  char *itfnm;

  double mu, xnreps, xN, xd, sm, smEV, xn2, nu_isim, r, s, vamx, sgn, xl, tstnu, nu_1, nu_2;
  double Top, stat1, stat2, stat3, stat4, pval1, pval2, pval3, pval4, zz, lambda, f_1, f_2;

  double *df, *pW, *SgmHlf, *xbuff, *muhat, *res, *Sighat, *WSSQ, *ptheta0, *rFDR;
  double *objv, *estimate, *estimater, *G, *H, *Y, *Sigma, *Lambda_isim, *LbdHlf, *sig;
  double *SigInv;

  Data *y;
  DataEV *yEV;
  gene *genelist;
  FILE *itfnm_ptr, *dtfl_ptr;
  char *ch;
  char *alnu;

  nsim = *pnsim;
  N = *pN;
  xN = (double)N;
  d = *pd;
  xd = (double) d;
  d2 = d*d;
  d4 = d2*d2;
  npar = d*(d+1)/2 + 1;
  npar2 = npar*npar;
  nFDRlist = *pnFDRlist;

  mxnreps=0;
  for(l=0;l<N;l++) if(mxnreps < *(pnreps+l)) mxnreps = *(pnreps+l);

  lbuff         = (long   *)S_alloc(        1,sizeof(long));
  pnpar         = (long   *)S_alloc(        1,   sizeof(long));

  df            = (double *)S_alloc(        1, sizeof(double));
  pW            = (double *)S_alloc(       d2, sizeof(double));
  xbuff         = (double *)S_alloc(        d, sizeof(double));
  SgmHlf        = (double *)S_alloc(       d2, sizeof(double));
  muhat         = (double *)S_alloc(      N*d, sizeof(double));
  res           = (double *)S_alloc(       d2, sizeof(double));
  Sighat        = (double *)S_alloc(     N*d2, sizeof(double));
  WSSQ          = (double *)S_alloc(        N, sizeof(double));
  ptheta0       = (double *)S_alloc(     npar, sizeof(double));
  objv          = (double *)S_alloc(        1, sizeof(double));
  estimate      = (double *)S_alloc(     npar, sizeof(double));
  estimater     = (double *)S_alloc(     npar, sizeof(double));
  G             = (double *)S_alloc(     npar, sizeof(double));
  H             = (double *)S_alloc(    npar2, sizeof(double));
  Y             = (double *)S_alloc(mxnreps*d, sizeof(double));
  Sigma         = (double *)S_alloc(       d2, sizeof(double));
  SigInv        = (double *)S_alloc(       d2, sizeof(double));
  Lambda_isim   = (double *)S_alloc(       d2, sizeof(double));
  LbdHlf        = (double *)S_alloc(       d2, sizeof(double));
  sig           = (double *)S_alloc(        d, sizeof(double));
  rFDR          = (double *)S_alloc(        N, sizeof(double));
  itfnm         = (char   *)S_alloc(      100, sizeof(char));
  alnu          = (char   *)S_alloc(       62, sizeof(char));
  ch            = (char   *)S_alloc(        2, sizeof(char));
  y             = (Data   *)S_alloc(        1, sizeof(Data));
  yEV           = (DataEV *)S_alloc(        1, sizeof(DataEV));
  genelist      = (gene   *)S_alloc(        N, sizeof(gene));

  alnu = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
  *(ch+1) = '\0';

  *pnpar = npar;

  f_1 = *f1f2;
  f_2 = *(f1f2+1);
  lambda = *nu/(2 * xd + 2.0);

  nu_1 = (2.0 * xd + 2.0)*(1.0 + (lambda - 1.0) *(1.0 - f_1)/(1-f_1/f_2));
  nu_2 = (2.0 * xd + 2.0)*(1.0 + (lambda - 1.0) *          f_2          );

  tstnu = f_1/(nu_2 - (2.0*xd+2.0)) + (1.0-f_1)/(nu_1 - (2.0*xd+2.0));
  tstnu = 1.0/tstnu + 2.0*xd + 2.0;
  Rprintf("nu_1=%g, nu_2=%g, nu ?= %g", nu_1, nu_2, tstnu);

  for (l=0;l<8*nFDRlist;l++) *(fdrtbl + l) = 0.0;

  ntruepos = 0;
  for (l=0;l<N;l++){
    xl = (double) (l+1);
    *(rFDR+l) = xl/xN;
    for(j=0;j<8;j++) *(roctbl + N*j + l) = 0.0;
    ntruepos += 1*(fabs(*(es+l)) > 0.001);
  }

  GetRNGstate();

  /*  name the 'iterno' file uniquely  */
  strcat(itfnm, "iterno");
  for(i=0;i<8;i++) {
    indx = (long) (62.0*unif_rand());
    *ch = *(alnu + indx);
    strncat(itfnm, ch, 1);
  }

  itfnm_ptr = fopen(itfnm, "w+");

  /* NOTE:                                                                            */
  /* this block computes the average std dev over genes from the model                */
  /* its diagonal elements, passed to the pointer, sig (of size 3)                    */
  /* are used for the purposes of assigning mean value to Y's under the alternative   */
  /*                                                                                  */
  for(i=0;i<d;i++)
    for(j=0;j<d;j++){
      sm = 0.0;
      for(k=0;k<d;k++) 
        sm += *(Lbdinvhlf + d*i + k) * (*(Lbdinvhlf + d*j + k));
      *(SigInv + d*j + i) = sm;
    }
  matinv(SigInv, Sigma, pd);
  for(i=0;i<d;i++) *(sig + i) = pow((*(Sigma + d*i + i))/(*nu - 2.0*xd - 2.0), 0.5);

  /*start of simulation loop: (no simulation loop just yet...just try one rep first   */
  for(isim=0;isim<nsim;isim++){

    if(isim>0){
      system("rm DATFILE");
    }
    dtfl_ptr = fopen("DATFILE", "w+");

    rewind(itfnm_ptr);

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
      nreps = *(pnreps+l);
      xnreps = (double) nreps;
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

      for(j=0;j<nreps*d;j++) fprintf(dtfl_ptr, "%g ", *(Y+j));
      fprintf(dtfl_ptr, "\n");

      /* compute per gene sample mean and residuals   */
      for(j=0;j<d;j++){
	sm = 0.0;
	for(i=0;i<nreps;i++) sm += *(Y + d*i + j);
	*(muhat + d*l + j) = sm/xnreps;
	for(i=0;i<nreps;i++) *(res+d*i+j) = *(Y +d*i +j) - (*(muhat + d*l + j));
      }

      /* compute per gene unbiased sample covariance matrix */
      smEV = 0.0;
      for(j=0;j<d;j++)
	for(k=0;k<d;k++){
	  sm = 0.0;
	  for(i=0;i<nreps;i++)
	    sm += *(res +d*i +j)*(*(res +d*i +k));
	  *(Sighat +d2*l +d*k +j) = sm/(xnreps -1.0);
          if(j==k) smEV += sm;
	}
      *(WSSQ + l) = smEV;
    }

    fclose(dtfl_ptr);

    /* Fit the model for the marginal distribution of Sighat under the Wishart/Inverse Wishart */
    y->MVM = Sighat;
    y->pN = pN;
    y->pd = pd;
    y->nreps = pnreps;

    for(l=0;l<npar;l++) *(ptheta0 + l) = 0.0;
    FitInvWish1(ptheta0, verb, y, objv, estimate, estimater, fail, fncnt, grcnt, 
                mask, usegr, G, H);

    /* Store the coefficient estimates for all simulation rounds:                              */
    for(l=0;l<npar;l++) *(coef + npar*isim + l) = *(estimate + l);

    /* Take the estimates from the model fit and use them to form 'nu_isim' and 'Lambda_isim'  */
    nu_isim = (exp(*estimate) + 1.0) * (2.0 * xd + 2.0);

    for(i=0;i<d;i++) *(LbdHlf + d*i + i) = exp(*(estimate + i + 1));

    k=0;
    for(i=0;i<d;i++)
      for(j=i+1;j<d;j++) {
        *(LbdHlf + d*j + i) = *(estimate + k + d + 1);
        k++;
      }

    for(i=0;i<d;i++)
      for(j=0;j<d;j++){
        sm = 0.0;
        for(k=0;k<d;k++) sm += *(LbdHlf + d*i + k) * (*(LbdHlf + d*j + k));
        *(Lambda_isim + d*j + i) = sm;
      }

    /* Fit the model for the marginal distribution of WSSQ                                        */
    /* under the scaled Chi-square/Inverse Gamma                                                  */
    yEV->S = WSSQ;
    yEV->pN = pN;
    yEV->pd = pd;
    yEV->nreps = pnreps;

    for(l=0;l<2;l++) *(ptheta0 + l) = 0.0;
    FitEqualVar1(ptheta0, verb, yEV, objv, estimate, fail, fncnt, grcnt, mask, usegr, G, H); 

    /* Store the coefficient estimates for all simulation rounds:                                 */
    for(l=0;l<2;l++) *(coefEV + 2*isim + l) = *(estimate + l);

    /* Take the estimates from the model fit and use them to form 's' and 'r'                     */
    s = exp(*estimate);
    r = exp(*(estimate+1));

    /* Now form the per gene statistics:                                                          */
    for(l=0;l<N;l++){
      nreps = *(pnreps + l);
      xnreps = (double) nreps;

      /*first calculate the per gene extreme mean value--store in vamx                            */
      vabsmax(muhat + d*l, d, k, vamx, sgn);

      /*shared variance Hotelling T^2 */
      for(i=0;i<d;i++)
        for(j=0;j<d;j++) 
          *(Sigma + d*j + i) = *(Lambda_isim + d*j + i) + (xnreps - 1.0)*(*(Sighat + d2*l + d*j + i));

      matinv(Sigma, SigInv, pd);
      for(j=0;j<d;j++){
        sm = 0.0;
        for(k=0;k<d;k++) sm += *(muhat + d*l + k) * (*(SigInv + d*j + k));
        *(xbuff + j) = sm;
      }
      sm=0.0;
      for(k=0;k<d;k++) sm += *(xbuff + k) * (*(muhat + d*l + k));
      xn2 = nu_isim + xnreps - 2.0*xd - 1.0;
      stat1 = xnreps * sm * xn2/xd;
      pval1 = pf(stat1, xd, xn2, 0, 0);

      /*ordinary Hotelling T^2        */
      for(i=0;i<d;i++)
        for(j=0;j<d;j++) 
          *(Sigma + d*j + i) = (xnreps - 1.0)*(*(Sighat + d2*l + d*j + i));

      matinv(Sigma, SigInv, pd);
      for(j=0;j<d;j++){
        sm = 0.0;
        for(k=0;k<d;k++) sm += *(muhat + d*l + k) * (*(SigInv + d*j + k));
        *(xbuff + j) = sm;
      }
      sm=0.0;
      for(k=0;k<d;k++) sm += *(xbuff + k) * (*(muhat + d*l + k));
      xn2 = xnreps - xd + 1.0;
      stat2 = xnreps * sm * xn2/xd;
      pval2 = pf(stat2, xd, xn2, 0, 0);

      /*shared variance univariate T2 */
      sm = 0.0;
      for(k=0;k<d;k++) sm+= (*(muhat + d*l + k))*(*(muhat + d*l + k));
      Top = (xnreps*xd-1.0)*(*(WSSQ+l)) + xnreps*sm;
      stat3 = Top/(2.0*r + (xnreps*xd-1.0)*(*(WSSQ+l)));
      stat3 = (2.0*s+xd*(xnreps-1.0))/(xnreps*xd)*stat3;
      pval3 = pf(stat3, xnreps*xd, 2.0*s + xd*(xnreps-1.0), 0, 0);

      /*ordinary univariate T2 */
      stat4 = Top/(*(WSSQ+l));
      stat4 = (xd*(xnreps-1.0))/(xnreps*xd)*stat4;
      pval4 = pf(stat4, xnreps*xd, xd*(xnreps-1.0), 0, 0);

      (genelist + l)->id = l;
      (genelist + l)->mean = vamx*sgn;
      (genelist + l)->ShHT2 = stat1;
      (genelist + l)->ShHT2pval = pval1;
      (genelist + l)->HT2 = stat2;
      (genelist + l)->HT2pval = pval2;
      (genelist + l)->ShUT2 = stat3;
      (genelist + l)->ShUT2pval = pval3;
      (genelist + l)->UT2 = stat4;
      (genelist + l)->UT2pval = pval4;
    }

    /*---ShHT2---BLOCK--------------------------------------------------------------------------------*/
    cmpr = &cmprShHT2;                                                                                
    qsort(genelist, N, sizeof(gene), cmpr);                                                           
                                      
    /* update ShHT2 portion of fdrtbl */                                                                
    for(j=0;j<nFDRlist;j++){                                                                          
      l=N;                                                                                            
      flagsig=0;                                                                                      
      while(flagsig==0&&l>0){                                                                         
        flagsig = 1*((genelist+l-1)->ShHT2pval <= *(FDRlist+j)*(*(rFDR+l-1)));                        
        l--;                                                                                          
      }                                                                                               
      Nsig = flagsig*(l++);                                                                           
      nTP = 0;                                                                                        
      nFP = 0;                                                                                        
      for(l=0;l<Nsig;l++) nTP += (fabs(*(es+(genelist+l)->id)) > 0.001);                              
      nFP = Nsig - nTP;                                                                               
      *(fdrtbl + (2*0+0)*nFDRlist + j) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);  
      *(fdrtbl + (2*0+1)*nFDRlist + j) += (Nsig > 0     ? ((double) nFP)/((double) Nsig) : 0.0);      
    }                                                                                                 
                                                                                                      
    /* update ShHT2 portion of roctbl  */                                                              
    nTP = 0;                                                                                          
    nFP = 0;                                                                                          
    for(l=0;l<N;l++){                                                                                 
      nTP += (fabs(*(es+(genelist+l)->id)) > 0.001);                                                  
      nFP = l+1 - nTP;                                                                                
      *(roctbl + N*0 + l) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);               
      *(roctbl + N*1 + l) += (nTP + nFP> 0 ? ((double) nFP)/((double) (nTP + nFP)) : 0.0);            
    }                                                                                                 
                                                                                                      
    /*---HT2---BLOCK----------------------------------------------------------------------------------*/
    cmpr = &cmprHT2;                                                                                  
    qsort(genelist, N, sizeof(gene), cmpr);                                                           
                                                                                                      
    /* update HT2 portion of fdrtbl */
    for(j=0;j<nFDRlist;j++){                                                                          
      l=N;                                                                                            
      flagsig=0;                                                                                      
      while(flagsig==0&&l>0){                                                                         
        flagsig = 1*((genelist+l-1)->HT2pval <= *(FDRlist+j)*(*(rFDR+l-1)));                          
        l--;                                                                                          
      }                                                                                               
      Nsig = flagsig*(l++);                                                                           
      nTP = 0;                                                                                        
      nFP = 0;                                                                                        
      for(l=0;l<Nsig;l++) nTP += 1*(fabs(*(es+l)) > 0.001);                                           
      nFP = Nsig - nTP;                                                                               
      *(fdrtbl + (2*1+0)*nFDRlist + j) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);  
      *(fdrtbl + (2*1+1)*nFDRlist + j) += (Nsig > 0     ? ((double) nFP)/((double) Nsig) : 0.0);      
    }                                                                                                 
                                                                                                      
    /* update HT2 portion of roctbl */
    nTP = 0;                                                                                          
    nFP = 0;                                                                                          
    for(l=0;l<N;l++){                                                                                 
      nTP += (fabs(*(es+(genelist+l)->id)) > 0.001);                                                  
      nFP = l+1 - nTP;                                                                                
      *(roctbl + N*2 + l) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);               
      *(roctbl + N*3 + l) += (nTP + nFP> 0 ? ((double) nFP)/((double) (nTP + nFP)) : 0.0);            
    }                                                                                                 
                                                                                                      
    /*---ShUT2---BLOCK--------------------------------------------------------------------------------*/
    cmpr = &cmprShUT2;                                                                                
    qsort(genelist, N, sizeof(gene), cmpr);                                                           
                                                                                                      
    /* update ShUT2 portion of fdrtbl */
    for(j=0;j<nFDRlist;j++){                                                                          
      l=N;                                                                                            
      flagsig=0;                                                                                      
      while(flagsig==0&&l>0){                                                                         
        flagsig = 1*((genelist+l-1)->ShUT2pval <= *(FDRlist+j)*(*(rFDR+l-1)));                        
        l--;                                                                                          
      }                                                                                               
      Nsig = flagsig*(l++);                                                                           
      nTP = 0;                                                                                        
      nFP = 0;                                                                                        
      for(l=0;l<Nsig;l++) nTP += 1*(fabs(*(es+l)) > 0.001);                                           
      nFP = Nsig - nTP;                                                                               
      *(fdrtbl + (2*2+0)*nFDRlist + j) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);  
      *(fdrtbl + (2*2+1)*nFDRlist + j) += (Nsig > 0     ? ((double) nFP)/((double) Nsig) : 0.0);      
    }                                                                                                 
                                                                                                      
    /* update ShUT2 portion of roctbl */
    nTP = 0;                                                                                          
    nFP = 0;                                                                                          
    for(l=0;l<N;l++){                                                                                 
      nTP += (fabs(*(es+(genelist+l)->id)) > 0.001);                                                  
      nFP = l+1 - nTP;                                                                                
      *(roctbl + N*4 + l) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);               
      *(roctbl + N*5 + l) += (nTP + nFP> 0 ? ((double) nFP)/((double) (nTP + nFP)) : 0.0);            
    }                                                                                                 
                                                                                                      
    /*---UT2---BLOCK----------------------------------------------------------------------------------*/
    cmpr = &cmprUT2;                                                                                  
    qsort(genelist, N, sizeof(gene), cmpr);                                                           

    /* update UT2 portion of fdrtbl */
    for(j=0;j<nFDRlist;j++){                                                                          
      l=N;                                                                                            
      flagsig=0;                                                                                      
      while(flagsig==0&&l>0){                                                                         
        flagsig = 1*((genelist+l-1)->UT2pval <= *(FDRlist+j)*(*(rFDR+l-1)));                          
        l--;                                                                                          
      }                                                                                               
      Nsig = flagsig*(l++);                                                                           
      nTP = 0;                                                                                        
      nFP = 0;                                                                                        
      for(l=0;l<Nsig;l++) nTP += 1*(fabs(*(es+l)) > 0.001);                                           
      nFP = Nsig - nTP;                                                                               
      *(fdrtbl + (2*3+0)*nFDRlist + j) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);  
      *(fdrtbl + (2*3+1)*nFDRlist + j) += (Nsig > 0     ? ((double) nFP)/((double) Nsig) : 0.0);      
    }                                                                                                 
                                                                                                      
    /* update UT2 portion of roctbl */
    nTP = 0;                                                                                          
    nFP = 0;                                                                                          
    for(l=0;l<N;l++){                                                                                 
      nTP += (fabs(*(es+(genelist+l)->id)) > 0.001);                                                  
      nFP = l+1 - nTP;                                                                                
      *(roctbl + N*6 + l) += (ntruepos > 0 ? ((double) nTP)/((double) ntruepos) : 0.0);               
      *(roctbl + N*7 + l) += (nTP + nFP> 0 ? ((double) nFP)/((double) (nTP + nFP)) : 0.0);            
    }                                                                                                 
    /*------------------------------------------------------------------------------------------------*/

    fprintf(itfnm_ptr, "isim=%d, nu=%g, L=",isim, nu_isim);
    for(i=0;i<d2;i++) fprintf(itfnm_ptr, "%g, ",*(Lambda_isim+i));
    fprintf(itfnm_ptr,"s=%g, r=%g\n",s,r);
  }
  fclose(itfnm_ptr);
  PutRNGstate();
}
