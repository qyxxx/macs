#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <Rcpp.h>
#include "R_ext/Print.h"

using namespace std;
using namespace Rcpp;

struct floatlink {
  float         present;
  struct floatlink *next;
};

typedef struct floatlink FELEMENT;
typedef FELEMENT *FLINK;

void   forward_step(double**, double**, double***, double***, int**, int*,
                    int**, int&, double**, double**, double**, double**,
                    double**, double**, double**, double**, double**, int,
                    int, int, int*, double*, NumericMatrix, double**, double*,
                    int, int**, int*, int*, int*, int, int, IntegerVector);
void   setup_1st_term(double**, double***, double*, NumericMatrix, double**,
                      int, int, int*, int, int);
NumericVector backward(double**, double**, double**, double***, double**,
                       double**, int, double*, NumericMatrix, int, int*, int*,
                       int, int,IntegerVector);
void   decomposite_residuals(double**, double*, double**, double**, double*, int,
                             int, double**, int*, int*, int, int);
void   output_decp_residuals(char*, double*, double**, double**, double**, double**,
                             int*, int*, int);
void   get_sqrt_mat(double**, double***, double***, double*, int*, int*, int, int);
void   deleteMatrix(double**, int);
void   deleteMatrix3(double***, int, int*);
void   copy_A_to_B3(double**, double**, int, int*);
void   distinct(double**, int**, int**, int*, int, int);
void   input(double*, int, int, int*, double**);
void   sorting_data(int*, int*, int, int, int, double**, int&);
void   greetings(void);
double **myMatrix(int, int);
double **myMatrix2(int, int*);
double ***myMatrix3(int, int*);
int    find_max(int*, int);
int    **myImatrix(int, int);
int    **myImatrix2(int, int*);


double (*nrfunc)(double*, double*, int, double**, int*, int*, int, int);

// [[Rcpp::export]]
List cmasal(DataFrame dataset, int MAXTERMN, int MAXTERM, int MAXITER, IntegerVector varhat, int enforce, IntegerVector number){
  //MAXTERMN: maximally allowed number of terms
  //MAXTERM:  maximally allowed number of interactions
  //varhat:   const1 + const2*time + const3*time^2
  //enforce:  Do you want to enforce certain terms into the model? (0 for no, 1 for yes)
  //number:   the corresponding covariate numbers
  int    n_covs, nT, i;
  int    *varNo;
  double *scan_datafile(DataFrame, int&, int&), *data;
  void copy_to_ntime(int *, IntegerVector, int);

  data = scan_datafile(dataset, nT, n_covs);

  if(enforce){
    varNo = (int *)malloc((number.size()+1)*sizeof (int));
    for(i=1; i<=number.size(); i++) varNo[i] = number[i-1];
    varNo[0] = number.size();
  }
  else{
    varNo = (int *)malloc((unsigned) 1*sizeof(int));
    varNo[0]=0;
  }

  int    n_sub;
  int    yes = (varNo[0] > 0);
  int    nterms;
  int    iter=1;
  int    MEMTERM;
  int    Mntime;
  int    j, k;
  int    *ntime;
  int    regression;
  int    *var_inclusion;
  double RSS;
  NumericMatrix varList(MAXTERMN+2, 3*MAXTERM+1);  // record bases in each term

  double **Pmat;
  //double **in_covariate = myMatrix(n_covs+2, nT+1);
  double **in_covariate = myMatrix(n_covs+3, nT+1);
  NumericMatrix indexpara(nT, n_covs+3);
  /*add one for the covariance column 03/03/04*/
  double **in_response;

  input(data, nT, n_covs, &n_sub, in_covariate);
  ntime = (int *) malloc((unsigned) (n_sub+1)*sizeof(int));
  sorting_data(ntime, &Mntime, nT, n_covs, n_sub, in_covariate, regression);

  var_inclusion= (int*) malloc ((unsigned) 4*sizeof(int));
  if(regression){
    var_inclusion[0] = 0;
    var_inclusion[1] = 0;
    var_inclusion[2] = 0;
  }
  //*blocked out to use princinpal components of the covariance matrix, 11/11/01
  else{
    var_inclusion[0] = varhat[0];
    var_inclusion[1] = varhat[1];
    var_inclusion[2] = varhat[2];
  }//*/

  in_response = myMatrix2(n_sub, ntime);
  int **maps = myImatrix2(n_sub, ntime);
  for(i=1, k=1; i<= n_sub; i++){
    for(j=1; j<=ntime[i]; j++){
      //in_response[i][j] = in_covariate[n_covs+1][k];
      in_response[i][j] = in_covariate[n_covs+2][k];
      /*add one for the covariance column*/
      maps[i][j] = k;
      k++;
    }
  }
  double ***V2 = myMatrix3(n_sub, ntime);
  double ***V = myMatrix3(n_sub, ntime);
  int  **rank = myImatrix(n_covs+1, nT);
  int  **ni   = myImatrix(n_covs+1, nT);
  //int  *mi    = (int *) malloc((unsigned) (n_covs+1)*sizeof(int));
  int  *mi    = (int *) malloc((unsigned) (n_covs+2)*sizeof(int));
  /*add one for the covariance column*/
  void   distinct(double**, int**, int**, int*, int, int);
  void   covmat(double*, double***, double*, int*, int*, int);
  distinct(in_covariate, rank, ni, mi, nT, n_covs);
  double **B = myMatrix2(n_sub, ntime);
  double **X = myMatrix2(n_sub, ntime);
  double **BX = myMatrix2(n_sub, ntime);
  double **BX1 = myMatrix2(n_sub, ntime);
  double **B1  = myMatrix2(n_sub, ntime);
  double **BX2 = myMatrix2(n_sub, ntime);
  double **B2  = myMatrix2(n_sub, ntime);
  double **newR = myMatrix2(n_sub, ntime);
  double **SP   = myMatrix(MAXTERMN+4, nT);
  double *Ui = (double *) malloc((unsigned) (n_sub+1)*sizeof(double));
  double **Vi = myMatrix(2, n_sub+1);
  double *sigma = (double *) malloc((unsigned) (5)*sizeof(double));
  void   estimate_sigmas(double*, int, double**, int*, int*, int, int);
  double *v4 = (double *) malloc((unsigned) (MAXTERMN+3)*sizeof(double));
  int    *rowANDcol = (int *) malloc((unsigned) (nT*2+1)*sizeof(int));
  IntegerVector outntime(n_sub);
  List outV;
  List copy_to_outV(double***, int, int *);
  IntegerVector finTerm;
  IntegerVector signal;
  signal = IntegerVector(100);

  for(i=0; i<=99; i++){
    signal[i]=100;
  }
  NumericVector beta;
  for(i=1, k=1; i<=n_sub; i++){
    for(j=1; j<=ntime[i]; j++, k++){
      rowANDcol[2*k-1]=i;
      rowANDcol[2*k]=j;
    }
  }
  double **oriR = myMatrix2(n_sub, ntime);
  int max_iter= (regression)? 1:MAXITER;
  for(iter=1; iter<=max_iter; iter++){
    if(iter==1){ //assume independence
      sigma[0]=1.;
      sigma[1]=0., sigma[2]=0., sigma[3]=0.;
    }
    get_sqrt_mat(in_covariate, V, V2, sigma, var_inclusion, ntime,
                 n_sub, n_covs);
    copy_A_to_B3(in_response, oriR, n_sub, ntime); //preserve original response
    //varList = (double **) malloc((unsigned) (MAXTERMN+2)*sizeof(double*));
      // if(iter == 2){
      // //varList = NumericMatrix(MAXTERMN+2, 3*MAXTERM+1);
      // for(i = 0;i < MAXTERM+2;i++){
      //   for(nnumber = 0;nnumber<3*MAXTERM+1;nnumber++){
      //     printf("(%d, %d) - %d\n",i,nnumber,varList(i,nnumber));
      //   }
      // }
      // }
    Pmat = (double **) malloc((unsigned) (MAXTERMN+2)*sizeof(double*));
    //response is modified from here
    setup_1st_term(in_response, V, &RSS, varList, Pmat, MAXTERMN, nT,
                   ntime, n_sub, n_covs);
    forward_step(in_covariate, in_response, V, V2, rank, mi, ni, nterms,
                 X, B, BX, B1, BX1, B2, BX2, newR, SP, yes, MAXTERMN,
                 MAXTERM, &MEMTERM, &RSS, varList, Pmat, v4, nT, maps,
                 varNo, rowANDcol, ntime, n_sub, n_covs, signal);
    copy_A_to_B3(oriR, in_response, n_sub, ntime);
    //move original response back
    //oriR will host final residuals after backward deletion
    finTerm = IntegerVector(nterms + 1);
    beta = backward(Pmat, B, in_covariate, V, oriR, X, nterms, &RSS,
                    varList, nT, varNo, ntime, n_sub, n_covs, finTerm);

    //add the following two lines at home 11/21/01
    for(k=1; k<=n_sub; k++) for(j=1; j<=ntime[k]; j++)
      //in_covariate[n_covs+1][maps[k][j]]=oriR[k][j];
      in_covariate[n_covs+2][maps[k][j]]=oriR[k][j];
    /*add one for the covarance column 03/03/04 */
    estimate_sigmas(sigma, Mntime, in_covariate, var_inclusion, ntime,
                    n_sub, n_covs);
    decomposite_residuals(oriR, Ui, Vi, B, sigma, Mntime, nT, in_covariate,
                          var_inclusion, ntime, n_sub, n_covs);
    copy_to_ntime(ntime, outntime, n_sub);
    //deleteMatrix(varList, nterms);
    deleteMatrix(Pmat, nterms);
  }

  covmat(in_covariate[n_covs+1], V, sigma, var_inclusion, ntime, n_sub);
  outV = copy_to_outV(V, n_sub, ntime);

  NumericVector sigma_out(4);
  for(i=0; i<4; i++) sigma_out[i] = sigma[i];

  for(i=0; i<n_covs+3; i++){
    for(j=0;j<nT; j++) indexpara(j, i) = in_covariate[i+1][j+1];
  }

  /*output_decp_residuals(outfile, Ui, Vi, B, in_response,
                          oriR, var_inclusion, ntime, n_sub);*/

  NumericVector eij(nT);
  for(i=1, k=0; i<=n_sub; i++){
    for(j=1; j<=ntime[i]; j++, k++){
      eij(k) = B[i][j];
    }
  }

  NumericMatrix coef(n_sub, 3);
  for(i=0; i<n_sub; i++){
    coef(i, 0) = Ui[i+1];
    coef(i, 1) = Vi[1][i+1];
    coef(i, 2) = Vi[2][i+1];
  }

  deleteMatrix3(V, n_sub, ntime);
  free((char*)v4);

  return List::create(
    Named("beta") = beta,
    Named("varList") = varList,
    Named("finTerm") = finTerm,
    Named("ntime") = outntime,
    Named("cov_matrix") = outV,
    Named("sigma") = sigma_out,
    Named("data_processed") = indexpara,
    // Named("rand_effects") = coef,
    // Named("meas_error") = eij,
    Named("nsub") = n_sub
  );
}

List copy_to_outV(double ***V, int n_sub, int *ntime)
{
  NumericMatrix convert_to_Matrix(double **, int, int);
  int i;
  List ls_V;
  for(i=1;i<=n_sub;i++)
  {
    ls_V.push_back(convert_to_Matrix(V[i],ntime[i],ntime[i]));
  }
  return ls_V;
}

NumericMatrix convert_to_Matrix(double **V,int col,int row)
{
  NumericMatrix out(col, row);
  int i,j;
  for(i=0;i<col;i++)
    for(j=0;j<row;j++)
      out(i,j) = V[i+1][j+1];
  return out;
}

void copy_to_ntime(int *ntime, IntegerVector outntime, int n_sub)
{
  int i;
  for(i=0;i<n_sub;i++)
    outntime[i] = ntime[i+1];
}

int find_max(int *vec, int num)
{
  int i,tempmax = 0;
  for(i=1;i<=num;i++)
  {
    if(vec[i]>tempmax)
      tempmax = vec[i];
  }
  return tempmax;
}

void forward_step(double **covariate, double **R, double ***V, double ***V2,
                  int **rank, int *mi, int **ni, int &nterms,
                  double **X, double **B, double **BX, double **B1,
                  double **BX1, double **B2, double **BX2, double **newR,
                  double **SP, int yes, int MAXTERMN, int MAXTERM, int *MEMTERM, double *RSS,
                  NumericMatrix varList, double **Pmat, double *v4, int nT, int **maps, int *varNo,
                  int *rowANDcol, int *ntime, int n_sub, int n_covs, IntegerVector signal)
{
  int    iter=3;
  int    *nonZeroi = (int *) malloc((unsigned) (3*nT+1)*sizeof(int));
  //new int[3*nT+1];
  int    i;
  void   forcedTerms(int, double**, int, NumericMatrix, double**, double**,
                     double***, double**, double**,double**, double*, double*, int,
                     int*, int, int);
  void   allocate_mem_bases(int, int, int*,  NumericMatrix, double**, int);
  int    addNewPairs(double**, int, NumericMatrix, double**, int**, int*, int**,
                     double**, double***, double***, double**, double**,
                     double**, double**, double**, double**, double**, double**,
                     double**, int*, int, double*, int, double*, int, int**, int*,
                     int*, int, int);
  nterms = *MEMTERM = 1;
  if(yes){
    for(i=1; i<=varNo[0]; i++){
      allocate_mem_bases(nterms, MAXTERM, MEMTERM, varList, Pmat, nT);
      forcedTerms(varNo[i], Pmat, nterms, varList, covariate,
                  R, V, B, B1, SP, RSS, v4, nT, ntime, n_sub, n_covs);
      nterms++;
    }
  }
  int t=0;
  while(iter > 0 && nterms < MAXTERMN){
    allocate_mem_bases(nterms, MAXTERM, MEMTERM, varList, Pmat, nT);
    iter=addNewPairs(Pmat, nterms, varList, covariate, rank, mi, ni, R,
                     V, V2, X, B, BX, B1, BX1, B2, BX2, newR, SP,
                     nonZeroi, MAXTERM, RSS, MAXTERMN, v4, nT, maps, rowANDcol,
                     ntime, n_sub, n_covs);
    if(iter) nterms++;
    nterms += (iter==3);
    signal[t]=iter;
    t++;
  }
  free((char*) nonZeroi);
}

void setup_1st_term(double **R, double ***V, double *RSS, NumericMatrix varList,
                    double **Pmat, int MAXTERMN, int nT, int *ntime, int n_sub, int n_covs)
{
  double *myVector(int);
  double *tmp = myVector(nT);
  double *tmp2 = myVector(nT);
  double firstResidual(double*, double**, int*, int, int);
  double a_dot_b(double*, double*, int);
  double  tmp1;
  void   copy_va_to_vb(double*, double*, int);
  void   A_mult_v(double**, double*, double*, int, int);
  void   set_v_to_a(double*, double, int);
  int     i, j, kk;

  *RSS=0.;
  for(j=1; j<=n_sub; j++) for(i=1; i<=ntime[j]; i++)
    *RSS += R[j][i]*R[j][i];

  Pmat[1] = (double *) malloc((unsigned) (nT+1)*sizeof(double));
  //new double[nT+1];
  tmp1 = 0.;
  for(j=1, kk=1; j<=n_sub; j++){
    set_v_to_a(tmp2, 1., ntime[j]);
    A_mult_v(V[j], tmp2, tmp, ntime[j], ntime[j]);
    tmp1 += a_dot_b(tmp, tmp, ntime[j]);
    for(i=1; i<= ntime[j]; i++, kk++)
      Pmat[1][kk] = tmp[i];
    A_mult_v(V[j], R[j], tmp2, ntime[j], ntime[j]);
    copy_va_to_vb(tmp2, R[j], ntime[j]);
  }
  tmp1 = sqrt(tmp1);
  //#ifdef DEBUG
  //cout << "normalizing constant " << tmp1 << endl;
  //#endif
  for(i=1; i<=nT; i++)  Pmat[1][i] /= tmp1;
  *RSS=0.;
  for(j=1; j<=n_sub; j++) for(i=1; i<=ntime[j]; i++)
    *RSS += R[j][i]*R[j][i];
  //#ifdef DEBUG
  //cout << "transformed RSS " << *RSS << endl;
  //#endif
  //  varList[1] = (double *) malloc((unsigned) 1*sizeof(double));
  //new double[1];
  varList(1,0)=0;
  *RSS=firstResidual(Pmat[1], R, ntime, n_sub, n_covs);

  free((char*) tmp);
  free((char*)tmp2);
}

void allocate_mem_bases(int nterms, int MAXTERM, int *MEMTERM, NumericMatrix varList,
                        double **Pmat, int nT)
{
  int itmp = *MEMTERM;
  if(itmp==nterms){
    Pmat[itmp+1]    = (double *) malloc((unsigned) (nT+1)*sizeof(double));
    Pmat[itmp+2]    = (double *) malloc((unsigned) (nT+1)*sizeof(double));
    itmp += 2;
  }
  else if(itmp==nterms+1){
    Pmat[itmp+1]    = (double *) malloc((unsigned) (nT+1)*sizeof(double));
    itmp++;
  }
  *MEMTERM = itmp;
}

double firstResidual(double *p1, double **R, int *ntime, int n_sub, int n_covs)
{
  double tmp=0., RSS1=0.;
  int    i, j, kk;

  for(j=1, kk=1; j<=n_sub; j++) for(i=1; i<=ntime[j]; i++, kk++)
    tmp += p1[kk]*R[j][i];
  for(j=1, kk=1; j<=n_sub; j++) for(i=1; i<=ntime[j]; i++, kk++){
    R[j][i] -= tmp*p1[kk];
    RSS1 += R[j][i]*R[j][i];
  }
  return RSS1;
}

void get_sqrt_mat(double **covariate, double ***V, double ***V2, double *sigma,
                  int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  void   covmat(double*, double***, double*, int*, int*, int);
  void   sqrtMat(double**, double**, int);
  void   invMat(double**, double**, int, int);
  int    i;

  //covmat(covariate[n_covs], V, sigma, var_inclusion, ntime, n_sub, n_covs);
  covmat(covariate[n_covs+1], V, sigma, var_inclusion, ntime, n_sub);
  /*add one for the covariance column 03/03/04 */
  for(i=1; i<=n_sub; i++){
    if(ntime[i] > 1){
      invMat(V[i], V2[i], ntime[i], 1);
      sqrtMat(V2[i], V[i], ntime[i]);
    }
    else{
      V2[i][1][1]=1./V[i][1][1];
      V[i][1][1] = sqrt(V2[i][1][1]);
    }
  }
}

void distinct(double **covariate, int **rank, int **ni, int *mi, int nT, int n_covs)
{
  int distinct_x(double*, int*, int*, int);

  for(int i=1; i<=n_covs; i++)
    mi[i]=distinct_x(covariate[i], rank[i], ni[i], nT);
}

int distinct_x(double *tmpx, int *rank, int *ni, int nT)
{
  int   i, j, k;
  void  indexx(double*, int*, int);
  double tmp;

  indexx(tmpx, rank, nT);
  k=1, j=rank[1], ni[k]=1;
  tmp = tmpx[j];            // current distinct value
  for(i=2; i<=nT; i++){
    j=rank[i];
    if(tmpx[j]==tmp)
      ++ni[k];
    else{
      tmp = tmpx[j];
      //       cout << tmp << " ";
      k = k+1;
      ni[k]=ni[k-1]+1;
    }
  }
  return k;
}

void indexx(double *arrin, int *indx, int n1)
{
  int    l,j,ir,indxt,i;
  double q;

  for (j=1;j<=n1;j++) indx[j]=j;
  l=(n1 >> 1) + 1;
  ir=n1;
  for (;;) {
    if (l > 1)
      q=arrin[(indxt=indx[--l])];
    else {
      q=arrin[(indxt=indx[ir])];
      indx[ir]=indx[1];
      if (ir-- == 1) {
        indx[1]=indxt;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
      if (q < arrin[indx[j]]) {
        indx[i]=indx[j];
        j += (i=j);
      }
      else j=ir+1;
    }
    indx[i]=indxt;
  }
}

void  covmat_1person(int i, int noccasion, double *time, double *sigma, double **V,
                     int kk, int *var_inclusion, int *ntime, int n_sub)
{
  int j, k;
  //sigma[0] = sigma^2, sigma[1] = sigma_u^2,
  //sigma[2] = sigma_v^2 for linear time
  //sigma[3] = sigma_v^2 for square time

  for(k=1; k<=noccasion; k++){
    V[k][k] = sigma[0] + var_inclusion[0]*sigma[1] +
      var_inclusion[1]*sigma[2]*time[kk+k]+
      var_inclusion[2]*sigma[3]*time[kk+k]*time[kk+k];
    for(j=k+1; j<=noccasion; j++){
      V[k][j] = var_inclusion[0]*sigma[1] +
        var_inclusion[1]*sigma[2]*sqrt(time[kk+j]*time[kk+k])+
        var_inclusion[2]*sigma[3]*time[kk+j]*time[kk+k];
      V[j][k] = V[k][j];
    }
  }
}

void  covmat(double *time, double ***V, double *sigma, int *var_inclusion,
             int *ntime, int n_sub)
{
  int i, j, k, kk;
  //sigma[0] = sigma^2, sigma[1] = sigma_u^2,
  //sigma[2] = sigma_v^2 for linear time
  //sigma[3] = sigma_v^2 for square time

  //    << "inside covmat " << sigma[0] << " " << sigma[1];
  //cout << " " << sigma[2] << endl;
  for(i=1, kk=1; i<=n_sub; i++){
    for(k=1; k<=ntime[i]; k++, kk++)
      V[i][k][k] = sigma[0] + var_inclusion[0]*sigma[1] +
        var_inclusion[1]*sigma[2]*time[kk]+
        var_inclusion[2]*sigma[3]*time[kk]*time[kk];
  }
  for(i=1, kk=1; i<=n_sub; i++){
    for(k=1; k<=ntime[i]; k++, kk++)
      for(j=k+1; j<=ntime[i]; j++){
        V[i][k][j] = var_inclusion[0]*sigma[1] +
          var_inclusion[1]*sigma[2]*sqrt(time[kk-k+j]*time[kk])+
          var_inclusion[2]*sigma[3]*time[kk-k+j]*time[kk];
        V[i][j][k] = V[i][k][j];
      }
  }
}

int addNewPairs(double **Pm, int k, NumericMatrix tmpList, double **covariate,
                int **rank, int *mi, int **ni, double **R, double ***V,
                double ***V2, double **X, double **B, double **BX,
                double **B1, double **BX1, double **B2, double **BX2,
                double **newR, double **SP, int *nonZeroi, int MAXTERM, double *RSS,
                int MAXTERMN, double *v4, int nT, int **maps, int *rowANDcol,
                int *ntime, int n_sub, int n_covs)
{
  double kn=0., nX, nPX, itmp, bestK=0.;
  double decreR=0., h1, minh=0.;
  int    maxiter=0;
  int    i, j, bestV=1, bestP=1;
  double bestKnot(double&, double**, double**, double**, double***,
                  int*, int, int*, double**, double**, int, double**,
                  double**, double**, double**, double, int*, int, int, int**, int*,
                  int*, int, int);
  void   newRes_and_Proj(int, int, NumericMatrix, double**, double**,
                         double***, double**, double**, double**, double*, double*, int,
                         int*, int, int);
  double linearTerm(double**, int, double**, double**, double***, double**,
                    double&, double&, double*, int, int*, int, int);
  double RSS1 = *RSS;
  void   update_resid_with_oneVec(double**, double*, int*, int, int);
  void   get_some_matrices(double***, double**, double**,
                           double**, double**, int, int, int*, int, int);
  void   extends(double*, double**, int*, int);
  void   A_dot_B3(double**, double**, double**, int, int*);
  void   copy_A_to_B3(double**, double**, int, int*);
  void   previous_basis(NumericVector, double**, double**, int*, int, int);
  void   update_varlist(NumericMatrix, int, int, int, double, int);
  int    checkPreviousTerms(int, NumericVector, int);

  for(i=1; i<=k; i++){   // loop for each of previous k bases
    previous_basis(tmpList(i, _), covariate, B, ntime, n_sub, n_covs); //reconstruct a previous basis
    itmp = tmpList(i,0)+1.;
    for(j=1; j<=n_covs; j++){  // loop for each of the covariates
      if(checkPreviousTerms(j, tmpList(i,_), MAXTERM)){
        extends(covariate[j], X, ntime, n_sub);
        copy_A_to_B3(R, newR, n_sub, ntime);
        A_dot_B3(B, X, BX, n_sub, ntime);
        h1=linearTerm(Pm, k, BX, R, V, B1, nX, nPX, v4, nT, ntime, n_sub, n_covs);
        if(h1 < RSS1*0.005 || nPX*nPX <= nX*0.015*itmp || nPX <= 0.) h1=0.;
        decreR = h1;
        if(h1 > 0.) update_resid_with_oneVec(newR, Pm[k+1], ntime, n_sub, n_covs);
        get_some_matrices(V, newR, Pm, B1, SP, k, (h1 > 0.), ntime, n_sub, n_covs);
        copy_A_to_B3(B1, newR, n_sub, ntime);
        kn=bestKnot(decreR, X, B, BX, V2, rank[j], mi[j], ni[j], SP,
                    newR, k+(h1>0.), BX1, BX2, B1, B2, itmp, nonZeroi, MAXTERMN, nT, maps,
                    rowANDcol, ntime, n_sub, n_covs);
        if(decreR > minh){
          minh = decreR;
          maxiter = (h1 > 0.) + 2*(decreR > h1+RSS1*0.01);
          bestV = j;
          bestK = kn;
          bestP = i;
        }
      }
    }
  }

  if(maxiter){
    update_varlist(tmpList, k, bestP, bestV, bestK, maxiter);
    newRes_and_Proj(k, maxiter, tmpList, Pm, R, V, covariate,
                    B, B1, RSS, v4, nT, ntime, n_sub, n_covs);\
  }
  return maxiter;
}

double bestKnot(double &DecreR, double **X, double **B, double **BX,
                double ***S2, int *rank, int mi, int *ni,
                double **SP, double **SR, int nterm, double **BX1,
                double **BX2, double **B1, double **B2,
                double terms, int *nonZeroi, int MAXTERMN, int nT, int **maps,
                int *rowANDcol, int *ntime, int n_sub, int n_covs)
{
  double tmp, kn1, newh=0.;
  double lowk, highk;
  double *ci = (double *) malloc((unsigned) 13*sizeof(double));
  // new double[13];
  double **myMatrix(int, int);
  double **ci35 = myMatrix(2, MAXTERMN+3);
  double candiateKnot(double*, double, double, double&, double);
  void   nonZeroIndex(int*, int*, int, int*, int, int**, int*);
  void   new_ci(double***, double**, double**, double**,
                double**, double**, double**, int*, int, double*, double**, int, int*);
  void   deleteMatrix(double**, int);
  void   set_A_to_a(double**, double, int, int);
  void   set_A_to_a3(double**, double, int, int*);
  void   set_v_to_a(double*, double, int);
  int    i, j, i1=1, i2, j1=1, j2;
  set_A_to_a3(BX1, 0., n_sub, ntime);
  set_A_to_a3(BX2, 0., n_sub, ntime);
  set_A_to_a3(B1, 0., n_sub, ntime);
  set_A_to_a3(B2, 0., n_sub, ntime);
  set_v_to_a(ci, 0., 11);         // initiate all c's to zero, or the
  // knots starts from the largest x
  set_A_to_a(ci35, 0., 2, MAXTERMN+2);
  j2 = rank[ni[mi]];
  i2 = rowANDcol[2*j2-1], j2 = rowANDcol[2*j2];
  double kn=X[i2][j2];    // giving the largest x
  highk = kn;
  for(i=mi-1; i>=1; i--){
    nonZeroIndex(rank, ni, i, nonZeroi, nT, maps, rowANDcol); // locate nonzeros in (X-kn)+
    j2 = ni[i+1]-ni[i];
    for(j=1; j<=j2; j++){
      i1=nonZeroi[3*j-2], j1=nonZeroi[3*j-1];
      B1[i1][j1]=B[i1][j1];
      BX1[i1][j1]=BX[i1][j1];
    }
    new_ci(S2, SR, SP, B1, B2, BX1, BX2, nonZeroi, j2, ci, ci35, nterm, ntime);
    j1 = rank[ni[i]];
    i1 = rowANDcol[2*j1-1], j1 = rowANDcol[2*j1];
    lowk = X[i1][j1];
    tmp=candiateKnot(ci, lowk, highk, kn1, terms);
    if(ni[mi]-ni[i] < 5 || ni[i] < 5) //avoiding small number of samples on
      tmp = 0.; //one side
    highk = lowk;
    if(tmp > newh){
      newh = tmp;
      kn = kn1;
    }
    for(j=1; j<=ni[i+1]-ni[i]; j++){
      i1=nonZeroi[3*j-2], j1=nonZeroi[3*j-1];
      B1[i1][j1]=0.;
      BX1[i1][j1]=0.;
      B2[i1][j1]=B[i1][j1];
      BX2[i1][j1]=BX[i1][j1];
    }
  }
  DecreR += newh;
  free((char*)ci);
  deleteMatrix(ci35, 2);
  return kn;
}

void nonZeroIndex(int *rank, int *ni, int i, int *nonZeroi, int nT, int **maps,
                  int *rowANDcol)
{
  int j;
  for(j=1; j<=ni[i+1]-ni[i]; j++){
    nonZeroi[3*j-2]=rowANDcol[2*rank[j+ni[i]]-1];
    nonZeroi[3*j-1]=rowANDcol[2*rank[j+ni[i]]];
    nonZeroi[3*j] = maps[nonZeroi[3*j-2]][nonZeroi[3*j-1]];
  }
}

double candiateKnot(double *ci1, double xi, double xi1, double &kn,
                    double terms)
{
  int i;
  double tmp, newh, tau;
  double ht(double*, double, double);

  newh=ht(ci1, xi, terms);
  kn = xi;
  tmp=ht(ci1, xi1, terms);
  if(tmp > newh){
    newh = tmp;
    kn = xi1;
  }
  tau = ci1[2]*ci1[4]-ci1[1]*ci1[5];
  if( (tau > 0. && ci1[2]*ci1[3]-ci1[1]*ci1[4] > xi*tau &&
      ci1[2]*ci1[3]-ci1[1]*ci1[4] < xi1*tau) ||
      (tau < 0. && ci1[2]*ci1[3]-ci1[1]*ci1[4] < xi*tau &&
      ci1[2]*ci1[3]-ci1[1]*ci1[4] > xi1*tau) ){
    i=(tau < 1.e-7 && tau > -1.e-7);
    tau = (ci1[2]*ci1[3]-ci1[1]*ci1[4])/tau;
    tmp = ht(ci1, tau, terms);
    if(tmp > newh){
      kn=tau;
      newh = tmp;
    }
  }
  (void) i;
  return newh;
}

double ht(double *ci1, double tau, double terms)
{
  double tmp = ci1[3]-2*ci1[4]*tau+ci1[5]*tau*tau;
  double tmp1 = ci1[8]-2*ci1[7]*tau+ci1[6]*tau*tau;

  if(tmp1 <= 0. || tmp <= 0. || tmp <= tmp1*0.015*terms) tmp=0.;
  else  tmp = (ci1[1]-ci1[2]*tau)*(ci1[1]-ci1[2]*tau)/tmp;
  return tmp;
}

void new_ci(double ***S2, double **SR, double **SP, double **B,
            double **B1, double **BX, double **BX1, int *index, int ni,
            double *ci1, double **ci35, int k, int *ntime)
{
  int    i;
  double recurrence_3to5(double*, double**, int*, int);
  double recurrence(double***, double**, double**, double**, double**,
                    int*, int, int*);
  double recurrence_1to2(double**, double**, int*, int);

  ci1[1] += recurrence_1to2(SR, BX, index, ni);
  ci1[2] += recurrence_1to2(SR, B, index, ni);
  ci1[6] += recurrence(S2, B, B1, B, B1, index, ni, ntime);
  ci1[7] += recurrence(S2, B, B1, BX, BX1, index, ni, ntime);
  ci1[8] += recurrence(S2, BX, BX1, BX, BX1, index, ni, ntime);
  ci1[3] = ci1[4] = ci1[5] = 0.;
  for(i=1; i<=k; i++){
    ci35[1][i] += recurrence_3to5(SP[i], BX, index, ni);
    ci35[2][i] += recurrence_3to5(SP[i], B, index, ni);
    ci1[3] += ci35[1][i]*ci35[1][i];
    ci1[4] += ci35[1][i]*ci35[2][i];
    ci1[5] += ci35[2][i]*ci35[2][i];
  }
  ci1[3] = ci1[8] - ci1[3];
  ci1[4] = ci1[7] - ci1[4];
  ci1[5] = ci1[6] - ci1[5];
  if(ci1[3] <= 0.) ci1[3] = ci1[4] = 0.;
  if(ci1[5] <= 0.) ci1[5] = ci1[4] = 0.;
}

double recurrence_1to2(double **SR, double **U, int *index, int ni)
{
  int   i, j, nsub;
  double tmp;

  tmp = 0.;
  for(nsub=1; nsub<=ni; nsub++){
    i=index[3*nsub-2], j=index[3*nsub-1];
    tmp += SR[i][j]*U[i][j];
  }
  return tmp;
}

double recurrence_3to5(double *p1, double **U, int *index, int ni)
{
  int   i, j, k, nsub;
  double tmp;

  tmp = 0.;
  for(nsub=1; nsub<=ni; nsub++){
    i=index[3*nsub-2], j=index[3*nsub-1];
    k = index[3*nsub];
    tmp += p1[k]*U[i][j];
  }
  return tmp;
}

double recurrence(double ***S, double **U, double **U1,
                  double **V, double **V1, int *index, int ni, int *ntime)
{
  int   i, j, t, nsub;
  double tmp;

  tmp = 0.;
  for(nsub=1; nsub<=ni; nsub++){
    i=index[3*nsub-2], j=index[3*nsub-1];
    for(t=1; t<=ntime[i]; t++){
      tmp += S[i][t][j]*(U1[i][t]*V[i][j]+V1[i][t]*U[i][j]+U[i][t]*V[i][j]);
    }
  }
  return tmp;
}

void  update_varlist(NumericMatrix tmpList, int nterms, int bestP, int bestV,
                     double kn, int iter)
{
  double ftmp =  tmpList(bestP,0);
  int    itmp=(int) ftmp;
  int    i, i1, jj;

  itmp *= 3;
  if(iter){
    if(iter==1){
      i = nterms+1;
      tmpList(i,0)=ftmp+1.;
      for(jj=1; jj<=itmp; jj++)
      tmpList(i,jj)=tmpList(bestP,jj);
      tmpList(i,itmp+2)=0.;
      tmpList(i,itmp+1)=(double) bestV;
    }
    else if(iter==2){
      i = nterms+1;
      tmpList(i,0)=ftmp+1.;
      for(jj=1; jj<=itmp; jj++)
        tmpList(i,jj)=tmpList(bestP,jj);
      tmpList(i,itmp+2)=1.;
      tmpList(i,itmp+1)=(double) bestV;
      tmpList(i,itmp+3) = kn;
    }
    else if(iter==3){
      i = nterms+1, i1 = i+1;
      tmpList(i,0)=ftmp+1.;
      tmpList(i1,0)=ftmp+1.;
      for(jj=1; jj<=itmp; jj++)
        tmpList(i,jj)=tmpList(i1,jj)=tmpList(bestP,jj);
      tmpList(i,itmp+2)=0.;
      tmpList(i,itmp+1)=(double) bestV;
      tmpList(i1,itmp+2)=1.;
      tmpList(i1,itmp+1)=(double) bestV;
      tmpList(i1,itmp+3) = kn;
    }
    else{
      // Rcpp::stop("no such iter");
      // exit( -1 );
    }
  }

}

double linearTerm(double **Pmat, int nterm, double **BX, double **R,
                  double ***V, double **X2, double &nX, double &nPX, double *v4, int nT,
                  int *ntime, int n_sub, int n_covs)
{
  double tmp;
  double proMatToMat1(double**, int, double***, double**, double**, double&, double*,
                      int, int*, int);
  double matVecNorm(double**, double*, int, int*);

  nPX=proMatToMat1(Pmat, nterm, V, BX, X2, nX, v4, nT, ntime, n_sub);
  if(nPX) tmp = matVecNorm(R, Pmat[nterm+1], n_sub, ntime);
  else tmp=0.;
  return tmp;
}

void update_resid_with_oneVec(double **newR, double *newPm, int *ntime, int n_sub, int n_covs)
{
  int    i, j, jn;
  double tmp=0.;

  for(i=1, jn=1; i<=n_sub; i++)
    for(j=1; j<=ntime[i]; j++, jn++)
      tmp += newR[i][j]*newPm[jn];
  for(i=1, jn=1; i<=n_sub; i++)
    for(j=1; j<=ntime[i]; j++, jn++)
      newR[i][j] -= tmp*newPm[jn];
}

void get_some_matrices(double ***S, double **R, double **Pmat,
                       double **SR, double **SP, int nterm, int yes, int *ntime, int n_sub, int n_covs)
{
  void   set_A_to_a3(double**, double, int, int*);
  int i, j, k, t, jn, nt=nterm+yes;
  set_A_to_a3(SR, 0., n_sub, ntime);

  for(i=1, jn=1; i<=n_sub; i++){
    for(j=1; j<=ntime[i]; j++, jn++){
      SP[nt][jn]=SP[nterm][jn]=0.;
      if(nterm > 1) SP[nterm-1][jn]=0.;
      for(k=1; k<=ntime[i]; k++){
        SR[i][j] += S[i][j][k]*R[i][k];
        for(t=nterm; t>=1 && t>=nterm-1; t--)
          SP[t][jn] += S[i][j][k]*Pmat[t][jn-j+k];
        if(yes) SP[nt][jn] += S[i][j][k]*Pmat[nt][jn-j+k];
      }
    }
  }
}

int checkPreviousTerms(int j, NumericVector tmpList, int MAXTERM)
{
  int k;

  if(tmpList[0] >= MAXTERM) return 0;
  for(k=1; k<=(int) tmpList[0]; k++)
    if(j==int(tmpList[3*k-2])) return 0;
    return 1;
}

void previous_basis(NumericVector varl, double **covariate, double **B, int *ntime,
                    int n_sub, int n_covs)
{
  void   set_A_to_a3(double**, double, int, int*);
  int varNo1, nterms, nt;
  int j, k, jn;

  set_A_to_a3(B, 1., n_sub, ntime);
  for(nterms=1; nterms<= (int) varl[0]; nterms++){
    nt = 3*nterms;
    varNo1 = (int) varl[nt - 2];
    for(j=1, jn=1; j<=n_sub; j++){
      for(k=1; k<=ntime[j]; k++, jn++){
        if(varl[nt - 1]){
          if(covariate[varNo1][jn]-varl[nt] > 0.)
            B[j][k] *= covariate[varNo1][jn]-varl[nt];
          else B[j][k] = 0.0;
        }
        else
          B[j][k] *= covariate[varNo1][jn];
      }
    }
  }
}

void  extends(double *allX, double **X, int *ntime, int n_sub)
{
  int jn, j, k;

  for(j=1, jn=1; j<=n_sub; j++)
    for(k=1; k<=ntime[j]; k++, jn++)
      X[j][k] = allX[jn];
}

double proMatToMat(double **P1, int nterm, double ***V1, double **X1,
                   double **X2, double *v4, int nT, int *ntime, int n_sub)
{
  int    i;
  double proVecToMat(double**, int, double**, double*, double*, int, int*, int);
  void   A_mult_v(double**, double*, double*, int, int);

  for(i=1; i<=n_sub; i++)
    A_mult_v(V1[i], X1[i], X2[i], ntime[i], ntime[i]);
  return proVecToMat(P1, nterm, X2, P1[nterm+1], v4, nT, ntime, n_sub);
}

double proVecToMat(double **P1, int nterm, double **X2, double *v2, double *v4, int nT,
                   int *ntime, int n_sub)
{
  int    i, j, k, jn;
  double pnorm=0.;
  void mat_mult_mat3(double**, int, double**, double*, int, int*);

  double tmp=0.;
  mat_mult_mat3(P1, nterm, X2, v4, n_sub, ntime);
  for (i=1, jn=1; i<= n_sub; i++){
    for(j=1; j<=ntime[i]; j++, jn++){
      v2[jn] = X2[i][j];
      tmp += X2[i][j]*X2[i][j];
      for (k=1; k<= nterm; k++)
        v2[jn] -= P1[k][jn]*v4[k];
      pnorm += v2[jn]*v2[jn];
    }
  }
  if(tmp <= 0.) tmp = 0.;
  if(tmp > 0. && pnorm > tmp*1.e-10){
    pnorm = sqrt(pnorm);
    for (i=1; i<= nT; i++)
      v2[i] /= pnorm;
  }
  else   pnorm = 0.;
  return pnorm;
}


double proMatToMat1(double **P1, int nterm, double ***V1, double **X1,
                    double **X2, double &nX, double *v4, int nT, int *ntime, int n_sub)
{
  int    i;
  double proVecToMat1(double**, int, double**, double*, double&, double*, int, int*,
                      int);
  void   A_mult_v(double**, double*, double*, int, int);

  for(i=1; i<=n_sub; i++)
    A_mult_v(V1[i], X1[i], X2[i], ntime[i], ntime[i]);
  return proVecToMat1(P1, nterm, X2, P1[nterm+1], nX, v4, nT, ntime, n_sub);
}

// return a residual double * of v1 onto double ** P1

double proVecToMat1(double **P1, int nterm, double **X2, double *v2,
                    double &tmp, double *v4, int nT, int *ntime, int n_sub)
{
  int    i, j, k, jn;
  double pnorm=0.;
  void   mat_mult_mat3(double**, int, double**, double*, int, int*);

  tmp=0.;
  mat_mult_mat3(P1, nterm, X2, v4, n_sub, ntime);
  for (i=1, jn=1; i<= n_sub; i++){
    for(j=1; j<=ntime[i]; j++, jn++){
      v2[jn] = X2[i][j];
      tmp += X2[i][j]*X2[i][j];
      for (k=1; k<= nterm; k++)
        v2[jn] -= P1[k][jn]*v4[k];
      pnorm += v2[jn]*v2[jn];
    }
  }
  if(tmp <= 0.) tmp = 0.;
  if(tmp > 0. && pnorm > tmp*1.e-10){
    pnorm = sqrt(pnorm);
    for (i=1; i<= nT; i++)
      v2[i] /= pnorm;
  }
  else   pnorm = 0.;
  return pnorm;
}

void newRes_and_Proj(int nterms, int iter, NumericMatrix tmpList, double **Pm,
                     double **R, double ***V, double **covariate,
                     double **B, double **X, double *RSS, double *v4, int nT,
                     int *ntime, int n_sub, int n_covs)
{
  int    i, j;
  double tmp;
  double proMatToMat(double**, int, double***, double**, double**, double*, int,
                     int*, int);
  void   previous_basis(NumericVector , double**, double**, int*, int, int);

  if(iter > 0){
    previous_basis(tmpList(nterms+1,_), covariate, B, ntime, n_sub, n_covs);
    tmp=proMatToMat(Pm, nterms, V, B, X, v4, nT, ntime, n_sub);
    if(tmp <= 0.){
      //Rcpp::stop("cannot add linear term");
      // cout << "cannot add linear term " << nterms+1 << endl;
      // cout << "nP " << tmp << endl;
      // cerr << "must be positive\n";
      //exit( -1 );
    }
    update_resid_with_oneVec(R, Pm[nterms+1], ntime, n_sub, n_covs);
    if(iter == 3){
      previous_basis(tmpList(nterms+2,_), covariate, B, ntime, n_sub, n_covs);
      tmp=proMatToMat(Pm, nterms+1, V, B, X, v4, nT, ntime, n_sub);
      if(tmp <= 0.){
        //Rcpp::stop("cannot add linear term");
        // cout << "cannot add term " << nterms+2 << endl;
        // cout << "nP " << tmp << endl;
        // cerr << "must be positive\n";
        //exit( -1 );
      }
      update_resid_with_oneVec(R, Pm[nterms+2], ntime, n_sub, n_covs);
    }
  }
  *RSS=0.;
  for(i=1; i<=n_sub; i++){
    for(j=1; j<=ntime[i]; j++){
      (*RSS) += R[i][j]*R[i][j];
    }
  }
}

void forcedTerms(int varNo1, double **Pm, int nterm, NumericMatrix tmpList,
                 double **covariate, double **R, double ***V,
                 double **B, double **B1, double **SP, double *RSS, double *v4, int nT,
                 int *ntime, int n_sub, int n_covs)
{
  int i, j, k, t, jn;

  for(i=1, jn=1; i<=n_sub; i++){
    for(j=1; j<=ntime[i]; j++, jn++){
      SP[nterm+1][jn]=SP[nterm][jn]=0.;
      if(nterm > 1) SP[nterm-1][jn]=0.;
      for(k=1; k<=ntime[i]; k++){
        for(t=nterm; t>=1 && t>=nterm-1; t--)
          SP[t][jn] += V[i][j][k]*Pm[t][jn-j+k];
        SP[nterm+1][jn] += V[i][j][k]*Pm[nterm+1][jn-j+k];
      }
    }
  }
  update_varlist(tmpList, nterm, 1, varNo1, 0., 1);
  newRes_and_Proj(nterm, 1, tmpList, Pm, R, V, covariate,
                  B, B1, RSS, v4, nT, ntime, n_sub, n_covs);
}

void sqrtMat(double **A, double **sA, int n)
{
  void svdcmp(double**, double*, double**, int, int);
  void A_mult_B(double**, double**, double**, int, int, int);
  void A_mult_Bt(double**, double**, double**, int, int, int);
  void   deleteMatrix(double**, int);
  void   copy_A_to_B(double**, double**, int, int);
  void   set_A_to_a(double**, double, int, int);
  double *myVector(int);
  double **myMatrix(int, int);
  double *w = myVector(n);
  double **B = myMatrix(n, n);
  double **B1 = myMatrix(n, n);
  double **D = myMatrix(n,n);
  int i;

  copy_A_to_B(A, B, n, n);
  svdcmp(B, w, sA, n, n);
  set_A_to_a(D, 0., n, n);
  for(i=1; i<=n; i++){
    D[i][i] = sqrt(w[i]);
  }
  A_mult_B(B, D, B1, n, n, n);
  A_mult_Bt(B1, sA, B, n, n, n);
  copy_A_to_B(B, sA, n, n);
  deleteMatrix(D, n);
  deleteMatrix(B, n);
  deleteMatrix(B1, n);
  free((char*) w);
}

void invMat(double **A, double **sA, int n, int enforcing)
{
  void A_mult_B(double**, double**, double**, int, int, int);
  void A_mult_Bt(double**, double**, double**, int, int, int);
  void   deleteMatrix(double**, int);
  if(n==1){
    if(A[1][1]==0.){
      //Rcpp::stop("Singular Covariance Matric: invMat");
      // cout << "Singular Covariance Matric: invMat" << endl;
      //exit(-1);
    }
    else sA[1][1]=1./A[1][1];
    return;
  }
  void svdcmp(double**, double*, double**, int, int);
  void   copy_A_to_B(double**, double**, int, int);
  void   set_A_to_a(double**, double, int, int);
  double *myVector(int);
  double **myMatrix(int, int);
  double *w = myVector(n);
  double **B = myMatrix(n, n);
  double **B1 = myMatrix(n, n);
  double **D = myMatrix(n,n);
  //double delta;
  int i;

  copy_A_to_B(A, B, n, n);
  svdcmp(B, w, sA, n, n);
  set_A_to_a(D, 0., n, n);
  for(i=1; i<=n; i++){
    D[i][i] = 1./w[i];
  }
  A_mult_Bt(B, sA, B1, n, n, n); //check the inverse matrix
  /*for(i=1; i<=n; i++) for(j=1; j<=n; j++){
    delta = (i==j)? 1-B1[i][j] : B1[i][j];
    if(delta > 1.e-7 || delta < -1.e-7){
      cout << "Inverse nearly singular matrix: invMat" << endl;
      cout << "delta: "<< delta << endl;
      if(!enforcing) Rcpp::stop("Inverse nearly singular matrix: invMat");
    }
  }*/
  A_mult_B(B, D, B1, n, n, n);
  A_mult_Bt(B1, sA, B, n, n, n);
  copy_A_to_B(B, sA, n, n);
  deleteMatrix(D, n);
  deleteMatrix(B, n);
  deleteMatrix(B1, n);
  free((char*)w);
}

NumericVector backward(double **iPP, double **B, double **covariate, double ***V1,
                       double **R, double **X2, int nterm, double *RSS, NumericMatrix varList,
                       int nT, int *varNo, int *ntime, int n_sub, int n_covs,IntegerVector finTerm)
{
  void  Final_formula(int, char *, int*, NumericVector, double**);
  double Gcv(double, int, int);
  double *myVector(int);
  void   deleteMatrix3(double***, int, int*);

  int    deli, i, j, j1, k, jn, bnterm;
  int    RSS_del_one(double**, double*, double*, int, double&, int*);
  //new int[nterm+1];
  int    *delTerm=(int *) malloc((unsigned) (nterm+1)*sizeof(int));
  //new int[nterm+1];
  double gcv, mgcv, tmp;
  double **myMatrix(int, int);
  double ***myMatrix3(int, int*);
  double *vtmp = myVector(nT);
  double **ixx = myMatrix(nterm, nterm);
  double **bxx = myMatrix(nterm, nterm);
  double *yx  =  myVector(nterm);
  double *byx =  myVector(nterm);
  double *x1  =  myVector(nterm);
  void   deleteMatrix(double**, int);
  void   copy_A_to_B(double**, double**, int, int);
  void   copy_va_to_vb(double*, double*, int);
  void   A_mult_v(double**, double*, double*, int, int);
  void   A_mult_v1(double**, double*, NumericVector, int, int);
  void   invMat(double**, double**, int, int);

  for(i=1; i<=nterm; i++) finTerm[i]=delTerm[i]=i;
  tmp=0.;
  for(i=1; i<=n_sub; i++){
    A_mult_v(V1[i], R[i], vtmp, ntime[i], ntime[i]);
    copy_va_to_vb(vtmp, R[i], ntime[i]);
  }
  for(k=1; k<=nterm; k++){
    yx[k]=0.;
    previous_basis(varList(k,_), covariate, B, ntime, n_sub, n_covs);
    for(i=1, jn=1; i<=n_sub; i++){
      A_mult_v(V1[i], B[i], vtmp, ntime[i], ntime[i]);
      for(j=1; j<=ntime[i]; j++, jn++){
        iPP[k][jn]=vtmp[j];
        yx[k] += R[i][j]*vtmp[j];
      }
    }
  }

  for(i=1; i<=nterm; i++) for(j=1; j<=nterm; j++){
    bxx[i][j]=0.;
    for(k=1; k<=nT; k++)
      bxx[i][j] += iPP[i][k]*iPP[j][k];
  }
  invMat(bxx, ixx, nterm, 1);
  *RSS=0.;
  for(i=1; i<=n_sub; i++)
    for(j=1; j<=ntime[i]; j++)
      *RSS += R[i][j]*R[i][j];
  for(i=1; i<=nterm; i++)
    for(j=1; j<=nterm; j++)
      *RSS -= yx[i]*ixx[i][j]*yx[j];
  mgcv = Gcv(*RSS, nterm, nT);
  bnterm = nterm;
  copy_A_to_B(ixx, bxx, nterm, nterm);
  copy_va_to_vb(yx, byx, nterm);
  double bestRSS= *RSS;
  for(i=nterm; i>2+varNo[0]; i--){
    deli = RSS_del_one(ixx, yx, x1, i, tmp, varNo);
    *RSS += tmp;
    gcv = Gcv(*RSS, i-1, nT);
    for(j=deli; j<i; j++)
      delTerm[j]=delTerm[j+1];
    if(gcv <= mgcv){
      mgcv = gcv;
      for(j=1; j<i; j++)
        finTerm[j]=delTerm[j];
      bestRSS = *RSS;
      copy_A_to_B(ixx, bxx, i-1, i-1);
      copy_va_to_vb(yx, byx, i-1);
      bnterm = i-1;
    }
  }

  NumericVector beta(bnterm+1);
  A_mult_v1(bxx, byx, beta, bnterm, bnterm);

  tmp = 0.;
  for(j=1, jn=1; j<=n_sub; j++) for(j1=1; j1<=ntime[j]; j1++, jn++){
    for(i=1; i<=bnterm; i++)
      R[j][j1] -= iPP[finTerm[i]][jn]*beta[i];
    tmp += R[j][j1]*R[j][j1];
  }
  if((bestRSS-tmp)/bestRSS > 1.e-5 || (bestRSS-tmp)/bestRSS < -1.e-5){
    //Rcpp::stop("RSS is not consistent");
    // cout << "RSS is not consistent " << bestRSS << " " << tmp << endl;
    //exit(-1);
  }
  double ***V2=myMatrix3(n_sub, ntime);
  for(i=1; i<=n_sub; i++){
    invMat(V1[i], V2[i],  ntime[i], 1);
    A_mult_v(V2[i], R[i], vtmp, ntime[i], ntime[i]);
    copy_va_to_vb(vtmp, R[i], ntime[i]);
  }
  //output residuals
  //void output_residuals(char*, double*, double**, int*, int);
  //output_residuals(outfile, covariate[n_covs], R, ntime, n_sub, n_covs);
  //output_residuals(outfile, covariate[n_covs+1], R, ntime, n_sub);
  /*add one for the covariance column 03/03/04 */

  //Final_formula(bnterm, outfile, finTerm, beta, varList);
  deleteMatrix3(V2, n_sub, ntime);
  deleteMatrix(bxx, nterm);
  deleteMatrix(ixx, nterm);
  free((char*)vtmp);
  free((char*)delTerm);
  free((char*)yx);
  free((char*)byx);
  free((char*) x1);
  return beta;
}

int RSS_del_one(double **ixx, double *yx, double *x1,
                int nterm, double &miny, int *varNo)
{
  int    i, k, deli;
  double tmp;
  void   copy_va_to_vb(double*, double*, int);

  tmp=0.;
  for(k=1; k<=nterm; k++) tmp += ixx[ 2+varNo[0] ][k]*yx[k];
  miny = tmp*tmp/ixx[ 2+varNo[0] ][ 2+varNo[0] ];
  deli = 2+varNo[0];
  //forced terms not be deleted
  for(i=3+varNo[0]; i<=nterm; i++){
    tmp = 0.;
    for(k=1; k<=nterm; k++) tmp += ixx[i][k]*yx[k];
    if(tmp*tmp/ixx[i][i] < miny){
      miny = tmp*tmp/ixx[i][i];
      deli = i;
    }
  }
  copy_va_to_vb(ixx[deli], x1, nterm);
  for(k=deli; k<nterm; k++){
    for(i=1; i<=nterm; i++)
      ixx[k][i] = ixx[k+1][i];
    yx[k] = yx[k+1];
  }
  copy_va_to_vb(x1, ixx[nterm], nterm);
  for(i=1; i<=nterm; i++) x1[i] = ixx[i][deli];
  for(k=deli; k<nterm; k++)
    for(i=1; i<=nterm; i++)
      ixx[i][k] = ixx[i][k+1];
  for(k=1; k<nterm; k++)
    for(i=1; i<nterm; i++)
      ixx[k][i] -= x1[k]*x1[i]/x1[nterm];
  return deli;
}


double  Gcv(double rss, int nterm, int nT)
{
  double tmp = 1.-(1.+(nterm-1.)*5.)/nT;
  //  double tmp = 1.-(1.+(nterm-1.)*4)/nT;
  return rss/(tmp*tmp);
}

void output_residuals(char *outfile1, double *lastvar, double **R,
                      int *ntime, int n_sub)
{
  char outfile[100];
  char datafile[100];
  int ind=0;
  int i=0, j, j1, jn;
  void trans1(char*);

  strcpy(datafile, outfile1);
  trans1(datafile);
  do{
    if(datafile[ind]=='.')
      i=0;
    else if((datafile[ind]!='\\')&&(datafile[ind]!='/'))
      outfile[i++] = datafile[ind];
  } while((datafile[ind]!='\0')&&(datafile[ind]!='\\')
            &&(datafile[ind++]!='/'));
  outfile[i--]='\0';
  trans1(outfile);
  strcat(outfile,".res");
  ofstream residFile(outfile, ios::out);
  if (!residFile){
    //cerr << "cannot open file: " << outfile << endl;
    //exit( -1 );
  }
  for(j=1; j<=n_sub; j++) residFile << ntime[j] << " ";
  residFile << endl;
  for(j=1, jn=1; j<=n_sub; j++) for(j1=1; j1<=ntime[j]; j1++, jn++){
    residFile << lastvar[jn] << " " << R[j][j1] << endl;
  }
}

void trans1(char *s)
{
  int c,i,j;
  for(i=0, j=strlen(s) -1; i<j;i++,j--){
    c=s[i];
    s[i] = s[j];
    s[j] = c;
  }
}

void estimate_sigmas(double *sigma, int Mntime, double **in_covariate,
                     int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  double ftmp;
  double minization(double**, double*,
                    int, double (*)(double*, double*, int,
                                 double**, int*, int*, int, int),
                                 double*, double*, int, double*, int, double**,
                                 int*, int*, int, int);
  double neg_likelihood_for_sigma(double*, double*, int, double**, int*,
                                  int*, int, int);
  double *q1;
  double *psum, *ptry;
  double **myMatrix(int, int);
  void   deleteMatrix(double**, int);
  int    i, j, method1 = 0;

  q1 = (double *) malloc((unsigned) 6*sizeof(double));
  //new double[6];
  psum = (double *) malloc((unsigned) 5*sizeof(double));
  //new double[5];
  ptry = (double *) malloc((unsigned) 5*sizeof(double));
  //new double[5];


  q1[1]=sqrt(sigma[0]);
  for(i=2, j=2; i<=4; i++)
    if(var_inclusion[i-2]) q1[j++]=sqrt(sigma[i-1]);
    double **h1 = myMatrix(6,5);
    ftmp = minization(h1, q1, j-1, neg_likelihood_for_sigma, psum, ptry, method1,
                      sigma, Mntime, in_covariate, var_inclusion, ntime, n_sub, n_covs);
    (void) ftmp;
    sigma[0] = q1[1]*q1[1];

    for(i=2, j=2; i<=4; i++){
      if(var_inclusion[i-2]){
        sigma[i-1] = q1[j]*q1[j];
        j++;
      }
    }

    free((char*)q1);
    free((char*) psum);
    free((char*) ptry);
    deleteMatrix(h1, 6);
}

double neg_likelihood_for_sigma(double *sigma1, double *sigma, int Mntime,
                                double **in_covariate, int *var_inclusion,
                                int *ntime, int n_sub, int n_covs){
  void   covmat_1person(int, int, double*, double*, double**, int, int*,
                        int*, int);
  void   svdcmp(double**, double*, double**, int, int);
  void   A_mult_B(double**, double**, double**, int, int, int);
  void   A_mult_Bt(double**, double**, double**, int, int, int);
  void   deleteMatrix(double**, int);
  void   A_mult_v(double**, double*, double*, int, int);
  void   set_A_to_a(double**, double, int, int);
  double **myMatrix(int, int);
  //change it at home 11/21/01
  double *w =  (double *) malloc((unsigned) (Mntime+1)*sizeof(double));
  double **B1 = myMatrix(Mntime, Mntime);
  double **D = myMatrix(Mntime, Mntime);
  double **covM1 = myMatrix(Mntime, Mntime);
  double **covM2 = myMatrix(Mntime, Mntime);
  double lkh;
  int    i, j, kk;

  // cout << "sigma: ";
  sigma[0]=sigma1[1]*sigma1[1];
  for(i=2, j=2; i<=4; i++){
    if(var_inclusion[i-2]){
      sigma[i-1] =sigma1[j]*sigma1[j];
      j++;
    }
    else sigma[i-1]=0.;
    // cout << sigma[i-1] << " ";
  }
  lkh = 0.;
  for(j=1, kk=0; j<=n_sub; j++){
    //covmat_1person(j, ntime[j], in_covariate[n_covs], sigma, covM1, kk,
    covmat_1person(j, ntime[j], in_covariate[n_covs+1], sigma, covM1, kk,
                   var_inclusion, ntime, n_sub);
    /*add one for the covariance column */
    svdcmp(covM1, w, covM2, ntime[j], ntime[j]);
    set_A_to_a(D, 0., ntime[j], ntime[j]);
    for(i=1; i<=ntime[j]; i++){
      D[i][i] = 1./w[i];
      lkh -= log(w[i]);
    }
    A_mult_Bt(covM1, covM2, B1, ntime[j], ntime[j], ntime[j]);
    //check the inverse matrix
    /*for(i=1; i<=ntime[j]; i++) for(k=1; k<=ntime[j]; k++){
      delta = (i==k)? 1-B1[i][k] : B1[i][k];
      if(delta > 1.e-7 || delta < -1.e-7){
        Rcpp::stop("Inverse nearly singular matrix: neg_likelihood_for_sigm");
        //cout << "Inverse nearly singular matrix: neg_likelihood_for_sigma" << endl;
        //cout << "delta: " << delta << endl;
        //exit(-1);
      }
    }*/
    A_mult_B(covM1, D, B1, ntime[j], ntime[j], ntime[j]);
    A_mult_Bt(B1, covM2, covM1, ntime[j], ntime[j], ntime[j]);
    for(i=1; i<=ntime[j]; i++) //D[1][i] = in_covariate[n_covs+1][kk+i];
      D[1][i] = in_covariate[n_covs+2][kk+i];
    /*add one for the covariance column */
    A_mult_v(covM1, D[1], w, ntime[j], ntime[j]);
    for(i=1; i<=ntime[j]; i++){
      lkh -= w[i]*D[1][i];
    }
    kk += ntime[j];
  }
  deleteMatrix(D, Mntime);
  deleteMatrix(B1, Mntime);
  deleteMatrix(covM1, Mntime);
  deleteMatrix(covM2, Mntime);
  free((char*) w);
  //cout << sigma[0] << " " << sigma[1] << " " << lkh << endl;
  return -lkh;
}

void decomposite_residuals(double **resid, double *Ui, double **Vi,
                           double **Eij, double *sigma, int Mntime, int nT, double **in_covariate,
                           int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  //resid[i] is the residual vector for the ith subject, and is fitted
  //against time with an intercept and slope, which in turn become the
  //estimates of Ui and Vi as the random effects. The estimation is based
  //on the formula from Crowder and Hand (1991)

  void   covmat(double*, double***, double*, int*, int*, int);
  void   covmat_1person(int, int, double*, double*, double**, int, int*,
                        int*, int);
  void   svdcmp(double**, double*, double**, int, int);
  void   A_mult_B(double**, double**, double**, int, int, int);
  void   A_mult_Bt(double**, double**, double**, int, int, int);
  void   deleteMatrix(double**, int);
  void   A_mult_v(double**, double*, double*, int, int);
  void   set_A_to_a(double**, double, int, int);
  double **myMatrix(int, int);
  double *w = (double *) malloc((unsigned) (Mntime+1)*sizeof(double));
  //new double[Mntime+1];
  double **B1 = myMatrix(Mntime, Mntime);
  double **D = myMatrix(Mntime, Mntime);
  double **covM1 = myMatrix(Mntime, Mntime);
  double **covM2 = myMatrix(Mntime, Mntime);
  double variance(double*, int);
  double fxx, fx;
  int    i, j, k, kk;

  for(i=1, k=1, kk=0; i<=n_sub; i++){
    //covmat_1person(i, ntime[i], in_covariate[n_covs+1], sigma, covM1, kk,
    covmat_1person(i, ntime[i], in_covariate[n_covs+2], sigma, covM1, kk,
                   var_inclusion, ntime, n_sub);
    /*add one for the covariance column */
    //if(in_covariate[n_covs][k] == in_covariate[n_covs][k+ntime[i]-1]){
    if(in_covariate[n_covs+1][k] == in_covariate[n_covs+1][k+ntime[i]-1]){
      /*add one for the covariance column */
      Ui[i]=resid[i][1];
      Vi[1][i]=0.;
      Vi[2][i]=0.;
      k += ntime[i];
    }
    else{
      svdcmp(covM1, w, covM2, ntime[i], ntime[i]);
      set_A_to_a(D, 0., ntime[i], ntime[i]);
      for(j=1; j<=ntime[i]; j++){
        D[j][j] = 1./w[j];
      }
      A_mult_Bt(covM1, covM2, B1, ntime[i], ntime[i], ntime[i]);
      A_mult_B(covM1, D, B1, ntime[i], ntime[i], ntime[i]);
      A_mult_Bt(B1, covM2, covM1, ntime[i], ntime[i], ntime[i]);
      for(j=1; j<=ntime[i]; j++) //D[1][j] = in_covariate[n_covs+1][kk+j];
        D[1][j] = in_covariate[n_covs+2][kk+j];
      /*add one for the covariance column */
      A_mult_v(covM1, D[1], w, ntime[i], ntime[i]);
      Ui[i]=0.;
      Vi[1][i]=0.;
      Vi[2][i]=0.;
      for(j=1; j<=ntime[i]; j++){
        Ui[i] += var_inclusion[0]*sigma[1]*w[j];
        Vi[1][i] += var_inclusion[1]*sigma[2]*
          //sqrt(in_covariate[n_covs][k])*w[j];
          sqrt(in_covariate[n_covs+1][k])*w[j];
        /*add one for the covariance column */
        //Vi[2][i] += var_inclusion[2]*sigma[3]*in_covariate[n_covs][k]*w[j];
        Vi[2][i] += var_inclusion[2]*sigma[3]*in_covariate[n_covs+1][k]*w[j];
        /*add one for the covariance column */
        k++;
      }
    }
    kk += ntime[i];
  }
  fx=0.;
  fxx=0.;
  for(i=1, k=1; i<=n_sub; i++){
    for(j=1; j<=ntime[i]; j++){
      Eij[i][j] = resid[i][j] - Ui[i]*var_inclusion[0]
      //	- sqrt(in_covariate[n_covs][k])*Vi[1][i]*var_inclusion[1]
      //   - in_covariate[n_covs][k]*Vi[2][i]*var_inclusion[2];
      - sqrt(in_covariate[n_covs+1][k])*Vi[1][i]*var_inclusion[1]
      - in_covariate[n_covs+1][k]*Vi[2][i]*var_inclusion[2];
      /*add one for the covariance column */
      k++;
      fxx += Eij[i][j]*Eij[i][j];
      fx += Eij[i][j];
    }
  }

  deleteMatrix(D, Mntime);
  deleteMatrix(B1, Mntime);
  deleteMatrix(covM1, Mntime);
  deleteMatrix(covM2, Mntime);
  free((char*) w);
}

void output_decp_residuals(char *datafile, double *Ui, double **Vi,
                           double **Eij, double **in_response, double **oriR1, int *var_inclusion,
                           int *ntime, int n_sub)
{
  char outfile[100];
  char outfile1[100];
  int ind=0;
  int i=0, j, j1;
  void trans1(char*);

  trans1(datafile);
  do{
    if(datafile[ind]=='.')
      i=0;
    else if((datafile[ind]!='\\')&&(datafile[ind]!='/'))
      outfile[i++] = datafile[ind];
  } while((datafile[ind]!='\0')&&(datafile[ind]!='\\')
            &&(datafile[ind++]!='/'));
  outfile[i--]='\0';
  trans1(outfile);

  strcpy(outfile1, outfile);
  strcat(outfile1,".s");
  strcat(outfile,".uve");

  ofstream residFile(outfile, ios::out);
  if (!residFile){
    //cerr << "cannot open file: " << outfile << endl;
    //exit( -1 );
  }
  if(var_inclusion[0]){
    for(j=1; j<=n_sub; j++){
      residFile << Ui[j] << " ";
      if(j %5 == 0) residFile << endl;
    }
    residFile << endl;
  }
  if(var_inclusion[1]){
    for(j=1; j<=n_sub; j++){
      residFile << Vi[1][j] << " ";
      if(j %5 == 0) residFile << endl;
    }
    residFile << endl;
  }
  if(var_inclusion[2]){
    for(j=1; j<=n_sub; j++){
      residFile << Vi[2][j] << " ";
      if(j %5 == 0) residFile << endl;
    }
    residFile << endl;
  }
  for(j=1; j<=n_sub; j++) for(j1=1; j1<=ntime[j]; j1++){
    residFile << in_response[j][j1]-oriR1[j][j1];
    residFile << " " << Eij[j][j1] << " ";
    if(j %3 == 0) residFile << endl;
  }
  residFile << endl;
  ofstream SplusFile(outfile1, ios::out);
  if (!SplusFile){
    //cerr << "cannot open file: " << outfile1 << endl;
    //exit( -1 );
  }
  SplusFile << "all1 <- scan(\"" << outfile << "\")" << endl;
    j=1;
    if(var_inclusion[0]){
      SplusFile << "ui1 <- all1[" << j << ":" << j+n_sub-1 << "]" << endl;
      SplusFile << "plot(ui1)" << endl;
      SplusFile << "hist(ui1)" << endl;
      j += n_sub;
    }
    if(var_inclusion[1]){
      SplusFile << "vi1 <- all1[" << j << ":" << j+n_sub-1 << "]" << endl;
      SplusFile << "plot(vi1)" << endl;
      SplusFile << "hist(vi1)" << endl;
      j += n_sub;
    }
    if(var_inclusion[2]){
      SplusFile << "vi2 <- all1[" << j << ":" << j+n_sub-1 << "]" << endl;
      SplusFile << "plot(vi2)" << endl;
      SplusFile << "hist(vi2)" << endl;
      j += n_sub;
    }
    SplusFile << "ye <- matrix(all1[-c(1:" << j-1 << ")], ncol=2,";
    SplusFile << " byrow=T)" << endl;
    SplusFile << "plot(ye)" << endl;
    SplusFile << "hist(ye[,2])" << endl;
}

void  Final_formula(int nterms, char* outfile, int* Fterms,
                    NumericVector beta, NumericMatrix varList)
{
  ofstream outFile(outfile, ios::app);
  if(!outFile){
    //cerr << "cannot open " << outfile << endl;
    //exit( -1 );
  }
  int      i, i1, j, k;

  outFile << "Final model in TeX format:\n";
  outFile << beta[1];
  for(i1=2; i1<=nterms; i1++){
    i = Fterms[i1];
    if(beta[i1] <= 0.) outFile << " " << beta[i1];
    else outFile << " + " << beta[i1];
    k = (int) varList(i,0);
    for(j=1; j<=k; j++){
      if(varList(i,3*j-1)> 0.){
        outFile << "(x_{" << varList(i,3*j-2)<< "}";
        if(varList(i, 3*j) > 0)
          outFile << "-" << varList(i,3*j)<< ")^+";
        else outFile << "+" << -varList(i, 3*j) << ")^+";
      }
      else
        outFile << "x_{" << varList(i,3*j-2)<< "}";
    }
    outFile << endl;
  }
}

double minization(double **h1, double *q1,int p, double (*find_lkh)(double*, double*,
                                                         int, double**, int*, int*, int, int), double *psum, double *ptry,
                                                         int method1, double *sigma,
                                                         int Mntime, double **in_covariate, int *var_inclusion,
                                                         int *ntime, int n_sub, int n_covs)
{
  void amoeba(double**, double*, int, double, double (*)(double*, double*, int, double**),
              int*, int, double*, double*, int, double**);
  void powell(double*, double**, int, double, int*, double*,
              double (*)(double*, double*, int, double**, int*,
                      int*, int, int), double*, int, double**, int*, int*, int, int);
  int   i, j, iter;//NMAX=200;
  double ftol=0.00001, fret;

  /*
  for(i=1; i<=p; i++)
  h1[1][i] = q1[i];
  for(i=2; i<=p+1; i++){
  for(j=1; j<=p; j++)
  h1[i][j]=h1[1][j];
  h1[i][i-1] += 0.1;
  }
  for(i=1; i<=p+1; i++)
  q1[i] = find_lkh(h1[i], sigma, Mntime, in_covariate);
  printf("Initial function value %f\n", q1[1]);
  amoeba(h1, q1, p, ftol, find_lkh , &iter, NMAX, psum, ptry, sigma, Mntime, in_covariate);
  fret = q1[1];
  for(i=1; i<=p; i++) q1[i] = h1[1][i];
  */
  for(i=1; i<=p; i++){
    for(j=1; j<=p; j++){
      h1[i][j]=0.;
    }
    h1[i][i]=1.;
  }
  powell(q1,h1,p, ftol, &iter, &fret, find_lkh, sigma, Mntime, in_covariate, var_inclusion,
         ntime, n_sub, n_covs);
  return fret;
}

void amoeba(double **p,double *y, int ndim, double ftol,
            double (*funk)(double*),
            int *nfunk, int NMAX, double *psum, double *ptry)
{
  /*y[ndim+1] is initially the starting point and function value*/
  /*p[ndim+1][ndim] is initally defined from y*/
  /*p then hosts the output points and y the values within FTOL of a
  minimal value*/
  int i,j,ilo,ihi,inhi,mpts=ndim+1;
  double ytry,ysave,sum,rtol;
  double amotry(double**, double*, double*, int, double (*)(double*),
                int, int*, double, double*);

  *nfunk=0;
  for (j=1;j<=ndim;j++)
    for (i=1,sum=0.0;i<=mpts;i++){
      sum += p[i][j];
      psum[j]=sum;
    }
    for (;;) {
      ilo=1;
      ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
      for (i=1;i<=mpts;i++) {
        if (y[i] < y[ilo]) ilo=i;
        if (y[i] > y[ihi]) {
          inhi=ihi;
          ihi=i;
        } else if (y[i] > y[inhi])
          if (i != ihi) inhi=i;
      }
      rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
      if (rtol < ftol) break;
      if (*nfunk >= NMAX){
        Rprintf("Warning: maximum iterations reached\n");
        break;
      }
      ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,-1.0, ptry);
      if (ytry <= y[ilo])
        ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,2.0, ptry);
      else if (ytry >= y[inhi]) {
        ysave=y[ihi];
        ytry=amotry(p,y,psum,ndim,funk,ihi,nfunk,0.5, ptry);
        if (ytry >= ysave) {
          for (i=1;i<=mpts;i++) {
            if (i != ilo) {
              for (j=1;j<=ndim;j++) {
                psum[j]=0.5*(p[i][j]+p[ilo][j]);
                p[i][j]=psum[j];
              }
              y[i]=(*funk)(psum);
            }
          }
          *nfunk += ndim;
          for (j=1;j<=ndim;j++)
            for (i=1,sum=0.0;i<=mpts;i++){
              sum += p[i][j];
              psum[j]=sum;
            }
        }
      }
    }
}

double amotry(double **p, double *y, double *psum, int ndim,
              double (*funk)(double*),
              int ihi, int *nfunk, double fac, double *ptry)
{
  int j;
  double fac1,fac2,ytry;

  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=(*funk)(ptry);
  ++(*nfunk);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=1;j<=ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  return ytry;
}


void powell(double *p, double **xi, int n, double ftol, int *iter,
            double *fret, double (*func)(double*, double*, int, double**, int*,
                                  int*, int, int),
                                  double *sigma, int Mntime, double **in_covariate, int *var_inclusion,
                                  int *ntime, int n_sub, int n_covs)
{
  /*p[n] is initially the starting point*/
  /*xi[n][n] is the search directions*/
  /*p[n] then hosts the minimizer and fret the minimal value*/
  int i,ibig,j;
  double t,fptt,fp,del;
  double *pt,*ptt,*xit;
  void   linmin(double*, double*, int, double*, double (*)(double*, double*, int,
                                                        double**, int*, int*, int, int), double*, int, double**, int*, int*, int, int);
  int ITMAX = 200;

  pt = (double *) malloc((unsigned) (n+1)*sizeof(double));
  //new double[n+1];
  ptt = (double *) malloc((unsigned) (n+1)*sizeof(double));
  //new double[n+1];
  xit = (double *) malloc((unsigned) (n+1)*sizeof(double));
  //new double[n+1];
  *fret=(*func)(p, sigma, Mntime, in_covariate, var_inclusion, ntime, n_sub, n_covs);
  for (j=1;j<=n;j++) pt[j]=p[j];  //save the initial point
  for (*iter=1;;(*iter)++) {
    fp=(*fret);
    ibig=0;
    del=0.0;                     //Will be the biggest function decrease.
    for (i=1;i<=n;i++) { //In each iteration, loop over all directions.
      for (j=1;j<=n;j++) xit[j]=xi[j][i]; //Copy the direction,
      fptt=(*fret);
      linmin(p,xit,n,fret,func, sigma, Mntime, in_covariate, var_inclusion,
             ntime, n_sub, n_covs);//minimize along it,
      if (fabs(fptt-(*fret)) > del) {     // and record it if it is the
        del=fabs(fptt-(*fret));           // largest decrease so far.
        ibig=i;
      }
    }
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
      free((char*)xit);                      //Termination criterion
      free((char*)ptt);
      free((char*)pt);
      return;
    }
    if(*iter == ITMAX){
      //cout << "POWELL exceeding maximum iteration of " << ITMAX << endl;
      //Rcpp::stop("POWELL exceeding maximum iteration");
      //exit(-1);
    }
    for (j=1;j<=n;j++) { //Construct the extrapolated point and the average
      ptt[j]=2.0*p[j]-pt[j]; //direction moved. Save the old starting point.
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt, sigma, Mntime, in_covariate, var_inclusion,
          ntime, n_sub, n_covs);  //Function value at extrapolated point.
    if (fptt < fp) {    //One reason not to use new direction.
      t=2.0*(fp-2.0*(*fret)+fptt)*(fp-(*fret)-del)*(fp-(*fret)-del)-del*(fp-fptt)*(fp-fptt);
      if (t < 0.0) {    //Other reason no to use new direction.
        linmin(p,xit,n,fret,func, sigma, Mntime, in_covariate, var_inclusion,
               ntime, n_sub, n_covs);
        //Move the minimum of the new direction
        for (j=1;j<=n;j++) xi[j][ibig]=xit[j]; //and save the new direction.
      }
    }
  }
}

void linmin(double *p, double *xi, int n, double* fret,
            double (*func)(double*, double*, int, double**, int*, int*, int, int),
            double *sigma, int Mntime,
            double **in_covariate, int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  double brent(double, double, double, double (*)(double, double*, int, double**,
                                               int, double*, double*, int*, int*, int, int),
                                               double, double*, double*, int, double**, int, double*, double*, int*,
                                               int*, int, int);
  double f1dim(double, double*, int, double**, int, double*, double*, int*, int*, int, int);
  void   mnbrak(double*, double*, double*, double*,
                double*, double*, double (*)(double, double*, int, double**, int, double*,
                                          double*, int*, int*, int, int), double*, int, double**, int, double*, double*,
                                          int*, int*, int, int);
  double  *pcom, *xicom;
  double TOL = 2.0e-4;

  pcom = (double *) malloc((unsigned) (n+1)*sizeof(double));
  //new double[n+1];
  xicom = (double *) malloc((unsigned) (n+1)*sizeof(double));
  //new double[n+1];
  nrfunc=func;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  bx=2.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim, sigma, Mntime, in_covariate, n, pcom,
         xicom, var_inclusion, ntime, n_sub, n_covs);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin, sigma, Mntime, in_covariate, n, pcom,
              xicom, var_inclusion, ntime, n_sub, n_covs);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free((char*) xicom);
  free((char*) pcom);
}

double brent(double ax, double bx, double cx, double (*f)(double, double*, int,
                                                      double**, int, double*, double*, int*, int*, int, int), double tol, double *xmin,
                                                      double *sigma, int Mntime, double **in_covariate, int ncom,
                                                      double *pcom, double *xicom, int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  int iter;
  double a,b,d = 0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  double CGOLD = 0.3819660;
  double SIGN(double, double);
  int    ITMAX = 200;
  void   SHFT(double&, double&, double&, double);

  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
            ntime, n_sub, n_covs);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+1.0e-10);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
        d=p/q;
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
        ntime, n_sub, n_covs);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }

  *xmin=x;
  return fx;
}

double f1dim( double x, double *sigma, int Mntime, double **in_covariate, int ncom,
              double *pcom, double *xicom, int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  int j;
  double f,*xt;

  xt = (double *) malloc((unsigned) (ncom+1)*sizeof(double));
  //new double[ncom+1];
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt, sigma, Mntime, in_covariate, var_inclusion, ntime, n_sub, n_covs);
  free((char*) xt);
  return f;
}

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
            double *fc, double (*func)(double, double*, int, double**, int, double*,
                                double*, int*, int*, int, int),
                                double *sigma, int Mntime, double **in_covariate, int ncom, double *pcom,
                                double *xicom, int *var_inclusion, int *ntime, int n_sub, int n_covs)
{
  double ulim,u,r,q,fu,dum=0.;
  double GOLD = 1.618034;
  double MAX(double, double);
  double SIGN(double, double);
  void SHFT(double&, double&, double&, double);

  *fa=(*func)(*ax, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
       ntime, n_sub, n_covs);
  *fb=(*func)(*bx, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
       ntime, n_sub, n_covs);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum);
    SHFT(dum,*fb,*fa,dum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
       ntime, n_sub, n_covs);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(MAX(fabs(q-r),1.0e-20),q-r));
    ulim=(*bx)+100.0*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
          ntime, n_sub, n_covs);
      if (fu < *fc) {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
          ntime, n_sub, n_covs);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
          ntime, n_sub, n_covs);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
        SHFT(*fb,*fc,fu,(*func)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom,
                         var_inclusion, ntime, n_sub, n_covs));
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
          ntime, n_sub, n_covs);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u, sigma, Mntime, in_covariate, ncom, pcom, xicom, var_inclusion,
          ntime, n_sub, n_covs);
    }
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
  }
}

void sort_masal(int n, double *ra, int *rb)
{
  int l,j,ir,i;
  int rrb;
  double rra;

  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
      rrb=rb[l];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      ra[ir]=ra[1];
      rb[ir]=rb[1];
      if (--ir == 1) {
        ra[1]=rra;
        rb[1]=rrb;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir)	{
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
  }
}

double **myMatrix(int nrow, int ncol){
  int     i;
  double  **tmp_matrix;

  tmp_matrix = (double **) malloc((unsigned) (nrow+1)*sizeof(double*));
  //new double*[nrow+1];
  for(i=1; i<=nrow; i++){
    tmp_matrix[i] = (double *) malloc((unsigned) (ncol+1)*sizeof(double));
    //new double[ncol+1];
  }
  return tmp_matrix;
}

double **myMatrix2(int nrow, int *ncol){
  int     i;
  double  **tmp_matrix;

  tmp_matrix = (double **) malloc((unsigned) (nrow+1)*sizeof(double*));
  //new double*[nrow+1];
  for(i=1; i<=nrow; i++){
    tmp_matrix[i] = (double *) malloc((unsigned) (ncol[i]+1)*sizeof(double));
    //new double[ncol[i]+1];
  }
  return tmp_matrix;
}

double ***myMatrix3(int nrow, int *ncol){
  int     i, j;
  double  ***tmp_matrix;

  tmp_matrix = (double ***) malloc((unsigned) (nrow+1)*sizeof(double**));
  //new double**[nrow+1];
  for(i=1; i<=nrow; i++){
    tmp_matrix[i] = (double **) malloc((unsigned) (ncol[i]+1)*sizeof(double*));
    //new double*[ncol[i]+1];
    for(j=1; j<=ncol[i]; j++)
      tmp_matrix[i][j] = (double *) malloc((unsigned) (ncol[i]+1)*sizeof(double));
    //new double[ncol[i]+1];
  }
  return tmp_matrix;
}

int **myImatrix(int nrow, int ncol){
  int     i;
  int   **tmp_matrix;

  tmp_matrix = (int **) malloc((unsigned) (nrow+1)*sizeof(int*));
  //new int*[nrow+1];
  for(i=1; i<=nrow; i++){
    tmp_matrix[i] = (int *) malloc((unsigned) (ncol+1)*sizeof(int));
    //new int[ncol+1];
  }
  return tmp_matrix;
}

int **myImatrix2(int nrow, int *ncol){
  int     i;
  int   **tmp_matrix;

  tmp_matrix = (int **) malloc((unsigned) (nrow+1)*sizeof(int*));
  //new int*[nrow+1];
  for(i=1; i<=nrow; i++){
    tmp_matrix[i] = (int *) malloc((unsigned) (ncol[i]+1)*sizeof(int));
    //new int[ncol[i]+1];
  }
  return tmp_matrix;
}

void A_dot_B3(double **mA, double **mB, double **mC, int nrow, int *ncol)
{
  int  i, j;

  for(i=1; i<=nrow; i++) for(j=1; j<=ncol[i]; j++)
    mC[i][j] = mA[i][j]*mB[i][j];
}

void copy_A_to_B(double **mA, double **mB, int nrow, int ncol)
{
  int i, j;
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol; j++)
    mB[i][j] = mA[i][j];
}

void copy_A_to_B3(double **mA, double **mB, int nrow, int *ncol)
{
  int i, j;
  for(i=1; i<=nrow; i++) for(j=1; j<=ncol[i]; j++)
    mB[i][j] = mA[i][j];
}

double matVecNorm(double **mA, double *v1, int nrow, int *ncol)
{
  int    i, j, k;
  double tmp=0.;

  for(i=1, k=1; i<=nrow; i++)
    for(j=1; j<=ncol[i]; j++, k++)
      tmp += mA[i][j]*v1[k];
  return (tmp*tmp);
}

void mat_mult_mat3(double **A1, int rowHi, double **B1, double *b2,
                   int nrow, int *ncol)
{
  register int i, j, jn, k;
  void mat_mult_mat3(double**, int, double**, double*, int, int*);

  for (i=1; i<= rowHi; i++){
    b2[i] = 0.;
    for (j=1, jn=1; j<=nrow; j++)
      for(k=1; k<=ncol[j]; k++, jn++)
        b2[i] += A1[i][jn]*B1[j][k];
  }
}

void A_mult_B(double **A1, double **B1, double **C1, int rowA, int colA,
              int colB)
{
  int i, j, k;

  for (i=1; i<= rowA; i++)
    for (j=1; j<= colB; j++) {
      C1[i][j] = 0.;
      for(k=1; k<=colA; k++)
        C1[i][j] += A1[i][k]*B1[k][j];
    }
}

void At_mult_B(double **A1, double **B1, double **C1, int rowA, int colA,
               int colB)
{
  int i, j, k;

  for (i=1; i<= colA; i++)
    for (j=1; j<= colB; j++) {
      C1[i][j] = 0.;
      for(k=1; k<=rowA; k++)
        C1[i][j] += A1[k][i]*B1[k][j];
    }
}

void A_mult_Bt(double **A1, double **B1, double **C1, int rowA, int colA,
               int rowB)
{
  int i, j, k;

  for (i=1; i<= rowA; i++)
    for (j=1; j<= rowB; j++) {
      C1[i][j] = 0.;
      for(k=1; k<=colA; k++)
        C1[i][j] += A1[i][k]*B1[j][k];
    }
}

void A_mult_v(double **mA, double *b1, double *b2, int rowHi, int colHi)
{
  int i, j;

  for (i=1; i<= rowHi; i++){
    b2[i] = 0.;
    for (j=1; j<= colHi; j++)
      b2[i] += mA[i][j]*b1[j];
  }
}
void A_mult_v1(double **mA, double *b1, NumericVector b2, int rowHi, int colHi)
{
  int i, j;

  for (i=1; i<= rowHi; i++){
    b2[i] = 0.;
    for (j=1; j<= colHi; j++)
      b2[i] += mA[i][j]*b1[j];
  }
}

double a_dot_b(double *va, double *vb, int n1)
{
  int    i;
  double tmp=0.;

  for(i=1; i<=n1; i++) tmp += va[i]*vb[i];
  return tmp;
}

double *myVector(int nrow){
  double  *tmp_vec;

  tmp_vec = (double *) malloc((unsigned) (nrow+1)*sizeof(double));
  //new double[nrow+1];
  return tmp_vec;
}

void set_A_to_a(double **mA, double a, int nrow, int ncol){
  int i, j;

  for(i=1; i<=nrow; i++)
    for(j=1; j<=ncol; j++) mA[i][j]=a;
}

void set_A_to_a3(double **mA, double a, int nrow, int *ncol){
  int i, j;

  for(i=1; i<=nrow; i++)
    for(j=1; j<=ncol[i]; j++) mA[i][j]=a;
}

void set_v_to_a(double *ma, double a, int n1){
  int  i;
  for(i=1; i<=n1; i++) ma[i]=a;
}

void copy_va_to_vb(double *va, double *vb, int n1){
  int  i;
  for(i=1; i<=n1; i++) vb[i]=va[i];
}

double variance(double *va, int n1){
  int  i;
  double tmp=0., tmp1;
  double mean(double*, int);

  for(i=1; i<=n1; i++) tmp += va[i]*va[i];
  tmp1 = mean(va, n1);
  tmp -= tmp1*tmp1*n1;
  return tmp/(n1-1.);
}

double mean(double *va, int n1){
  int  i;
  double tmp=0.;

  for(i=1; i<=n1; i++) tmp += va[i];
  tmp /= n1;
  return tmp;
}

void deleteMatrix(double **mA, int nn)
{
  int i;
  for(i=nn; i>=1; i--) free((char*) mA[i]);
  free((char*) mA);
}

void deleteMatrix3(double ***mA, int n1, int *n2)
{
  int i, j;

  for(i=n1; i>=1; i--){
    for(j=n2[i]; j>=1; j--) free((char*)mA[i][j]);
    free((char*)mA[i]);
  }
  free((char*) mA);
}

void svdcmp(double **a, double *w, double **v, int m, int n)
{
  int    flag,i,its,j,jj,k,l=1,nm=1;
  double c,f,h,s,x,y,z;
  double anorm=0.0,g=0.0,scale=0.0;
  double *rv1 = myVector(n);
  double MAX(double, double);
  double SIGN(double, double);
  double PYTHAG(double, double);

  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        if (i != n) {
          for (j=l;j<=n;j++) {
            for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
            f=s/h;
            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
          }
        }
        for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<=n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
        if (i != m) {
          for (j=l;j<=m;j++) {
            for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
          }
        }
        for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
        for (j=l;j<=n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n;i>=1;i--) {
    l=i+1;
    g=w[i];
    if (i < n)
      for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n) {
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
          f=(s/a[i][i])*g;
          for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
        }
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    }
    else {
      for (j=i;j<=m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
        nm=l-1;
        if (fabs(rv1[l])+anorm == anorm) {
          flag=0;
          break;
        }
        if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i<=k;i++) {
          f=s*rv1[i];
          if (fabs(f)+anorm != anorm) {
            g=w[i];
            h=PYTHAG(f,g);
            w[i]=h;
            h=1.0/h;
            c=g*h;
            s=(-f*h);
            for (j=1;j<=m;j++) {
              y=a[j][nm];
              z=a[j][i];
              a[j][nm]=y*c+z*s;
              a[j][i]=z*c-y*s;
            }
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
        }
        break;
      }
      /*if(its==30){
        Rcpp::stop("reached 30 iterations in svdcmp");
        //cout << "Warning: reached 30 iterations in svdcmp" << endl;
        // exit(-1);
      }*/
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=PYTHAG(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y=y*c;
        for (jj=1;jj<=n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=PYTHAG(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=(c*g)+(s*y);
        x=(c*y)-(s*g);
        for (jj=1;jj<=m;jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free((char*) rv1);
}

double MAX(double a, double b){
  return (a > b) ? a : b;
}

double SIGN(double a, double b){
  return (b > 0.0) ? a : -a;
}

void SHFT(double &a, double &b, double &c, double d){
  a=b;
  b=c;
  c=d;
}

double PYTHAG(double a, double b){
  double at,bt,ct;
  return ((at=fabs(a)) > (bt=fabs(b)) ?
            (ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0));
}


void input(double *data, int nT, int n_covs, int *n_sub, double **in_covariate)
{
  int   i, j=1, k, kk;
  int   tmpID;

  //in_covariate[n_covs+2][j] = 1.;
  in_covariate[n_covs+3][j] = 1.;
  /*add one for the covariance column */
  tmpID = (int) data[0];
  //for(kk=1; kk<=n_covs+1; kk++){
  for(kk=1; kk<=n_covs+2; kk++){
    /*add one for the covariance column */
    in_covariate[kk][1] = data[kk];
  }
  for(i=2; i<=nT; i++){
    if(tmpID == (int) data[kk])
      //in_covariate[n_covs+2][j] +=1.;
      in_covariate[n_covs+3][j] +=1.;
    /*add one for the covariance column */
    else{
      j++;
      //in_covariate[n_covs+2][j]=1.;
      in_covariate[n_covs+3][j]=1.;
      /*add one for the covariance column */
      tmpID = (int) data[kk];
    }
    kk++;
    //for(k=1; k<=n_covs+1; k++){
    for(k=1; k<=n_covs+2; k++){
      /*add one for the covariance column */
      in_covariate[k][i] = data[kk++];
    }
  }
  *n_sub = j;
}

void sorting_data(int *ntime, int *Mntime, int nT, int n_covs, int n_sub,
                  double **in_covariate, int &regression)
{
  int   i, j=1, k, kk;
  int   *order0;
  void  sort_masal(int, double*, int*);
  double *tmpT;

  regression = 1;
  for(i=1; i<=n_sub; i++){
    //ntime[i] = (int) in_covariate[n_covs+2][i];
    ntime[i] = (int) in_covariate[n_covs+3][i];
    /*add one for the covariance column */
    if(ntime[i] > 1)
      regression=0;
  }
  *Mntime = 0;
  for(i=1; i<=n_sub; i++)
    if(ntime[i] > (*Mntime) ) (*Mntime) = ntime[i];
    order0 = (int *) malloc((unsigned) (nT+1)*sizeof(int));
    tmpT = (double *) malloc((unsigned) (nT+1)*sizeof(double));
    for(i=1, kk=1; i<=n_sub; i++){
      for(j=1; j<=ntime[i]; j++){
        tmpT[j] = in_covariate[n_covs][kk];
        order0[j]=j;
        kk++;
      }
      if(ntime[i] > 1){
        sort_masal(ntime[i], tmpT, order0);
        //for(k=1; k<=n_covs+1; k++){
        for(k=1; k<=n_covs+2; k++){
          /*add one for the covariance column */
          for(j=1; j<=ntime[i]; j++)
            tmpT[j]=in_covariate[k][kk-ntime[i]+j-1];
          for(j=1; j<=ntime[i]; j++)
            in_covariate[k][kk-ntime[i]+j-1] = tmpT[order0[j]];
        }
      }
    }
    free((char*)order0);
    free((char*) tmpT);
}

double *scan_datafile(DataFrame dataset, int &nT, int &n_covs){
  int    i,j;
  double *data;
  int    ncol, nsub;

  NumericVector coltemp;
  coltemp = dataset[0];
  nT = coltemp.size();
  ncol = dataset.length();
  nsub = nT*ncol;
  //cout << nsub << endl;
  data = (double*) malloc ((unsigned) nsub*sizeof(double));
  //n_covs=ncol-2;                    /*save # of covariates*/
  n_covs=ncol-3;                      /*excludes the column for covariance 3/3/04 */
  for(i=0;i<ncol;i++)
  {
    coltemp = dataset[i];
    for(j=0;j<nT;j++)
    {
      data[i+j*ncol] = coltemp[j];
    }
  }
  return data;
}
