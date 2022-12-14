#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "head.h"
#include "R_ext/Print.h"

#define MNNODES 200
using namespace Rcpp;
using namespace std;

void   r_nrerror(const char *);
void   free_matrix(double **,int);
void   free_dmatrix(double **, int);
void   free();
void   exit(int);
int    **imatrix2(int, int *);
int    *r_ivector(int);
double *r_dvector(int);
double ran1(int *);
double **dmatrix();
double **dmatrix2(int , int *);
double *ddvector(int);
int    *icopy_point_array(int *, int);
double *dcopy_point_array(double *, int);
int    max_index(int *, int); // the index for the maximum class
int    exist_in_randomlist(int, int, int *); // decide whether i_temp in the first i-th of index_random
int    sum(int *, int); // sum of an array
void   sort_array(int *, int *, int *, int);
// sort array data into sort_data, and put the initial index for sort_data in index
void   free_linklist(struct linklist *);

#define REP 100
#define true 1
#define false 0

struct linklist
{
  int col;
  int row;
  struct linklist *next;
};

int    r_n_obs;             /*sample size or number of observations or cases*/
int    r_n_meas;            /*number of measurements taken for each case*/
int    r_n_class;           /*number of classes in y*/
double **r_X;             //X[i][j]: jth measurement for ith case; oriX: record the original data
                        // if ordinal, then the X[i][j] is the true value
                        // if nominal, then X[i][j] is the the index in the ori_catg
int    *r_y, *vartype; int r_beforeprune=0;
int    *r_catg; // numbers of catogories for variablefor for nominal
int    *nodeclass;

int    *odds;

// [[Rcpp::export]]
List cctree(DataFrame dataset, IntegerVector colname, int srule, double pvalue)
{
  IntegerMatrix grow_large_tree(double **,int ,int *,IntegerVector ,NumericVector ,
                                IntegerVector ,IntegerVector ,int *,int *, int *,
                                int *, int , double , NumericVector );
  void   find_dist(int , int *, int *);
  void   find_dist1(int , int *, IntegerVector);
  void   input_rtree(DataFrame , IntegerVector );
  void   trans1( char *);
  int    keyin2(int, char *);
  int    nnd, *npl;
  int    *ncases, *dist;
  double *dp_temp;

  input_rtree(dataset, colname);

  IntegerVector spv(MNNODES), dt(MNNODES), pt(MNNODES);
  NumericVector spvl(MNNODES), chi(MNNODES);

  npl = r_ivector(r_n_obs);
  ncases = r_ivector(MNNODES);
  dist = r_ivector(r_n_class);
  nodeclass = r_ivector(MNNODES);
  dp_temp = r_dvector(r_n_meas);
  IntegerMatrix final_counts_matrix = grow_large_tree(r_X, r_n_obs, &nnd, spv, spvl, pt, dt,
                                                      npl, ncases, dist, r_catg, srule, pvalue, chi);

  //	free
  free((char*)dp_temp); free((char*)npl); free((char*)ncases);
  free((char*)r_catg); free((char*)vartype); free((char*)nodeclass);
  free((char*)r_y);free_matrix(r_X, 2*r_n_meas);

  return List::create(
    Named("nnd") = nnd,
    Named("final_counts") = final_counts_matrix,
    Named("spv") = spv,
    Named("spvl") = spvl,
    Named("dt") = dt,
    Named("pt") = pt,
    Named("chi") = chi,
    Named("varcatg") = colname
  );
}

struct linklist * l_alloc(void )
{
  struct linklist *v;

  v=(struct linklist *)malloc(sizeof(struct linklist));
  if(!v)
  {
    r_nrerror("allocation failure in l_alloc()");
  }
  return v;
}

struct linklist * addmiss(struct linklist *p1,int cl,int rw)
{
  struct linklist *temp;
  temp = p1;
  if(temp == NULL)  /* first element */
  {
    p1 = l_alloc();
    p1->col = cl;
    p1->row=rw;
    p1->next=NULL;
  }
  else
  {
    while(temp->next!=NULL)
      temp = temp->next;
    temp->next = l_alloc();
    temp->next->col = cl;
    temp->next->row = rw;
    temp->next->next= NULL;
  }
  return (p1);
}

IntegerMatrix grow_large_tree(double **xx,int n,int *nnd,IntegerVector spv,NumericVector spvl,
                              IntegerVector pt,IntegerVector dt,int *npl,int *ncases, int *dist,
                              int *r_catg, int srule, double pvalue, NumericVector chi)
{
  int   t, *bodds;
  void  initialize_tree(int , int* , IntegerVector , IntegerVector , NumericVector );
  void  prune_tree(IntegerVector ,int *, int *, int *);
  void  split_a_node(int ,double **,int ,int *,IntegerVector , NumericVector ,
                     IntegerVector , IntegerVector , int *,int *,int *,
                     int , NumericVector);
  void  split_a_node1(int, double **,int ,int *,int *,double *,int *,
                      int *,int *,int *,int *,int *,int );
  void  print_tree(int, int *, double *, int *, int *, int *,
                   int *, int *, FILE *);
  void  print_draw_tree(int,int *,int *,int *,double *, int *, int **, int);
  void  reorganize_tree(int *, IntegerVector , IntegerVector , int *, int *,
                        IntegerVector , NumericVector , IntegerMatrix , int *);

  void  pause(char *);
  void  ps_format(FILE *);
  void  end_ps(FILE *, int);
  int   *chi1, *neworder;

  initialize_tree(n, npl, dt, spv, spvl);
  chi1=r_ivector(MNNODES);
  bodds=r_ivector(2*r_n_class);
  odds=r_ivector(2*r_n_class);
  neworder = r_ivector(MNNODES);
  for(t=0, *nnd=1; t < *nnd; t++)
  {
    split_a_node(t, xx, n, nnd, spv, spvl, pt, dt,
                 npl, ncases, bodds, srule, chi);
  }
  IntegerMatrix final_counts(*nnd, r_n_class);

  for(t=0; t < *nnd; t++)
  {
    chi1[t] = 1*(chi[t] >= R::qchisq(1 - pvalue, r_n_class - 1, TRUE, FALSE));
  }
  prune_tree(dt, npl, chi1, nnd);
  reorganize_tree(nnd, dt, pt, npl, dist, spv, spvl, final_counts,  neworder);
  for(t=0; t < *nnd; t++)
  {
    ncases[t] = ncases[neworder[t]];
    chi[t] = chi[neworder[t]];
  }

  free((char*)chi1);
  free((char*)neworder); 	free((char*)bodds);free((char*)odds);
  return final_counts;
}

void split_a_node(int t,double **xx,int n,int *nnd,IntegerVector spv, NumericVector spvl,
                  IntegerVector pt,IntegerVector dt, int *npl,int *ncases,int *bodds,
                  int srule, NumericVector chi)
{
  int    i, ii, j, k, kk;
  int    *cases_in_t,  n_cases;
  int    yes;
  double min_mr;
  double thresh, it;
  int    *find_cases_in_t(int,int,int *,int *);
  double ord_best_split(double *,int *,int,double *, int);
  double catg_best_split(double *,int,int *,int,double *,int);
  void   split_via(int ,double **,int *,int ,
                   int *,int *, NumericVector , int , int );
  int    check_if_needed_split(double **,int *,int);
  double find_chisq(int *,int);

  cases_in_t=find_cases_in_t(t, n, npl, &n_cases);
  if(n_cases <= 0)
  {
    j=pt[t];
    Rprintf("The split using variable %d at a value of %f results in",
           spv[j]+1, spvl[j]);
    Rprintf(" an empty node\n");
    Rprintf("Press Enter to end this program!");
    getchar();
    //exit(0);
  }
  ncases[t]=n_cases;

  if(*nnd < MNNODES)
  {
    yes=check_if_needed_split(xx, cases_in_t, n_cases);
    if(yes)
    {
      min_mr=1.*r_n_class;
      for(i=0; i<r_n_meas; i++)
      {
        if(r_catg[i] > 1)
          it=catg_best_split(xx[i], r_catg[i], cases_in_t, n_cases, &thresh, srule);
        else
          it=ord_best_split(xx[i], cases_in_t, n_cases, &thresh, srule);
        if(it < min_mr)
        {
          min_mr=it;
          spv[t]=i;
          spvl[t]=thresh;
        }
      } // for(i_temp=0; i_temp<i_numRanViar; i_temp++)
      if(min_mr < 1.*r_n_class)
      {
        k=spv[t];
        kk=(int) spvl[t];
        /* kk is the criterion of split for categorical variable, it doesn't matter if
        kk is not consistant with spvl in ordinal variable because kk is referenced
        in split_via and it is only meanful for categorical variable C-Y*/
        split_via(t, xx, cases_in_t, n_cases, npl, nnd, spvl, k, kk);
        for(ii=0; ii<2*r_n_class; ii++) bodds[ii]=odds[ii];
        pt[*nnd]=t, pt[*nnd+1]=t;
        dt[t]= *nnd;
        *nnd += 2;
        chi[t] = find_chisq(bodds, ncases[t]);
      } // if(min_mr < 1.*n_class)
      else
      {
        dt[t]=0;
        chi[t] = 0.;
      }
    } // if(yes)
    else
    {
      dt[t]=0;
      chi[t] = 0.;
    } // if(yes) else
  } // if (*nnd < MNNODES)
  free((char*)cases_in_t);
}

int *find_cases_in_t(int t,int n,int * npl,int * k)
{
  int *tmpi, i;
  tmpi=r_ivector(n);
  for(i=0, *k=0; i<n; i++)
  {
    if(npl[i]==t)
    {
      tmpi[*k]=i;
      *k += 1;
    }
  }
  return tmpi;
}

double ord_best_split(double *x,int *cases_in_t,int n_cases,double *thresh, int srule)
{
  double *tmpx;
  double min_mr, tmp, imp,impurity(int *,int,int *,int,int,int);
  int    *cp1, *cp2, i, j, jobs=0, *tmpc;
  double prob(int , int ,int ,int *);
  int    k, *disx;
  int    *distinct_x(double *, int, int *);
  void   sort2(int ,double *, int *);

  cp1  = r_ivector(r_n_class);
  cp2  = r_ivector(r_n_class);
  tmpx = r_dvector(n_cases+1);
  tmpc = r_ivector(n_cases+1);
  for(i=1; i<=n_cases; i++)
  {
    j=cases_in_t[i-1];
    tmpx[i]=x[j];
    tmpc[i]=r_y[j];
  }
  sort2(n_cases, tmpx, tmpc);
  for(i=0; i<n_cases; i++)
  {
    tmpx[i]=tmpx[i+1];
    tmpc[i]=tmpc[i+1];
  }
  disx=distinct_x(tmpx, n_cases, &k);
  min_mr=1.*r_n_class;
  if(k>1)
  {
    for(j=0; j<r_n_class; j++)
    {
      cp1[j]=0;
      cp2[j]= (int) prob(j, 0, n_cases, tmpc);
    }
    for(i=0; i<k-1; i++)
    {
      for(j=0; j<r_n_class; j++)
      {
        tmp=(i>0)? prob(j, disx[i-1], disx[i], tmpc):
        prob(j, 0, disx[0], tmpc);
        cp1[j] += (int) tmp, cp2[j] -= (int) tmp;
        if(cp2[j]<=0) cp2[j]=0;
      }
      if(disx[i] > r_n_obs/100 && n_cases-disx[i] > r_n_obs/100)
      {
        imp=impurity(cp1, disx[i], cp2, n_cases-disx[i], n_cases, srule);
        if(imp < min_mr)
        {
          min_mr = imp;
          jobs = disx[i]-1;
        }
      }
    }
  }
  *thresh=tmpx[jobs];
  free((char*)tmpx), free((char*)tmpc), free((char*)disx);
  free((char*)cp1), free((char*)cp2);
  return min_mr;
}

int *distinct_x(double * tmpx,int n_cases,int *k)
{
  int   i, j, *disx;
  disx=r_ivector(n_cases);
  *k=1; disx[0]=0, j=0;
  for(i=0; i<n_cases; i++)
  {
    if(tmpx[i]==tmpx[j])
      disx[*k-1] +=1;
    else
    {
      *k += 1;
      disx[*k-1] = disx[*k-2]+1;
      j=disx[*k-2];
    }
  }
  return disx;
}

void initialize_tree(int n, int* npl, IntegerVector dt, IntegerVector spv, NumericVector spvl)
{
  int i;
  for(i=0; i<n; i++) npl[i]=0;
  for(i=0; (i<2*n) && (i<MNNODES); i++)
  {
    dt[i]=0;
    spv[i]=0;
    spvl[i]=0.;
  }
}

double prob(int i, int i1,int i2,int *tmpc)
{
  int j, k;
  for(j=0, k=i1; k<i2; k++)
  j += (tmpc[k]==i);
  return (double)j;
}

void sort2(int n,double *ra, int *rb)
{
  int l,j,ir,i, rrb;
  double rra;
  l=(n >> 1)+1;
  ir=n;
  for(;;)
  {
    if(l > 1)
    {
      rra=ra[--l];
      rrb=rb[l];
    }
    else
    {
      rra=ra[ir];
      rrb=rb[ir];
      ra[ir]=ra[1];
      rb[ir]=rb[1];
      if (--ir == 1)
      {
        ra[1]=rra;
        rb[1]=rrb;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir)
    {
      if (j < ir && ra[j] < ra[j+1]) ++j;
        if (rra < ra[j])
        {
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

int check_if_needed_split(double **xx,int *cases,int n_cases)
{
  int yes1=0, yes2=0, i, j, i0, k;
  if(n_cases==1) return 0;
  i0=cases[0];
  for(i=1; i<n_cases; i++)
  {
    j=cases[i];
    if(r_y[j] != r_y[i0]) yes1=1;
  }
  if(yes1)
  {
    for(i=0; i<r_n_meas; i++)
    {
      for(k=1; k<n_cases; k++)
      {
        j=cases[k];
        if(xx[i][j] != xx[i][i0]) yes2=1;
      }
    }
  }
  return yes2;
}

double catg_best_split(double *x, int nxcat, int *p, int size,
                       double *lcutc, int srule)
{
  double  imp, timp=1.*r_n_class;
  double  impurity(int *,int ,int *,int ,int ,int );

  int    i, j, k, l, m, ucat=0, ulim, isave;
  int    *xu, *counts, *cntl, *cntr;
  int    nc, nx, ny, nl, nr;
  int    *xcat, *ynum, *ycnt, *ycat;

  isave=1;
  xu=r_ivector(nxcat), counts=r_ivector(nxcat*r_n_class);
  cntr=r_ivector(r_n_class), cntl=r_ivector(r_n_class);
  xcat=r_ivector(nxcat), ycat=r_ivector(nxcat*r_n_class);
  ynum=r_ivector(nxcat), ycnt=r_ivector(nxcat*r_n_class);

  for(k=0;k<nxcat;k++)
  {
    xu[k]=xcat[k]=ynum[k]=0;
    for(i=0;i<r_n_class;i++) counts[k*r_n_class+i]=0;
  }
  for(i=0;i<r_n_class;i++)
    cntl[i]=0, cntr[i]=0;
  for(i=0;i<size;i++)
  {
    j=p[i];
    nx = (int) x[j];
    ny = r_y[j];  /*nxcat*n_class cell */
    xu[nx]=1;
    /*categorical y: counts for each */
    counts[nx*r_n_class + ny]++;
    cntr[ny]++;
  }
  ucat=0;
  for (i=0; i<nxcat; i++) if(xu[i]) xcat[ucat++]=i;
  if (ucat<=1) return(1.*r_n_class);
  for(k=0;k<ucat;k++)
  {
    for(m=0;m<r_n_class;m++)
    {
      nc=counts[xcat[k]*r_n_class + m];
      if(nc)
      {
        ycnt[k*r_n_class+ynum[k]]=nc;
        ycat[k*r_n_class+ynum[k]]=m;
        ynum[k]++;
      }
    }
  }
  ulim = (((unsigned)1)<<(ucat-1)) - 1;
  nr=size, nl=0;
  for(i=1; i<=ulim; i++)
  {
    l=i; k=0;
    while(!(l&1))
    {
      l>>=1; k++;
    }
    for(j=0; j<ynum[k]; j++)
    {
      nc = ycnt[k*r_n_class+j];
      m = ycat[k*r_n_class+j];
      if(!(l&2))
      {
        /* bit goes 0 -> 1 */
        cntl[m]+=nc; nl+=nc;
        cntr[m]-=nc; nr-=nc;
      }
      else
      {
        /* bit goes 1 -> 0 */
        cntl[m]-=nc; nl-=nc;
        cntr[m]+=nc; nr+=nc;
      }
    }
    if(nr > r_n_obs/100 && nl > r_n_obs/100)
    {
      imp=impurity(cntr, nr, cntl, nl, size, srule);
      if(imp < timp)
      {
        /* save minimp and best split */
        isave = i ;
        timp = imp;
      }
    }
  }
  *lcutc = (double) isave;
  free((char*)xu),  free((char*)counts);
  free((char*)cntl), free((char*)cntr);
  free((char*)ynum), free((char*)ycnt);
  free((char*)ycat), free((char*)xcat);
  return(timp);
}


void r_grayrep(int i,int *v,int n)
  /* outputs n-digit gray code for i */
{
  int j,b,ob=0;
  for(j=0;j<n;j++)
  {
    /* binary rep of i */
    v[j] = i&1;
    i >>= 1;
  }

  for(j=n-1;j>=0;j--)
  {
    /* convert to gray rep */
    b=v[j];
    v[j] = ob ? !b : b;
    ob=b;
  }
}


int split_cases(double *x, int *xu,int nxcat,int *p,int size)
/*
  x: xx[k]
  xu: ivextor(r_catg[k])
  nxcat: r_catg[k]
  p: cases_in_t
  size: n_cases*/
{
  int    k, i, j, ucat;
  int    nx;
  int    *r_ivector(int );
  for(k=0;k<nxcat;k++)//initial xu to all 0
    xu[k]=0;
  for(i=0;i<size;i++)
  {
    j=p[i];//serial number of observation
    nx = (int) x[j];//kth variable of jth observation
    xu[nx]=1;
  }
  ucat=0;
  for(i=0; i<nxcat; i++) ucat += (xu[i]>0);
  if(ucat<=1) return(1);
  return ucat;
}

/* Chin-Pei's modified version
Note that original odd and bodd are index by row first,
in this version they are index by column first */
void split_via(int t,double **xx,int *cases_in_t,int n_cases,
                       int *npl,int *nnd, NumericVector spvl, int k, int kk)
{
  int    i, j, k1, j1;
  int    *bin, *xu, ucat;
  void r_grayrep(int ,int *,int );
  double thresh;
  for(i=0; i<2*r_n_class; i++)
  odds[i]=0;
  if(r_catg[k]>1 && kk)
  {
    bin=r_ivector(r_catg[k]);
    xu=r_ivector(r_catg[k]);
    ucat=split_cases(xx[k], xu, r_catg[k], cases_in_t, n_cases);
    r_grayrep(kk,bin,ucat-1);
    bin[ucat-1]=0;
    k1=0;
    for(j=0; j<r_catg[k]; j++)
    if(xu[j] != 0)
      xu[j] = bin[k1++];
    thresh=0.;
    for(j=0; j<r_catg[k]; j++)
    {
      if(xu[j] != 0)
      {
        thresh = thresh*10.+j+1;
      }
    }
    for(i=0; i<n_cases; i++)
    {
      j=cases_in_t[i];
      k1=(int)xx[k][j];
      if(xu[k1])
      {
        j1 = 2*(r_n_class-1-r_y[j]);
        /* =2*y[j] if odds[0] is the count for y=0 (smallest) */
        odds[j1] ++; /* odds[0] is the count for y=n_class-1 */
        npl[j] = *nnd+1;
      }
      else
      {
        j1 = 2*(r_n_class-1-r_y[j]) + 1;
        /* =2*y[j]+1 if odds[0] is the count for y=0 (smallest) */
        odds[j1] ++; /* odds[0] is the count for y=n_class-1 */
        npl[j]= *nnd;
      }
    }
    free((char*)bin), free((char*)xu);
    spvl[t]=thresh;
  }
  else
  {
    for(i=0; i<n_cases; i++)
    {
      j=cases_in_t[i];
      if(xx[k][j] <= spvl[t])
      {
        npl[j]= *nnd;
        j1 = 2*(r_n_class-1-r_y[j]) +1;
        /* =2*y[j] if odds[0] is the count for y=0 (smallest) */
        odds[j1] ++; /* odds[0] is the count for y=n_class-1 */
      }
      else
      {
        npl[j]= *nnd+1;
        j1 = 2*(r_n_class-1-r_y[j]) ;
        /* =2*y[j]+1 if odds[0] is the count for y=0 (smallest) */
        odds[j1] ++; /* odds[0] is the count for y=n_class-1 */
      }
    }
  }
}

void find_dist(int t, int *npl, int *dist)
{
  int i, j;
  for(j=0; j<r_n_class; j++) dist[j]=0;
  for(i=0; i<r_n_obs; i++)
  {
    if(npl[i]==t)
    {
      j=r_y[i];
      dist[j]++;
    }
  }
}

void find_dist1(int t, int *npl, IntegerVector dist)
{
  int i, j;
  for(j=0; j<r_n_class; j++) dist[j]=0;
  for(i=0; i<r_n_obs; i++)
  {
    if(npl[i]==t)
    {
      j=r_y[i];
      dist[j]++;
    }
  }
}

int find_min_y(int *x)
{
  int  i, minv;

  minv=x[0];
  for(i=1; i<r_n_obs; i++)
  if(x[i] < minv) minv=x[i];
  return minv;
}

double find_chisq(int *bodds, int nn)
{
  double o[20], r[2], c[10], d=0., d1, e[20];
  int i, j, tmp;
  /*calculate Chi-square independence test statistic: 2x(n_class) table*/

  tmp=0;
  for(i=0;i<2*r_n_class;i++)
  {
    o[i] = bodds[i];
    if (bodds[i]<5)
    tmp++;
  }

  if(tmp>0)
  {
    for(i=0;i<2*r_n_class;i++)
      o[i] += 0.5;
    nn += r_n_class;
  }
  for(i=0;i<r_n_class;i++)
    c[i]=0;
  for(j=0;j<2;j++)
  {
    r[j]=0;
    for(i=0;i<r_n_class;i++)
    {
      r[j] += o[2*i+j];
      c[i] += o[2*i+j];
    }
  }
  for(i=0;i<r_n_class;i++)
  {
    e[2*i] = c[i]*r[0]/nn;
    e[2*i+1] = c[i]*r[1]/nn;
  }
  for(i=0; i<2*r_n_class; i++)
  {
    d1 = o[i]-e[i];
    d1 *= d1/e[i];
    d += d1;
  }
  return d;
}

void prune_tree(IntegerVector dt,int *npl, int *chi, int *nnd)
{
  int i, j, k=0, t1, t2;
  r_beforeprune = *nnd;
  for(i=*nnd-1; i>=0; i--)
  {
    if(dt[i] && chi[i]==0)
    {
      t1 = dt[i];
      t2 = dt[i]+1;
      if(dt[t1]==0 && dt[t2]==0)
      {
        dt[i]=0;
        k++;
        for(j=0; j<r_n_obs; j++)
          if(npl[j]==t1 || npl[j]==t2) npl[j]=i;
      }
    }
  }
  *nnd -= 2*k;
}

double  impurity(int *cntr,int nr,int *cntl,int nl,int size,int srule)
{
  double imp;
  int   j;
  if(srule==1) imp = 0.;
  else imp = 1.;
  for(j=0; j<r_n_class; j++)
    if(srule==1)
    {
      if(cntr[j]>0 && nr > 0) imp -= cntr[j]*log(cntr[j]*1./nr)/size;
      if(cntl[j]>0 && nl > 0) imp -= cntl[j]*log(cntl[j]*1./nl)/size;
    }
    else
    {
      if(nr > 0.) imp -= 1.*cntr[j]*cntr[j]/(nr*size);
      if(nl < size) imp -= 1.*cntl[j]*cntl[j]/(nl*size);
    }
  return imp;
}

void check_imp(int *bodds)
{
  int *tmp1, *tmp2;
  int i;
  int r1=0, r2=0;
  double  impurity(int *,int,int *,int ,int ,int);
  tmp1=r_ivector(10), tmp2=r_ivector(10);

  for(i=0; i<r_n_class; i++)
  {
    tmp1[i] = bodds[2*i];
    tmp2[i] = bodds[2*i+1];
    r1 += tmp1[i];
    r2 += tmp2[i];
  }
  Rprintf("The impurity of the split=%f",
               impurity(tmp1, r1, tmp2, r2, r1+r2, 1));
  free((char*)tmp1), free((char*)tmp2);
}

void reorganize_tree(int *nnd, IntegerVector dt, IntegerVector pt, int *npl, int *dist,
                           IntegerVector spv, NumericVector spvl, IntegerMatrix fcounts, int *neworder)
{
  int i, j, k, yes, newsize=0, *oldorder;
  if(r_beforeprune < *nnd)
  {
    Rprintf("Error: a pruned tree cannot be larger than the original tree\n");
    //exit(-1);
  }
  for(i=0; i<r_beforeprune; i++)
    newsize += (dt[i] > 0);
  newsize = 2*newsize + 1;
  oldorder = r_ivector(r_beforeprune);
  for(k=0, i=0; i<r_beforeprune; i++)
  {
    if(dt[i])
    oldorder[i]=k, neworder[k++]=i;
    /*node i in the pre-pruned tree is node k in the post-pruned one*/
    else
    {
      find_dist(i, npl, dist);
      yes=0;
      for(j=0; (j<r_n_class && yes==0); j++) yes=(dist[j]>0);
      if(yes)
        oldorder[i]=k, neworder[k++]=i;
    }
  }
  if(k != newsize)
  {
    Rprintf("Cannot match the tree sizes before and after pruning\n");
    Rprintf("Press Enter to end this program!");
    getchar();
    //exit(-1);
  }
  else
  {
    *nnd = newsize;
    for(i=0; i< newsize; i++)
    {
      j = neworder[i];
      spv[i] = spv[j];
      spvl[i] = spvl[j];
      k = dt[j];
      dt[i] = oldorder[k];
      k = pt[j];
      pt[i] = oldorder[k];
    }
    for(i=0; i< r_n_obs; i++) npl[i]=oldorder[npl[i]];
    for(i=newsize-1; i>=0; i--)
    {
      if(dt[i]==0){
        IntegerVector Temp(r_n_class);
        for (int j = 0; j < r_n_class; j++) Temp[j] = fcounts(i,j);
        find_dist1(i, npl, Temp);
        for (int j = 0; j < r_n_class; j++) fcounts(i,j) = Temp[j];
      }
      else
      {
        for(j=0; j<r_n_class; j++)
        fcounts(i,j)=fcounts(dt[i],j)+fcounts(dt[i]+1,j);
      }
    }
  }
  free((char*)oldorder);
}

int *find_ori_catg(double * covar,int *k)
{
  void sort2(int ,double *, int *);
  int   *distinct_x(double *, int, int *);
  double *ftmp;
  int   *itmp, *distx, i;

  itmp = r_ivector(r_n_obs+1);
  ftmp = r_dvector(r_n_obs+1);
  for(i=0; i<r_n_obs; i++)
  {
    itmp[i+1] = 1;
    ftmp[i+1] = covar[i];
  }
  sort2(r_n_obs, ftmp, itmp);
  for(i=0; i<r_n_obs; i++)
  {
    ftmp[i]=ftmp[i+1];
    itmp[i]=itmp[i+1];
  }
  itmp=distinct_x(ftmp, r_n_obs, k);
  distx = r_ivector(*k);
  for(i=0; i<*k; i++)
  {
    distx[i] = (int) ftmp[itmp[i]-1];
  }
  free((char*)itmp);
  free((char*)ftmp);
  return distx;
}

void r_nrerror(const char *error_text)
/*numerical recipes standard error handler*/
{
  //fRprintf(stderr, "Numerical Recipes run-time error ...\n");
  //fRprintf(stderr, "%s\n", error_text);
  //fRprintf(stderr, "...now exiting to system...\n");
  //fRprintf(stderr, "Press Enter to end this program!");
  getchar();
  //exit(1);
}

void free_matrix(double **m,int nr)
{
  int i;
  for(i=nr-1;i>=0;i--) free((char*)m[i]);
  free((char*)m);
}

void free_dmatrix(double **m, int nr)
{
  int i;
  for(i=nr-1;i>=0;i--) free((char*)m[i]);
  free((char*)m);
}

double **dmatrix2(int nr, int *nc)
{
  int i;
  double **m;
  /*unsigned malloc();*/
  m=(double **) malloc((unsigned) nr*sizeof(double*));
  if(!m) r_nrerror("error in matrix");
  for(i=0;i<nr;i++)
  {
    m[i]=(double *) malloc((unsigned) (nc[i]+1)*sizeof(double));
    if(!m[i]) r_nrerror("error in dmatrix");
  }
  return m;
}

int **imatrix2(int nr,int *nc)
{
  int    i, **m;
  /*unsigned malloc();*/
  m=(int **) malloc((unsigned) nr*sizeof(int*));
  if(!m) r_nrerror("error in imatrix2");
  for(i=0;i<nr;i++)
  {
    m[i]=(int *) malloc((unsigned) (nc[i]+1)*sizeof(int));
    if(!m[i]) r_nrerror("error in imatrix2");
  }
  return m;
}

double *r_dvector(int nh)
/*allocate a double vector with range [nl, nh]*/
{
  double  *v;
  /*unsigned malloc();
  v=(double *)malloc((unsigned)nh*sizeof(double));
  */

  v=(double  *)calloc(nh, sizeof(double));
  /*
  Rprintf("\nr_dvector allocate nh=%d of double=%d of  %d times\n",nh,sizeof(double),++dvct);
  */
  if(!v)
  {
    r_nrerror("allocations failure in r_dvector");
  }
  return v;
}

double *ddvector(int nh)
/*allocate a double vector with range [nl, nh]*/
{
  double *v;
  v=(double *)malloc((unsigned)nh*sizeof(double));
  if(!v) r_nrerror("allocation failure in ddvector()");
  return v;
}

int *r_ivector(int nh)
/*allocate a int vector with range [nl, nh]*/
{
  int sz;
  int  *v;
  sz = sizeof(int);

  v=(int  *)calloc(nh, sz);
  if(!v)
  {
    r_nrerror("allocation failure in r_ivector()");
  }
  return v;
}

int co_ancestor(int i,IntegerVector pt)
{
  int j,k;
  j=pt[i];
  k=pt[i+1];

  while(j!=k)
  {
    j=pt[j];
    k=pt[k];
  }
  return(j);
}

int exist_in_randomlist(int i_temp, int i, int *index_random) // decide whether i_temp in the first i-th of index_random
{
  int j, flag = 0;
  for (j=0; j<i; j++)
  {
    if (i_temp == index_random[j])
      flag = 1;
  }
  return(flag);
}

int *icopy_point_array(int *oriarray, int len)
//copy a int vector with length number
{
  int i;
  int  *v;

  v=(int  *)calloc(len, sizeof(int));
  if(!v)
  {
    r_nrerror("allocation failure in icopy_point_array()");
  }
  for (i=0; i<len; i++)
  {
    v[i] = oriarray[i];
  } // for (i=0; i<len; i++)
  return v;
}

double *dcopy_point_array(double *oriarray, int len)
//copy a double vector with length number
{
  int i;
  double  *v;

  v=(double  *)calloc(len, sizeof(double));
  if(!v)
  {
    r_nrerror("allocation failure in dcopy_point_array()");
  }
  for (i=0; i<len; i++)
  {
    v[i] = oriarray[i];
  } // for (i=0; i<len; i++)
  return v;
}



int sum(int *data, int lenth)
{
  int i, summ = 0;
  for (i=0; i<lenth; i++)
  {
    summ += data[i];
  }
  return(summ);
}

void sort_array(int *data, int *index, int *sort_data, int lenth)
{
  int i, j, temp;
  // initial
  for(i=0; i<lenth; i++)
  {
    index[i] = i;
    sort_data[i] = data[i];
  } // for(i=0; i<lenth; i++)

  for (i=0; i<lenth; i++)
  {
    for (j=i+1; j<lenth; j++)
    {
      if (sort_data[j] > sort_data[i])
      {
        temp = sort_data[i];
        sort_data[i] = sort_data[j];
        sort_data[j] = temp;
        temp = index[i];
        index[i] = index[j];
        index[j] = temp;
      } // if (data[j] > data[i])
    } // for (j=i+1; j<lenth; j++)
  } // for (i=0; i<lenth; i++)

  Rprintf("data="); for(i=0;i<10;i++) Rprintf(" %d", data[i]); Rprintf("\n");
  Rprintf("sort_data="); for(i=0;i<10;i++) Rprintf(" %d", sort_data[i]); Rprintf("\n");
  Rprintf("index="); for(i=0;i<10;i++) Rprintf(" %d", index[i]+1); Rprintf("\n");
}

void free_linklist(struct linklist *list)
{
  struct linklist *temp;
  while (list!=NULL)
  {
    temp = list->next;
    free((char*)list);
    list = temp;
  }
}



void readdata(DataFrame dataset,int ncol, IntegerVector colname)
{
  int i, k, flag = 0, *yclass;
  unsigned int j;
  /*std::vector<double> coltemp;*/
  /*Rcpp::DataFrame coltemp;*/
  NumericVector coltemp;
  coltemp = dataset[0];
  r_n_obs = coltemp.size();
  r_y=r_ivector(r_n_obs);
  if (!r_y) r_nrerror("error in vector y");
  r_X = (double  **) calloc(2*ncol, sizeof(double *));
  if(!r_X) r_nrerror("error in matrix X");
  for(i=0;i<ncol;i++)
  {
    r_X[i]=(double  *) calloc(r_n_obs, sizeof(double));
    if (!r_X[i]) r_nrerror("error in matrix X");
  }

  for(j=0; j<coltemp.size(); j++)
  {
    r_y[j] = coltemp[j];
  }
  for(i=1; i<dataset.length(); i++)
  {
    coltemp = dataset[i];
    for(j=0; j<coltemp.size(); j++)
    {
        r_X[i-1][j] = coltemp[j];
    }
  }
  yclass = (int *)malloc(sizeof(int)*r_n_obs);
  for(i=0;i<r_n_obs;i++)
  {
    if(!i)
    {
      yclass[0] = r_y[0];
      r_n_class+=1;
    }
    else
    {
      flag = 0;
      for(k=0;k<r_n_class ;k++)
        if(r_y[i] == yclass[k])
        {
          flag = 1;
          break;
        }
      if(flag==0)
      {
        yclass[r_n_class] = r_y[i];
        r_n_class+=1;
      }
    }
  }
  free(yclass);
}

void input_rtree(DataFrame dataset, IntegerVector colname)
{
  int    i, j, ncol=0;
  int    miny;
  int    find_min_y(int *);
  struct linklist *icatgc;

  r_n_class = 0;
  r_n_meas=0;
  icatgc=NULL;
  for(i=0;i<colname.size();i++)
  {
    if(colname[i]==1)  //orderial
      icatgc=addmiss(icatgc,r_n_meas,1);
    else if(colname[i]>1)
      icatgc=addmiss(icatgc,r_n_meas,colname[i]);
    else if(colname[i]==0)
      icatgc=addmiss(icatgc,r_n_meas,0);
    else if(colname[i]==-1)
      icatgc=addmiss(icatgc,r_n_meas,-1);
    else
    {
      Rprintf("ERROR: The colname is invalidate\n");
      //exit(-1);
    }
    r_n_meas++;
  }
  r_n_meas--;
  if(r_n_meas<1)
  {
    Rprintf("ERROR: No covariate!\n");
    //exit(-1);
  }
  else
  {
    r_catg=r_ivector(2*r_n_meas+1);
    if(icatgc==NULL)
      Rprintf("ABNORMAL EXIT\n");
    while(icatgc!=NULL)
    {
      r_catg[icatgc->col] = icatgc->row;
      icatgc=icatgc->next;
    }
  }
  for(i=0;i<=r_n_meas;i++)
    ncol += (r_catg[i] > 0);
  for(i=0, j=0;i<=r_n_meas; i++)
  {
    if(r_catg[i] > 0)
      r_catg[j++]=r_catg[i];
  }
  i=0;
  readdata(dataset,ncol,colname);
  if(r_n_class==1) r_n_class=2;
  r_n_meas=ncol;
  miny=find_min_y(r_y);
  for(i=0; i<r_n_obs; i++) r_y[i] -= miny;
  for(i=0; i<ncol;i++)
  {
    if(colname[i]>1)
    {
      for(j=0;j<r_n_obs;j++)
        r_X[i-1][j]=r_X[i-1][j]-1;
    }
  }
  free_linklist(icatgc);
}
