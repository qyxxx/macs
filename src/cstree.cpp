/*the part is edited to compare several method for some real data*/
/*updated on December 10, 2009 from ver7. This was to correct the error in
printing the categorical variable. It involved two functions:
split_via_best() and print_int_string()
*/

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include "head.h"
#include "R_ext/Print.h"
using namespace Rcpp;

#define MAXDEP 10
#define MINNODE 15
#define MAXIMP 100000.0
#define rmiss "NA" /* missing value */
#define rmiss1 "." /* another possible missing value */
#define true 1
#define false 0

int n_obs; /*sample size or number of observations or cases*/
int n_meas; /*number of measurements taken for each cases*/
int n_out, *catg, *n_class, srule, *maxcatg;
int MAXNODES=40;
int **ori_catg;
int *var_to_var, *order;
int ***posi1, **casesInt, **rank1;
float **sumy1, **sums, *lrank;
float *cntl1, *cntr1, **cntl2, **cntr2;
float **X, *tmpx;
float **y, **ty, *unceny, *medy;
float ***cost, **prior, tradof, **yprop, *ssy;
float **cp, **tmp12;
int *disx, *c1, *c2, *xdist, *tl, *tr, *indc;
int *xu, *bin, *xcat, **counts, **ynum, **ycnt, **ycat;
int **spv1, **pt1, **indi, **indj;
float **spvl1;
int adp, missr, *icatg, beforeprune=0;
char outfile[100], prefix_file[50], figname[50];
IntegerVector death_catg(100);// number of death in each node

// [[Rcpp::export]]
List cstree(DataFrame dataset, IntegerVector colname, int srule)
{
  void input(DataFrame, IntegerVector);
  void define_cost();
  float *r,*dvector(int );
  NumericVector spvl(MAXNODES);
  IntegerVector spv(MAXNODES);
  IntegerVector dt(MAXNODES);
  IntegerVector pt(MAXNODES);
  IntegerVector ncases(MAXNODES);
  NumericVector median(MAXNODES);

  int *npl, **clm, nnd;
  int *ivector(int ), **imatrix(int , int );
  float **meany, **matrix(int , int );
  int similarity(float **, float **, int , float *, IntegerVector , NumericVector , IntegerVector , IntegerVector , int *, int **, IntegerVector , float **, NumericVector , int );
  input(dataset,colname);
  define_cost();
  r=dvector(MAXNODES);
  clm=imatrix(MAXNODES, n_out);
  meany=matrix(MAXNODES, n_out);
  lrank=dvector(MAXNODES);
  npl=ivector(n_obs);
  medy=dvector(n_obs);
  order=ivector(n_obs);
  spv1=imatrix(5, MAXNODES);
  spvl1=matrix(5, MAXNODES);
  pt1=imatrix(5, MAXNODES);

  nnd=similarity(X, y, n_obs, r, spv, spvl, pt, dt,npl, clm, ncases, meany, median, srule);

  return List::create(
    Named("nnd") = nnd,
    Named("dt") = dt,
    Named("pt") = pt,
    Named("spv") = spv,
    Named("spvl") = spvl,
    Named("ncases") = ncases,
    Named("death_catg") = death_catg,
    Named("median") = median,
    Named("colname") = colname
  );
}

int similarity(float **xx, float **yy, int n, float *r, IntegerVector spv, NumericVector spvl, IntegerVector pt, IntegerVector dt, int *npl, int **clm, IntegerVector ncases, float **meany, NumericVector median, int srule)
{
  int nnd = 0;
  int keyin2(int, char *);
  void grow_large_tree1(float **, float **, int , int *, float *, IntegerVector , NumericVector , IntegerVector , IntegerVector ,
                        int *, IntegerVector , float **, NumericVector );
  void grow_large_tree2(float **, float **, int , int *, float *, IntegerVector , NumericVector ,
                        IntegerVector , IntegerVector , int *, IntegerVector , NumericVector );
  void grow_large_tree3(float **, float **, int , int *, float *, IntegerVector , NumericVector , IntegerVector ,
                        IntegerVector , int *, IntegerVector , NumericVector );
  void grow_large_tree(float **, float **, int , int *, float *, IntegerVector , NumericVector , IntegerVector , IntegerVector ,
                       int *, int **, IntegerVector , float **, NumericVector );
  void clean_y(), specify_prior();
  void prune_tree(IntegerVector , int *, float *, int *);
  void nrerror(const char *);

  if (srule==1)
  {
    grow_large_tree1(xx, yy, n, &nnd, r, spv, spvl, pt, dt, npl, ncases, meany, median);
    prune_tree(dt, npl, lrank, &nnd);
  }
  else if (srule==2)
  {
    grow_large_tree2(xx, yy, n, &nnd, r, spv, spvl, pt, dt, npl, ncases, median);
    prune_tree(dt, npl, lrank, &nnd);
  }
  else if (srule==3)
  {
    grow_large_tree3(xx, yy, n, &nnd, r, spv, spvl, pt, dt,npl, ncases, median);
    prune_tree(dt, npl, lrank, &nnd);
  }
  else if (srule==4)
  {
    tradof=0.5, adp=1;
    clean_y(), specify_prior();
    grow_large_tree(xx, yy, n, &nnd, r, spv, spvl, pt, dt, npl, clm, ncases, meany, median);
    prune_tree(dt, npl, lrank, &nnd);
  }
  else if (srule==5)
  {
    tradof=1.0, adp=0;
    clean_y(), specify_prior();
    grow_large_tree(xx, yy, n, &nnd, r, spv, spvl, pt, dt,
                    npl, clm, ncases, meany, median);
    prune_tree(dt, npl, lrank, &nnd);
  }
  else {
    Rprintf("wrong input of split rule!\n");
    nrerror("exit from split");
  }
  return nnd;
}

void find_min_max(float *x, int n, float *minv, float *maxv)
{
  int i;

  *minv=*maxv=x[0];
  for(i=1; i<n; i++)
  {
    if(x[i] > *maxv) *maxv=x[i];
    if(x[i] < *minv) *minv=x[i];
  }
}

void define_cost()
{
  int j, k, kk;
  float **matrix(int, int);
  /* unsigned malloc(); */

  cost=(float***)malloc((unsigned)n_out*sizeof(float**));
  for(kk=0; kk<n_out; kk++) if(n_class[kk] > 1)
    cost[kk]=matrix(n_class[kk], n_class[kk]);
  for(kk=0; kk<n_out; kk++) if(n_class[kk] > 1)
    for(j=0; j<n_class[kk]; j++)
      for(k=0; k<n_class[kk]; k++)
        cost[kk][j][k]=1.;
  for(kk=0; kk<n_out; kk++) if(n_class[kk] > 1)
    for(j=0; j<n_class[kk]; j++) cost[kk][j][j]=0.;
  cost[1][1][0]=2.;
}

void clean_y()
{
  void find_min_max(float *, int , float *, float *);
  void find_yvar(float *, int , float *, float *);
  float minv, maxv, tmp, tmp1;
  int kk;

  for(kk=0; kk<n_out; kk++)
  {
    find_min_max(y[kk], n_obs, &minv, &maxv);
    find_yvar(y[kk], n_obs, &tmp, &tmp1);
    ssy[kk]=tmp1-tmp*tmp/n_obs;
  }
  yprop[1][1]=tmp;
  yprop[1][0]=n_obs-yprop[1][1];
}

void specify_prior()
{
  int i, kk;

  for(kk=0; kk<n_out; kk++)
    if(n_class[kk]>1)
      for(i=0; i<n_class[kk]; i++)
        prior[kk][i] = 1.*yprop[kk][i]/n_obs;
}

void find_yvar(float *xx1, int n, float *sumy, float *sum2y)
{
  int i;

  *sumy = *sum2y =0.;
  for(i=0; i<n; i++)
    *sumy += xx1[i], *sum2y += xx1[i]*xx1[i];
}

/*************grow.c*************/

float ord_best_split(float *x, float **yy, int *rank, int *cases_in_t, int n, float *r1, float *r2, float *r, float *thresh, float *sumy)
{
  float min_mr, imp;
  float impurity(float **, float **, int , int );
  int i, j, jobs;
  int k, kk;
  void distinct_x(float *, int , int *, int *);
  void prob(int , int , int , float *, float *), find_costs(float **, float **, float *, float *, float *, int , int );

  for(i=0; i<n; i++)
  {
    j=cases_in_t[i];
    k=rank[j];
    if(k < 0 || k >= n) Rprintf("rank %d out of bound\n", k);
    for(kk=0; kk<n_out; kk++)
      ty[kk][k]=yy[kk][j];
    tmpx[k]=x[j];
  }
  distinct_x(tmpx, n, &k, disx);
  min_mr=1.e10;
  jobs=0;
  if(k>1)
  {
    for(kk=0; kk<n_out; kk++)
    {
      for(j=0; j<= n_class[kk]; j++)
        cntl2[kk][j]=0.;
      if(n_class[kk] > 1)
        prob(0, n, kk, ty[kk], cntr2[kk]);
      else
        cntr2[kk][0]=sumy[kk], cntr2[kk][1]=sumy[n_out+kk];
    }
    cntl1[1]=0., cntr1[1]=sumy[1];   //C-Y
    for(i=1; i<k; i++)
    {
      for(kk=0; kk<n_out; kk++)
      {
        if(n_class[kk] > 1)
        {
          prob(disx[i-1], disx[i], kk, ty[kk], tmp12[kk]);
          for(j=0; j<n_class[kk]; j++)
            cntl2[kk][j] += tmp12[kk][j], cntr2[kk][j] -= tmp12[kk][j];
          for(j=disx[i-1]; j<disx[i]; j++)
            cntl1[1] += 1.-ty[1][j], cntr1[1] -= 1.- ty[1][j];
        }
        else
        {
          for(j=disx[i-1]; j<disx[i]; j++)
          {
            cntl2[kk][0] += ty[kk][j], cntr2[kk][0] -= ty[kk][j];
            cntl2[kk][1] += ty[kk][j]*ty[kk][j];
            cntr2[kk][1] -= ty[kk][j]*ty[kk][j];
          }
        }
      }
      imp = impurity(cntl2, cntr2, disx[i], n-disx[i]);

      if(imp < min_mr)
      {
        find_costs(cntl2, cntr2, r, r1, r2, disx[i], n-disx[i]);
        min_mr = imp;
        jobs = disx[i]-1;
        for(kk=0; kk<n_out; kk++)
        {
          if(n_class[kk]==1)
          {
            sums[kk][0]=cntl2[kk][0], sums[kk][1]=cntl2[kk][1];
            sums[kk][2]=cntr2[kk][0], sums[kk][3]=cntr2[kk][1];
          }
          else
            sums[kk][1]=cntl1[1], sums[kk][3]=cntr1[1];
        }
      }
    }
  }
  *thresh=tmpx[jobs];
  return min_mr;
}

float catg_best_split(float *x, float **yy, int nxcat, int *p, int size, float *r1, float *r2, float *r, float *lcutc, float *sumy)
{
  float imp, timp;
  int i, j, k, kk, l, m, ucat=0, ulim, isave;
  int nc, nx, ny, nl, nr;
  void find_costs(float **, float **, float *, float *, float *, int , int );
  float impurity(float **cl, float **cr, int nl, int nr);
  //   FILE *fpp;
  //fpp=fopen("debugstree.txt","a");

  timp=1.e10;
  isave=1;
  for(k=0;k<nxcat;k++)
  {
    xu[k]=xcat[k]=0; xdist[k]=0;
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk] > 1) /* censor */
      {
        ynum[kk][k]=0;
        for(i=0;i<n_class[kk];i++)
          counts[kk][k*n_class[kk]+i]=0;
        ty[kk][k]=0.;
      }
      else /* response */
ty[kk][k]=ty[kk][nxcat+k]=0.;
    }
  }
  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
      for(i=0;i<n_class[kk];i++)
        cntl2[kk][i]=cntr2[kk][i]=0.;
    else
    {
      cntl2[kk][0]=cntl2[kk][1]=0.;
      cntr2[kk][0]=sumy[kk], cntr2[kk][1]=sumy[n_out+kk];
    }
  }
  for(i=0;i<size;i++)
  {
    j=p[i];
    if(j >= n_obs || j < 0)
      Rprintf("error obs ind\n");
    nx = (int) x[j];
    xu[nx]=1;
    xdist[nx]++;
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk] > 1) /* censor */
      {
        ny = (int)yy[kk][j]; /* nxcat*n_class cell */
if(ny >= n_class[kk] || ny < 0)
  Rprintf("error in y %d %d\n", ny, n_class[kk]);
counts[kk][nx*n_class[kk] + ny]++;
cntr2[kk][ny]++;
      }
      else
      {
        ty[kk][nx] += yy[kk][j];
        ty[kk][nxcat+nx]+=yy[kk][j]*yy[kk][j];
      }
    }
  }
  for (i=0, ucat=0; i<nxcat; i++) /* not <= */
if(xu[i]==1) xcat[ucat++]=i; /* category index */

if (ucat<=1) return(MAXIMP);
for(k=0;k<ucat;k++)
{
  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
    {
      for(m=0;m<n_class[kk];m++)
      {
        nc=counts[kk][xcat[k]*n_class[kk] + m];
        /*nc is the counts in cell x:xcat[k] and y:m*/
        if(nc)
        {
          ycnt[kk][k*n_class[kk]+ynum[kk][k]]=nc;
          ycat[kk][k*n_class[kk]+ynum[kk][k]]=m;
          ynum[kk][k]++;
        }
      }
    }
  }
}

/*ynum counts non-empty cells for each level of x*/
/*ycat matches to the original levels of y*/
/*ycnt excludes empty cells from counts*/
ulim = (((unsigned)1)<<(ucat-1)) - 1;
nr=size, nl=0;
for(i=1; i<=ulim; i++) /* all possible splits */
{
  l=i; k=0;
  while(!(l&1))
  {
    l>>=1; k++;
  }
  if(!(l&2))
  { /* bit goes 0 -> 1 */
nl += xdist[xcat[k]], nr -= xdist[xcat[k]];
  }
  else
  {
    nl -= xdist[xcat[k]], nr += xdist[xcat[k]];
  }
  //	Rprintf("i=%d,l=%d,k=%d,nr=%d,nl=%d",i,l,k,nr,nl);
  cntl1[1]=0., cntr1[1]=0.;
  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
    {
      for(j=0; j<ynum[kk][k]; j++)
      {
        nc = ycnt[kk][k*n_class[kk]+j];
        m = ycat[kk][k*n_class[kk]+j];

        if(!(l&2))
        {
          /* bit goes 0 -> 1 */
          cntl2[kk][m]+=nc;
          cntr2[kk][m]-=nc;
          cntl1[1]+=ty[1][xcat[k]], cntr1[1] -= ty[1][xcat[k]];
        }
        else
        {
          /* bit goes 1 -> 0 */
          cntl2[kk][m]-=nc;
          cntr2[kk][m]+=nc;
          cntl1[1]-=ty[1][xcat[k]], cntr1[1] += ty[1][xcat[k]];
        }

      }
    }
    else
    {
      if(!(l&2))
      {
        /* bit goes 0 -> 1 */
        cntl2[kk][0]+= ty[kk][xcat[k]], cntr2[kk][0] -= ty[kk][xcat[k]];
        cntl2[kk][1] += ty[kk][nxcat+xcat[k]];
        cntr2[kk][1] -= ty[kk][nxcat+xcat[k]];
      }
      else
      {
        cntl2[kk][0] -= ty[kk][xcat[k]], cntr2[kk][0] += ty[kk][xcat[k]];
        cntl2[kk][1] -= ty[kk][nxcat+xcat[k]];
        cntr2[kk][1] += ty[kk][nxcat+xcat[k]];
      }
    }

  }
  imp = impurity(cntl2, cntr2, nl, nr);

  //	fprintf(fpp,"i=%d, nl=%d, nr=%d, imp=%f \n",i,nl,nr,imp);
  if(imp < timp)
  {
    /* save minimp and best split */

    find_costs(cntl2, cntr2, r, r1, r2, nl, nr);
    isave = i;
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk]==1)
      {
        sums[kk][0]=cntl2[kk][0], sums[kk][1]=cntl2[kk][1];
        sums[kk][2]=cntr2[kk][0], sums[kk][3]=cntr2[kk][1];
      }
      else
        sums[kk][1]=cntl1[1], sums[kk][3]=cntr1[1];
    }
    timp = imp;

  }
}
*lcutc = isave;
//   fprintf(fpp,"return index=%d, minimp=%f\n",isave, timp);
//   fclose(fpp);
return(timp);
}

float impurity(float **cl, float **cr, int nl, int nr)
{
  int kk, j;
  float imp, prtl, prtr, n1, d1;
  void nrerror(const char *);

  if(nl == 0 || nr == 0) nrerror("zero cases in one side");
  if(2*nl <= MINNODE || 2*nr <= MINNODE) return MAXIMP;
  imp=0.;
  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
    {
      prtl=prtr=0.;
      for(j=0; j<n_class[kk]; j++)
      {
        prtl += cl[kk][j]*yprop[kk][j];
        prtr += cr[kk][j]*yprop[kk][j];
      }
      for(j=0; j<n_class[kk]; j++)
      {
        if(cl[kk][j]*yprop[kk][j] > 0. && prtl > 0.)
          imp -= cl[kk][j]*yprop[kk][j]*
            log(cl[kk][j]*yprop[kk][j]/prtl)/(prtl+prtr);
        if(cr[kk][j]*yprop[kk][j] > 0. && prtr > 0.)
          imp -= cr[kk][j]*yprop[kk][j]*
            log(cr[kk][j]*yprop[kk][j]/prtr)/(prtl+prtr);
      }
    }
    else
    {
      d1 = cl[kk][1]+cr[kk][1] -
        (cl[kk][0]+cr[kk][0])*(cl[kk][0]+cr[kk][0])/(nl+nr);
      if(d1 > 0.)
      {
        n1 = cl[kk][1] - cl[kk][0]*cl[kk][0]/nl;
        if(n1 < 0.)
        {
          if(n1 < -1.e-6)  // printf("warning: n1=%f\n", n1);
            n1=0.;
        }
        if(adp) imp += tradof*n1/d1;
        else imp += tradof*n1/ssy[kk];
        n1 = cr[kk][1] - cr[kk][0]*cr[kk][0]/nr;
        if(n1 < 0.)
        {
          if(n1 < -1.e-6)  // printf("warning: n1=%f\n", n1);
            n1=0.;
        }
        if(adp) imp += tradof*n1/d1;
        else imp += tradof*n1/ssy[kk];
      }
    }
  }
  return imp;
}

void grow_large_tree(float **xx, float **yy, int n, int *nnd, float *r, IntegerVector spv, NumericVector spvl, IntegerVector pt, IntegerVector dt,
                     int *npl, int **clm, IntegerVector ncases, float **meany, NumericVector median)
{
  void initialize_tree(float** , int , int *, float *, IntegerVector , IntegerVector , NumericVector , int **);
  void split_a_node(int , float **, float **, int ***, int **, int *, float *, IntegerVector , NumericVector ,
                    IntegerVector , IntegerVector , int *, int **, int **, IntegerVector , float **, NumericVector );
  int i, kk, t;
  void map_distx(float **, int , int , int **, int **), find_yvar(float *, int , float *, float *);
  void adjust_prior(float **y, int n);

  for(kk=0; kk<n_out; kk++)
    if(n_class[kk]==1)
      find_yvar(yy[kk], n, &sumy1[0][kk], &sumy1[0][n_out+kk]);
    sumy1[0][1]=n-sumy1[0][1];
    adjust_prior(yy, n);
    initialize_tree(yy, n, npl, r, dt, spv, spvl, clm);
    *nnd=1;
    map_distx(xx, n_meas, n, posi1[0], rank1);
    for(i=0; i<n; i++) casesInt[0][i]=i;
    ncases[0]=n;
    for(t=0, *nnd=1; t < *nnd; t++)
    {
      split_a_node(t, xx, yy, posi1, rank1, nnd, r, spv, spvl, pt,
                   dt, npl, clm, casesInt, ncases, sumy1, median);
      for(kk=0; kk<n_out; kk++) if(n_class[kk]==1)
        meany[t][kk]=sumy1[t][kk]/ncases[t];
    }

}

void split_a_node(int t, float **xx, float **yy, int ***posi, int **rank, int *nnd, float *r, IntegerVector spv, NumericVector spvl,
                  IntegerVector pt, IntegerVector dt, int *npl, int **clm, int **cases_in_t, IntegerVector ncases, float **sumy, NumericVector median)
{
  float min_mr, *dvector(int ),  find_best_split(int , float **, float ** , int **, int , int *, float *, IntegerVector , NumericVector ,
                         int **, int *, float **);
  void split_via_best(int , float **, int **, int ***, int , int *, IntegerVector , NumericVector ,
                      IntegerVector , IntegerVector , int *, int **, IntegerVector ), sort2(int , float *, int *);
  int n, yes, i, j, kk, *ivector(int ), check_if_needed_split(float **, float **, int *, int , int , IntegerVector );
  int count_death=0;

  n=ncases[t];
  for(i=0; i<n; i++)
  {
    order[i]=i+1;
    j=cases_in_t[t][i];
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk]==1)
        medy[i+1]=yy[kk][j];
      else
        count_death+=yy[kk][j];
    }
  }

  death_catg[t]=n-count_death;
  sort2(n, medy, order);

  if(n%2==1)  median[t]=medy[n/2+1];
  else    median[t]=(medy[n/2]+medy[n/2+1])/2.0;

  yes = (*nnd < MAXNODES-1);
  if(yes)
  {
    yes=check_if_needed_split(xx, yy, cases_in_t[t], n, t, pt);
  }

  if(yes)
  {
    min_mr=find_best_split(t, xx, yy, rank, n, nnd, r, spv, spvl,
                           clm, cases_in_t[t], sumy);
    if(min_mr < MAXIMP)
    {
      split_via_best(t, xx, rank, posi, n, nnd, spv, spvl,
                     pt, dt, npl, cases_in_t, ncases);

    }
    else dt[t]=0;

  }
  else
    dt[t]=0;
  /*  free((char*)order);
  free((char*)medy);       */

}

int check_if_needed_split(float **xx, float **yy, int *cases, int n_cases, int t, IntegerVector pt)
{
  int yes1=0, yes2=0, i, j, i0, k;
  int depth, parentn;

  if(n_cases<=MINNODE) return 0;
  parentn=t, depth=0;
  while(parentn)
  {
    parentn=pt[parentn];
    depth++;
  }
  if(depth > MAXDEP) return 0;
  i0=cases[0];
  for(k=0; k<n_out; k++)
    for(i=1; i<n_cases; i++)
    {
      j=cases[i];
      if(yy[k][j] != yy[k][i0]) yes1=1;
    }
    if(yes1)
    {
      for(i=0; i<n_meas; i++)
        for(k=1; k<n_cases; k++)
        {
          j=cases[k];
          if(xx[i][j] != xx[i][i0]) yes2=1;
        }
    }
    return yes2;
}

void grayrep(int i, int *v, int n)
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

int split_cases(float *x, int *xu, int nxcat, int *p, int size)
{
  int k, i, j, ucat;
  int nx;

  for(k=0;k<nxcat;k++)
    xu[k]=0;
  for(i=0;i<size;i++)
  {
    j=p[i];
    nx = (int) x[j];
    xu[nx]=1;
  }
  ucat=0;
  for (i=0; i<nxcat; i++) ucat += (xu[i]>0);
  if (ucat<=1) return(1);
  return ucat;
}

void distinct_x(float *tmpx1, int n_cases, int *k, int *ddisx)

{
  int i, j;

  *k=1; ddisx[0]=ddisx[*k]=0, j=0;
  for(i=0; i<n_cases; i++)
  {
    if(tmpx1[i]==tmpx1[j])
      ddisx[*k] +=1;
    else
    {
      *k += 1;
      ddisx[*k] = ddisx[*k-1]+1;
      j=ddisx[*k-1];
    }
  }
}

void initialize_tree(float** y, int n, int *npl, float *r, IntegerVector dt, IntegerVector spv, NumericVector spvl, int **clm)
{
  int i;
  float misrate(int , int , float **, int *);

  for(i=0; i<n; i++) npl[i]=0;
  r[0]=misrate(0, n, y, clm[0]);
  for(i=0; i<MAXNODES; i++)
  {
    dt[i]=0;
    spv[i]=0, spvl[i]=0.;
  }
}

float misrate(int i1, int i2, float **tmpc, int *c0)
{
  int i, j, kk;
  float tmr, mr, min_mr;
  void prob(int , int , int , float *, float *);

  tmr=0.0;
  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
    {
      prob(i1, i2, kk, tmpc[kk], cp[kk]);
      c0[kk]=0;
      for(j=0, min_mr=0.; j<n_class[kk]; j++)
        min_mr += cp[kk][j]*cost[kk][j][0]*yprop[kk][j];
      for(i=1; i<n_class[kk]; i++)
      {
        for(j=0, mr=0.; j<n_class[kk]; j++)
          mr += cp[kk][j]*cost[kk][j][i]*yprop[kk][j];
        if(mr < min_mr)
        {
          min_mr=mr;
          c0[kk]=i;
        }
      }
      tmr += min_mr;
    }
    else tmr += 1.;
  }
  return tmr;
}

void prob(int i1, int i2, int kk, float *tmpc, float *prob1)

{
  int i, j;

  for(i=0; i<n_class[kk]; i++) prob1[i]=0.;
  for(j=i1; j<i2; j++) prob1[(int)tmpc[j]] += 1.;
}

void map_distx(float **x1, int nv, int n, int **posi, int **rank)
{
  int i, j;
  void sort2(int , float *, int *);

  for(j=0; j<nv; j++)
  {
    for(i=1; i<=n; i++)
    {
      posi[j][i]=i-1;
      tmpx[i]=x1[j][i-1];
    }
    sort2(n, tmpx, posi[j]);
    for(i=0; i<n; i++)
    {
      tmpx[i]=tmpx[i+1];
      posi[j][i]=posi[j][i+1];
    }
    for(i=0; i<n; i++) rank[j][posi[j][i]]=i;
  }
}

float find_best_split(int t, float **xx, float **yy, int **rank, int n, int *nnd, float *r, IntegerVector spv, NumericVector spvl,
                      int **clm, int *cases_in_t, float **sumy)
{
  int i, kk;
  float min_mr, r1, r2;
  float ord_best_split(float *, float **, int *, int *, int , float *, float *, float *, float *, float *), thresh, it;
  float catg_best_split(float *, float **, int , int *, int , float *, float *, float *, float *, float *);
  void find_min_max(float *, int , float *, float *);
  float log_rank_catg(float *, float **, int, int *, int,float), log_rank_ord(float *, float **, int *, int *, int, float);

  min_mr=MAXIMP;
  for(i=0; i<n_meas; i++)
  {
    if(catg[i])
    {
      if(catg[i] > 1)
      {
        //		  find_min_max(xx[i], n_obs, &minv, &maxv);
        //		  catg[i]=(int)maxv+1;  //Comment by C-Y why we need this here???
        //		  maxcat=(int)maxv+1; // won't change correct catg[i]
        it=catg_best_split(xx[i], yy, maxcatg[i], cases_in_t, n,
                           &r1, &r2, r, &thresh, sumy[t]);
      }
      else
      {

        it=ord_best_split(xx[i], yy, rank[i], cases_in_t, n,
                          &r1, &r2, r, &thresh, sumy[t]);
      }
      if(it < min_mr)
      {
        min_mr=it;
        spv[t]=i;
        spvl[t]=thresh;
        r[*nnd]=r1, r[*nnd+1]=r2;
        for(kk=0; kk<n_out; kk++)
        {
          if(n_class[kk] > 1)
          {
            clm[*nnd][kk]=c1[kk], clm[*nnd+1][kk]=c2[kk];
            sumy[*nnd][kk]=sums[kk][1], sumy[*nnd+1][kk]=sums[kk][3];
          }
          else
          {
            sumy[*nnd][kk]=sums[kk][0], sumy[*nnd][n_out+kk]=sums[kk][1];
            sumy[*nnd+1][kk]=sums[kk][2], sumy[*nnd+1][n_out+kk]=sums[kk][3];
          }
        }
      }
    }
  }
  if(catg[spv[t]]>1)
  { lrank[t]=log_rank_catg(xx[spv[t]], yy, maxcatg[spv[t]], cases_in_t, n,
                           spvl[t]);
  }
  else
  {
    lrank[t]=log_rank_ord(xx[spv[t]], yy, rank[spv[t]], cases_in_t, n, spvl[t]);
  }
  return min_mr;
}

void split_via_best(int t, float **xx, int **rank, int ***posi, int n, int *nnd, IntegerVector spv, NumericVector spvl,
                    IntegerVector pt, IntegerVector dt, int *npl, int **cases_in_t, IntegerVector ncases)
{
  int i, j, k, kk, keep_xu[50];
  int ucat, nl, nr;
  int split_cases(float *, int *, int , int *, int );
  void split_a_node(int , float **, float **, int ***, int **, int *, float *, IntegerVector , NumericVector ,
                    IntegerVector , IntegerVector , int *, int **, int **, IntegerVector , float **, NumericVector );
  void find_maps(int **, int **, int , int *, int , int **, int **);
  void grayrep(int , int *, int );
  float thresh;

  k=spv[t];
  kk=(int) spvl[t];
  nr=0, nl=0;
  if(catg[k]>1 && kk)
  {
    ucat=split_cases(xx[k], xu, maxcatg[k], cases_in_t[t], n);
    grayrep(kk,bin,ucat-1);
    bin[ucat-1]=0;
    // printf("%d %d %d", kk, ucat, maxcatg[k]);
    //for(j=0; j<maxcatg[k]; j++)
    // {
    //   printf("bin[%d]=%d ",j,bin[j]);
    //}
    kk=0;
    for(j=0; j<maxcatg[k]; j++)
    {
      keep_xu[j]=0;
      //	  printf("bin[%d]=%d ",j,bin[j]);
      if(xu[j] != 0)
      {
        keep_xu[j]=1;
        xu[j] = bin[kk++];
        //	      printf(" j=%d, xu[j]=%d ",j,xu[j]);
      }
    }
    //	for(j=0; j<maxcatg[k]; j++)
    //		printf(" xu[%d]=%d %d",j,xu[j], keep_xu[j]);

    thresh=0.;
    for(j=0; j<maxcatg[k]; j++)
    {
      //	     if(xu[j] != 0)
      if(xu[j] == 0 && keep_xu[j] == 1)  // take compliment
        thresh = thresh*10.+j+1;
    }
    //	printf("thresh=%f\n",thresh);
    for(i=0; i<n; i++)
    {
      j=cases_in_t[t][i];
      kk=(int)xx[k][j];
      if(xu[kk])
      {
        npl[j] = *nnd;
        tl[nl]=j;
        nl++;
      }
      else
      {
        npl[j]= *nnd+1;
        tr[nr]=j;
        nr++;
      }
    }
    spvl[t]=thresh;
  }
  else
  {
    for(i=0; i<n; i++)
    {
      j=cases_in_t[t][i];
      if(xx[k][j] <= spvl[t])
      {
        npl[j]= *nnd;
        tl[nl]=j;
        nl++;
      }
      else
      {
        npl[j]= *nnd+1;
        tr[nr]=j;
        nr++;
      }
    }
  }
  for(i=0; i<nl; i++) cases_in_t[*nnd][i]=tl[i];
  for(i=0; i<nr; i++) cases_in_t[*nnd+1][i]=tr[i];
  ncases[*nnd]=nl, ncases[*nnd+1]=nr;
  find_maps(posi[t], rank, n, cases_in_t[*nnd], nl, posi[*nnd], rank);
  find_maps(posi[t], rank, n, cases_in_t[*nnd+1], nr, posi[*nnd+1], rank);
  pt[*nnd]=t, pt[*nnd+1]=t;
  dt[t]= *nnd;
  *nnd += 2;
}

void find_maps(int **posi, int **rank, int n, int *cases_in_t, int nl, int **posi0, int **rank0)
{
  int i, j, jj, k;
  void nrerror(const char *);

  for(i=0; i<n_meas; i++)
  {
    if(catg[i] == 1)
    {
      for(j=0; j<n; j++) indc[j]=0;
      for(j=0; j< nl; j++)
      {
        k=cases_in_t[j];
        if(k < 0 || k >= n_obs) nrerror("out of range in find_maps:1");
        k=rank[i][k];
        if(k < 0 || k >= n_obs) nrerror("out of range in find_maps:2");
        indc[k]=1;
      }
      for(j=jj=0; j<n; j++)
        if(indc[j]) posi0[i][jj++]=posi[i][j];
        if(jj != nl)
          Rprintf("%d %d inconsistent numbers\n", jj, nl);
        for(j=0; j< nl; j++)
        {
          k = posi0[i][j];
          if(k < 0 || k >= n_obs) nrerror("out of range in find_maps:3");
          rank0[i][k]=j;
        }
    }
  }
}

void sort2(int n, float *ra, int *rb)
{
  int l,j,ir,i, rrb;
  float rra;

  l=(n >> 1)+1;
  ir=n;
  for (;;)
  {
    if (l > 1)
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

void adjust_prior(float **y, int n)
{
  int kk, i;

  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
    {
      for(i=0; i<n_class[kk]; i++) yprop[kk][i]=0.;
      for(i=0; i<n; i++) yprop[kk][(int)y[kk][i]] += 1.;
      for(i=0; i<n_class[kk]; i++)
      {
        if(yprop[kk][i] > 0)
          yprop[kk][i] = prior[kk][i]/yprop[kk][i];
      }
    }
  }
}

void find_costs(float **cl, float **cr, float *r, float *r1, float *r2, int nl, int nr)
{
  int kk, j1, j;
  float ml, mr, rr1, rr2;
  void nrerror(const char *);

  *r1=*r2=0.;
  for(kk=0; kk<n_out; kk++)
  {
    if(n_class[kk] > 1)
    {
      rr1=r[0]+1., rr2=r[0]+1.;
      for(j=0; j< n_class[kk]; j++)
      {
        ml=mr=0.;
        for(j1=0; j1< n_class[kk]; j1++)
        {
          ml += cl[kk][j1]*cost[kk][j1][j]*yprop[kk][j1];
          mr += cr[kk][j1]*cost[kk][j1][j]*yprop[kk][j1];
        }
        if(ml < rr1)
          rr1=ml, c1[kk]=j;
        if(mr < rr2)
          rr2=mr, c2[kk]=j;
      }
      *r1 += rr1, *r2 += rr2;
    }
    else
    {
      ml = cl[kk][1] - cl[kk][0]*cl[kk][0]/nl;
      if(ml < 0.)
      {
        if(ml < -1.e-6) // Rprintf("warning: ml=%f\n", ml);
          ml=0.;
      }
      *r1 += tradof*ml/ssy[kk];
      mr = cr[kk][1] - cr[kk][0]*cr[kk][0]/nr;
      if(mr < 0.)
      {
        if(mr < -1.e-6) // Rprintf("warning: mr=%f\n", mr);
          mr=0.;
      }
      *r2 += tradof*mr/ssy[kk];
    }
  }
}

/*****************surv.c***********/
void grow_large_tree1(float **xx, float **yy, int n, int *nnd, float *r, IntegerVector spv, NumericVector spvl, IntegerVector pt, IntegerVector dt,
                      int *npl, IntegerVector ncases, float **meany, NumericVector median)
{
  void initialize_tree1(int , int *, float *, IntegerVector , IntegerVector , NumericVector , float *);
  void split_a_node1(int , float **, float **, int ***, int **, int *, float *, IntegerVector , NumericVector, IntegerVector, IntegerVector, int *, int **, IntegerVector ,float **, NumericVector );
  int i, kk, t;
  void map_distx(float **, int , int , int **, int **), find_yvar(float *, int , float *, float *);
  void adjust_prior(float **y, int n);
  float hazard(float , float), tmp;
  void prune_tree(IntegerVector , int *, float *, int *);

  for(kk=0; kk<n_out; kk++)
    find_yvar(yy[kk], n, &sumy1[0][kk], &tmp);

  sumy1[0][1] = n - sumy1[0][1];
  //   adjust_prior(yy,n);
  initialize_tree1(n, npl, r, dt, spv, spvl, sumy1[0]);
  *nnd=1;
  map_distx(xx, n_meas, n, posi1[0], rank1);
  for(i=0; i<n; i++) casesInt[0][i]=i;
  ncases[0]=n;
  for(t=0, *nnd=1; t < *nnd; t++)
  {
    split_a_node1(t, xx, yy, posi1, rank1, nnd, r, spv, spvl, pt,
                  dt, npl, casesInt, ncases, sumy1, median);
    meany[t][0]=hazard(sumy1[t][0], sumy1[t][1]);
  }
}

void split_a_node1(int t, float **xx, float **yy, int ***posi, int **rank, int *nnd, float *r, IntegerVector spv, NumericVector spvl,
                   IntegerVector pt, IntegerVector dt, int *npl, int **cases_in_t, IntegerVector ncases, float **sumy, NumericVector median)
{
  float min_mr, *dvector(int ),  find_best_split1(int , float **, float **, int **, int , int *, float *, IntegerVector , NumericVector ,
                         int *, float **);
  void split_via_best(int , float **, int **, int ***, int , int *, IntegerVector , NumericVector ,
                      IntegerVector , IntegerVector , int *, int **, IntegerVector ), sort2(int , float *, int *);
  int n, yes, i, j, kk,  *ivector(int ), check_if_needed_split(float **, float **, int *, int , int , IntegerVector );
  int count_death=0;

  n=ncases[t];
  for(i=0; i<n; i++)
  {
    order[i]=i+1;
    j=cases_in_t[t][i];
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk]==1)
        medy[i+1]=yy[kk][j];
      else
        count_death+=yy[kk][j];
    }
  }

  death_catg[t]=n-count_death;
  sort2(n, medy, order);
  if(n%2==1)  median[t]=medy[n/2+1];
  else    median[t]=(medy[n/2]+medy[n/2+1])/2.0;

  yes = (*nnd < MAXNODES-1);
  if(yes) yes=check_if_needed_split(xx, yy, cases_in_t[t], n, t, pt);
  if(yes)
  {
    min_mr=find_best_split1(t, xx, yy, rank, n, nnd, r, spv, spvl,
                            cases_in_t[t], sumy);
    if(min_mr < MAXIMP)
    {
      split_via_best(t, xx, rank, posi, n, nnd, spv, spvl,
                     pt, dt, npl, cases_in_t, ncases);
    }
    else dt[t]=0;
  }
  else dt[t]=0;
  /*  free((char*)order);
  free((char*)medy);    */
}

float find_best_split1(int t, float **xx, float **yy, int **rank, int n, int *nnd, float *r, IntegerVector spv,NumericVector spvl,
                       int *cases_in_t, float **sumy)
{
  int i;
  float min_mr, r1, r2, minv, maxv;
  float ord_best_split1(float *, float **, int *, int *, int , float *, float *, float *, float *), thresh, it;
  float catg_best_split1(float *, float **, int , int *, int , float *, float *, float *, float *);
  void find_min_max(float *, int , float *, float *);
  float log_rank_catg(float *, float **, int, int *, int, float ), log_rank_ord(float *, float **, int *, int *, int, float );

  min_mr=MAXIMP;
  for(i=0; i<n_meas; i++)
  {
    if(catg[i])
    {
      if(catg[i] > 1)
      {
        find_min_max(xx[i], n_obs, &minv, &maxv);
        catg[i]=(int)maxv+1;

        it=catg_best_split1(xx[i], yy, catg[i], cases_in_t, n,
                            &r1, &r2, &thresh, sumy[t]);

      }
      else
      {
        it=ord_best_split1(xx[i], yy, rank[i], cases_in_t, n,
                           &r1, &r2, &thresh, sumy[t]);

      }
      if(it < min_mr)
      {
        min_mr=it;
        spv[t]=i;
        spvl[t]=thresh;
        r[*nnd]=r1, r[*nnd+1]=r2;
        sumy[*nnd][0]=sums[0][0], sumy[*nnd][1]=sums[0][1];
        sumy[*nnd+1][0]=sums[0][2], sumy[*nnd+1][1]=sums[0][3];
      }
    }
  }
  if(catg[spv[t]]>1)
  {
    lrank[t]=log_rank_catg(xx[spv[t]], yy, catg[spv[t]], cases_in_t, n, spvl[t]);
  }
  else
  {
    lrank[t]=log_rank_ord(xx[spv[t]], yy, rank[spv[t]], cases_in_t, n,spvl[t]);
  }
  return min_mr;
}

float ord_best_split1(float *x, float **yy, int *rank, int *cases_in_t, int n, float *r1, float *r2, float *thresh, float *sumy)
{
  float min_mr, imp;
  float loss(float *, float *, float *, float *, float );
  int i, j, jobs;
  int k, kk;
  void distinct_x(float *, int , int *, int *);
  void nrerror(const char *);

  for(i=0; i<n; i++)
  {
    j=cases_in_t[i];
    k=rank[j];
    tmpx[k]=x[j];
    if(k < 0 || k >= n)
    {
      Rprintf("rank %d out of bound\n", k);
      nrerror("exit from ord_best_split");
    }
    for(kk=0; kk<2; kk++)
    {
      ty[kk][k]=yy[kk][j];
    }
  }
  distinct_x(tmpx, n, &k, disx);
  min_mr=1.e10;
  jobs=0;
  if(k>1)
  {
    for(kk=0; kk<2; kk++)
      cntl1[kk]=0., cntr1[kk]=sumy[kk];
    for(i=1; i<k; i++)
    {
      for(j=disx[i-1]; j<disx[i]; j++)
      {
        cntl1[0] += ty[0][j], cntr1[0] -= ty[0][j];
        cntl1[1] += 1.-ty[1][j], cntr1[1] -= 1.-ty[1][j];
      }
      if(2*disx[i] > MINNODE && 2*(n-disx[i]) > MINNODE)
      {
        imp = loss(cntl1, cntr1, r1, r2, min_mr);
        if(imp < min_mr)
        {
          min_mr = imp;
          jobs = disx[i]-1;
          sums[0][0]=cntl1[0], sums[0][1]=cntl1[1];
          sums[0][2]=cntr1[0], sums[0][3]=cntr1[1];
        }
      }
    }
  }
  *thresh=tmpx[jobs];
  return min_mr;
}

float catg_best_split1(float *x, float **yy, int nxcat, int *p, int size, float *r1, float *r2, float *lcutc, float *sumy)

{
  float imp, timp;
  int i, j, k, kk, l, ucat=0, ulim, isave, nx, nl, nr;
  void find_costs(float **, float **, float *, float *, float *, int , int );
  float loss(float *, float *, float *, float *, float );

  timp=1.e10;
  isave=1;
  for(k=0;k<nxcat;k++)
  {
    xu[k]=xcat[k]=0; xdist[k]=0;
    ty[0][k]=ty[1][k]=0.;
  }
  for(kk=0; kk<n_out; kk++)
  { cntl1[kk]=0.;
    cntr1[kk]=sumy[kk];
  }
  nl=0;
  nr=size;
  for(i=0;i<size;i++)
  {
    j=p[i];
    if(j >= n_obs || j < 0) Rprintf("error obs ind\n");
    nx = (int) x[j];
    xu[nx]=1;
    xdist[nx]++;
    ty[0][nx] += yy[0][j], ty[1][nx] += 1.-yy[1][j];
  }
  for (i=0, ucat=0; i<nxcat; i++) if(xu[i]) xcat[ucat++]=i;
  if (ucat<=1)
    return(MAXIMP);
  ulim = (((unsigned)1)<<(ucat-1)) - 1;
  for(i=1; i<=ulim; i++)
  {
    l=i; k=0;
    while(!(l&1))
    {
      l>>=1; k++;
    }
    if(!(l&2))
    {
      nl+=xdist[xcat[k]], nr-=xdist[xcat[k]];
      cntl1[0]+= ty[0][xcat[k]], cntr1[0] -= ty[0][xcat[k]];
      cntl1[1]+= ty[1][xcat[k]], cntr1[1] -= ty[1][xcat[k]];
      /* nl++, nr--; */
    }
    else
    {
      nl-=xdist[xcat[k]], nr+=xdist[xcat[k]];
      cntl1[0] -= ty[0][xcat[k]], cntr1[0] += ty[0][xcat[k]];
      cntl1[1] -= ty[1][xcat[k]], cntr1[1] += ty[1][xcat[k]];
      /* nl--, nr++; */
    }
    if(2*nl > MINNODE && 2*nr > MINNODE)
    {
      imp = loss(cntl1, cntr1, r1, r2, timp);
      if(imp < timp)
      {
        isave = i;
        sums[0][0]=cntl1[0], sums[0][1]=cntl1[1];
        sums[0][2]=cntr1[0], sums[0][3]=cntr1[1];
        timp = imp;
      }
    }
  }
  *lcutc = isave;
  return(timp);
}

float loss(float *cl, float *cr, float *r1, float *r2, float min_mr)
{
  float lossl, lossr;
  void nrerror(const char *);

  lossl=lossr=0.;
  /*
  if(cl[1] <= 0. || cr[1] <= 0.){
  printf("%f %f\n", cl[0], cr[0]);
  nrerror("yt must be positive in loss");
  }
  */
  if(cl[1] > 0.)
    lossl = cl[1]-cl[1]*log(cl[1]/cl[0]);
  if(cr[1] > 0.)
    lossr = cr[1]-cr[1]*log(cr[1]/cr[0]);
  if(lossl+lossr < min_mr)
    *r1 = lossl, *r2 = lossr;
  return lossl+lossr;
}

void initialize_tree1(int n, int *npl, float *r, IntegerVector dt, IntegerVector spv, NumericVector spvl, float *sumy)
{
  int i;
  void nrerror(const char *);
  if(sumy[0] <= 0. || sumy[1] <= 0.)
  {
    Rprintf("%f %f\n", sumy[0], sumy[1]);
    nrerror("error in initialization");
  }
  else r[0]= sumy[1]-sumy[1]*log(sumy[1]/sumy[0]);
  for(i=0; i<n; i++) npl[i]=0;
  for(i=0; i<MAXNODES; i++)
  {
    dt[i]=0;
    spv[i]=0, spvl[i]=0.;
  }
}

float hazard(float yt, float dt)
{
  void nrerror(const char *);

  if(yt <= 0.) nrerror("yt must be positive in hazard");
  if(dt <= 0.) return .5/yt;
  return dt/yt;
}

/***********logrank*********/
void grow_large_tree2(float **xx, float **yy, int n, int *nnd, float *r, IntegerVector spv, NumericVector spvl,
                      IntegerVector pt, IntegerVector dt, int *npl, IntegerVector ncases, NumericVector median)
{
  void initialize_tree2(int , int *, IntegerVector , IntegerVector , NumericVector ), split_a_node2(int t, float **, float **, int ***, int **, int *, float *, IntegerVector , NumericVector ,
                        IntegerVector , IntegerVector , int *, int **, IntegerVector , float** , NumericVector );
  void find_yvar(float *, int , float *, float *);
  int i, t;
  float tmp;
  void map_distx(float **, int , int , int **, int **);

  find_yvar(yy[1],n,&sumy1[0][1],&tmp);
  sumy1[0][1]=n-sumy1[0][1];
  initialize_tree2(n, npl, dt, spv, spvl);
  r[0]=0.;
  *nnd=1;
  map_distx(xx, n_meas, n, posi1[0], rank1);
  for(i=0; i<n; i++) casesInt[0][i]=i;
  ncases[0]=n;
  for(t=0, *nnd=1; t < *nnd; t++)
  {
    split_a_node2(t, xx, yy, posi1, rank1, nnd, r, spv, spvl, pt,
                  dt, npl, casesInt, ncases, sumy1, median);

  }
}

void split_a_node2(int t, float **xx, float **yy, int ***posi, int **rank, int *nnd, float *r, IntegerVector spv, NumericVector spvl,
                   IntegerVector pt, IntegerVector dt, int *npl, int **cases_in_t, IntegerVector ncases, float** sumy, NumericVector median)
{
  float min_mr, *dvector(int );
  float find_best_split2(int , float **, float **, int **, int , int *, float *, IntegerVector , NumericVector, int *, float **);
  void split_via_best(int , float **, int **, int ***, int , int *, IntegerVector , NumericVector ,IntegerVector , IntegerVector , int *, int **, IntegerVector );
  void sort2(int , float *, int *);
  int n, yes, i, j, kk,  *ivector(int ), check_if_needed_split(float **, float **, int *, int , int , IntegerVector );
  int count_death=0;

  n=ncases[t];

  for(i=0; i<n; i++)
  {
    order[i]=i+1;
    j=cases_in_t[t][i];
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk]==1)
        medy[i+1]=yy[kk][j];
      else
        count_death+=yy[kk][j];
    }
  }
  death_catg[t]=n-count_death;
  sort2(n, medy, order);
  if(n%2==1)  median[t]=medy[n/2+1];
  else    median[t]=(medy[n/2]+medy[n/2+1])/2.0;

  yes = (*nnd < MAXNODES-1);
  if(yes) yes=check_if_needed_split(xx, yy, cases_in_t[t], n, t, pt);
  if(yes)
  {
    min_mr=find_best_split2(t, xx, yy, rank, n, nnd, r, spv, spvl, cases_in_t[t], sumy);
    if(min_mr > 0.)
    {
      split_via_best(t, xx, rank, posi, n, nnd, spv, spvl,
                     pt, dt, npl, cases_in_t, ncases);
      r[t]=min_mr;
    }
    else dt[t]=0, r[t]=0.;
  }
  else dt[t]=0, r[t]=0.;
  /*  free((char*)order);
  free((char*)medy);    */
}

float find_best_split2(int t, float **xx, float **yy, int **rank, int n, int *nnd, float *r, IntegerVector spv,
                       NumericVector spvl, int *cases_in_t, float **sumy)
{
  int i;
  float min_mr, minv, maxv, r1, r2;
  float ord_best_split2(float *, float **, int *, int *, int, float *, float *, float *, float *), thresh, it;
  float catg_best_split2(float *, float **, int, int *, int, float *, float *,  float *, float *);
  void find_min_max(float *, int , float *, float *);
  float log_rank_catg(float *, float **, int, int *, int, float),log_rank_ord(float *, float **, int *, int *, int, float);

  min_mr=0.;
  for(i=0; i<n_meas; i++)
  {
    if(catg[i])
    {
      if(catg[i] > 1)
      {
        find_min_max(xx[i], n_obs, &minv, &maxv);
        catg[i]=(int)maxv+1;
        it=catg_best_split2(xx[i], yy, catg[i], cases_in_t, n, &r1, &r2, &thresh, sumy[t]);
      }
      else
        it=ord_best_split2(xx[i], yy, rank[i], cases_in_t, n, &r1, &r2, &thresh, sumy[t]);
      if(it > min_mr)
      {
        min_mr=it;
        spv[t]=i;
        spvl[t]=thresh;
        sumy[*nnd][1]=sums[0][1];
        sumy[*nnd+1][1]=sums[0][3];
      }
    }
  }
  float tempvalue= spvl[t];
  if(catg[spv[t]]>1)
  { lrank[t]=log_rank_catg(xx[spv[t]], yy, catg[spv[t]], cases_in_t, n, tempvalue);
  }
  else
  {
    lrank[t]=log_rank_ord(xx[spv[t]], yy, rank[spv[t]], cases_in_t, n,tempvalue);
  }
  return min_mr;
}

void initialize_tree2(int n, int *npl, IntegerVector dt, IntegerVector spv, NumericVector spvl)
{
  int i;

  for(i=0; i<n; i++) npl[i]=0;
  for(i=0; i<MAXNODES; i++)
  {
    dt[i]=0;
    spv[i]=0, spvl[i]=0.;
  }
}

float ord_best_split2(float *x, float **yy, int *rank, int *cases_in_t, int n, float *r1, float *r2, float *thresh, float *sumy)
{
  float min_mr, imp;
  float logrank(int, int *, int *, int *, int, int);
  int i, ii, tt, t1, j, jobs;
  int dth, k, kk;
  void distinct_x(float *, int , int *, int *), sort2(int , float *, int *);
  void nrerror(const char *);

  for(i=0; i<n; i++)
  {
    j=cases_in_t[i];
    if(j < 0 || j >= n_obs)
    {
      Rprintf("wrong case %d\n", j);
      nrerror("exit from ord_best_split2");
    }
    k=rank[j];
    tmpx[k]=x[j];
    if(k < 0 || k >= n_obs)
    {
      Rprintf("rank %d out of bound for case %d\n", k, j);
      for(j=0; j<n_obs; j++) Rprintf("%d(%d, %d) ", rank[j], j, (int)x[j]);
      nrerror("exit from ord_best_split2");
    }
    for(kk=0; kk<2; kk++) ty[kk][k]=yy[kk][j];

  }
  distinct_x(tmpx, n, &k, disx);

  min_mr=0.;
  jobs=0;
  if(k>1)
  {
    for(j=t1=0; j<n; j++) if(ty[1][j] != 1.)
    {
      unceny[1+t1]=ty[0][j];
      indc[1+t1]=t1;
      t1++;
    }

    if(t1<=1) return min_mr;
    sort2(t1, unceny, indc);
    for(i=0; i<t1; i++) unceny[i] = unceny[i+1];
    distinct_x(unceny, t1, &tt, bin);

    if(tt<=1) return min_mr;
    for(j=0; j<tt; j++)
    {
      tr[j]=tl[j]=0;
      for(ii=0; ii<n; ii++) tr[j] += (ty[0][ii] >= unceny[bin[j]]);
    }
    dth=0;
    cntl1[1]=0., cntr1[1]=sumy[1];   //C-Y
    for(i=1; i<k; i++)
    {
      for(j=disx[i-1]; j<disx[i]; j++)
      {
        for(ii=0; ii<tt; ii++)
          tl[ii] += (ty[0][j] >= unceny[bin[ii]]);
        dth += (ty[1][j] == 0);
        cntl1[1] += 1.-ty[1][j], cntr1[1] -= 1.- ty[1][j];
      }

      if(2*disx[i] > MINNODE && 2*(n-disx[i]) > MINNODE)
      {

        imp = logrank(dth, tr, tl, bin, tt, t1);
        if(imp > min_mr)
        {
          min_mr = imp;
          jobs = disx[i]-1;
          sums[0][1]=cntl1[1]; /* should be sums[1] */
  sums[0][3]=cntr1[1];
        }
      }
    }
  }
  *thresh=tmpx[jobs];
  return min_mr;
}

float catg_best_split2(float *x, float **yy, int nxcat, int *p, int n, float *r1, float *r2,  float *lcutc, float *sumy)
{
  float min_mr, imp;
  float logrank(int, int *, int *, int *, int, int);
  int i, ii, l, tt, t1, j, *xnew, k, ucat=0, ulim, nx, isave;
  int *dth, death, **tm, nl, nr, *ivector(int ), **imatrix(int, int);
  void distinct_x(float *, int , int *, int *), sort2(int , float *, int *);
  void nrerror(const char *), free_imatrix(int **, int);

  min_mr=0.;
  isave=1;

  dth=ivector(nxcat);
  xnew=ivector(n);
  for(k=0; k<nxcat; k++)
  {
    xu[k]=xcat[k]=0; xdist[k]=dth[k]=0;
    ty[1][k]=0.; //C-Y
  }
  cntl1[1]=0.;
  cntr1[1]=0.;
  for(i=0, t1=0; i<n; i++)
  {
    j=p[i];
    if(j < 0 || j >= n_obs)
    {
      Rprintf("wrong case %d\n", j);
      nrerror("exit from catg_best_split2");
    }
    nx=(int)x[j];
    xu[nx]=1;
    xdist[nx]++;
    ty[1][nx] += 1.-yy[1][j];
    if(yy[1][j]!=1.) dth[nx]++;
    ty[1][i]=yy[1][j]; ty[0][i]=yy[0][j];
    xnew[i]=(int)x[j];
  }

  for (i=0, ucat=0; i<nxcat; i++)
    if(xu[i]) xcat[ucat++]=i;

    if(ucat<=1)
    {
      free((char*)xnew);
      free((char*)dth);
      return min_mr;
    }
    for(i=t1=0; i<n; i++)
    {
      if(ty[1][i]==0.)
      { unceny[1+t1]=ty[0][i];
        indc[1+t1]=t1;
        t1++;
      }

    }

    if(t1<=1)
    {
      free((char*)xnew);
      free((char*)dth);
      return min_mr;
    }
    sort2(t1, unceny, indc);
    for(i=0; i<t1; i++) unceny[i]=unceny[i+1];
    distinct_x(unceny, t1, &tt, bin);
    if(tt<=1)
    {
      free((char*)xnew);
      free((char*)dth);
      return min_mr;
    }
    tm=imatrix(nxcat, tt);

    for(j=0; j<tt; j++)
    {
      for(i=0; i<nxcat; i++) tm[i][j]=0;
    }

    for(j=0; j<tt; j++)
    {
      tr[j]=tl[j]=0;
      for(ii=0; ii<n; ii++)
      {
        for(i=0; i<ucat; i++) if(xnew[ii]==xcat[i])
          tm[xcat[i]][j]+=(ty[0][ii]>=unceny[bin[j]]);
      }
    }

    for(j=0; j<tt; j++)
      for(i=0; i<ucat; i++)
        tr[j]+=tm[xcat[i]][j];

    nl=0; nr=n;
    death=0;

    ulim=(((unsigned)1)<<(ucat-1))-1;

    for(i=1; i<=ulim; i++)
    {
      l=i; k=0;
      while(!(l&1))
      {
        l>>=1; k++;
      }

      if (!(l&2))
      {
        for(ii=0; ii<tt; ii++) tl[ii]+=tm[xcat[k]][ii];
        death+=dth[xcat[k]];
        nl+=xdist[xcat[k]], nr-=xdist[xcat[k]];
        cntl1[1]+=ty[1][xcat[k]], cntr1[1] -= ty[1][xcat[k]];
      }

      else
      {
        for(ii=0; ii<tt; ii++) tl[ii]-=tm[xcat[k]][ii];
        death-=dth[xcat[k]];
        nl-=xdist[xcat[k]], nr+=xdist[xcat[k]];
        cntl1[1]-=ty[1][xcat[k]], cntr1[1] += ty[1][xcat[k]];
      }

      if(2*nl > MINNODE && 2*nr > MINNODE)
      {
        imp = logrank(death, tr, tl, bin, tt, t1);
        if(imp > min_mr)
        {
          min_mr = imp;
          isave=i;
          sums[0][1]=cntl1[1], sums[0][3]=cntr1[1];
        }
      }
    }

    *lcutc=isave;
    free((char*)xnew);
    free((char*)dth);
    free_imatrix(tm, nxcat);
    return min_mr;
}

float log_rank_ord(float *x, float **yy, int *rank, int *cases_in_t, int n, float tempvalue)
{
  float logrank(int, int *, int *, int *, int, int), imp, min_mr;
  int i, ii, tt, t1, j, jobs;
  int dth, k, kk;
  void distinct_x(float *, int , int *, int *), sort2(int , float *, int *);
  void nrerror(const char *);

  for(i=0; i<n; i++)
  {
    j=cases_in_t[i];
    if(j < 0 || j >= n_obs)
    {
      Rprintf("wrong case %d\n", j);
      nrerror("exit from ord_best_split2");
    }
    k=rank[j];
    tmpx[k]=x[j];
    if(k < 0 || k >= n_obs)
    {
      Rprintf("rank %d out of bound for case %d\n", k, j);
      for(j=0; j<n_obs; j++) Rprintf("%d(%d, %d) ", rank[j], j, (int)x[j]);
      nrerror("exit from ord_best_split2");
    }
    for(kk=0; kk<2; kk++) ty[kk][k]=yy[kk][j];
  }

  distinct_x(tmpx, n, &k, disx);
  min_mr=0.;
  jobs=0;
  if(k>1)
  {
    for(j=t1=0; j<n; j++) if(ty[1][j] != 1.)
    {
      unceny[1+t1]=ty[0][j];
      indc[1+t1]=t1;
      t1++;
    }

    if(t1<=1) return min_mr;
    sort2(t1, unceny, indc);

    for(i=0; i<t1; i++) unceny[i] = unceny[i+1];

    distinct_x(unceny, t1, &tt, bin);
    if(tt<=1) return min_mr;

    for(j=0; j<tt; j++)
    {
      tr[j]=tl[j]=0;
      for(ii=0; ii<n; ii++) tr[j] += (ty[0][ii] >= unceny[bin[j]]);
    }
    dth=0;

    for(i=1; i<k; i++)
    {
      for(j=disx[i-1]; j<disx[i]; j++)
      {
        for(ii=0; ii<tt; ii++)
          tl[ii] += (ty[0][j] >= unceny[bin[ii]]);
        dth += (ty[1][j] == 0);
      }
      jobs = disx[i]-1;
      if(tempvalue==tmpx[jobs])
      {
        imp = logrank(dth, tr, tl, bin, tt, t1);
        return imp;
      }
    }
  }
  return min_mr;
}

float log_rank_catg(float *x, float **yy, int nxcat, int *p, int n, float tempvalue)
{
  float min_mr, imp;
  float logrank(int, int *, int *, int *, int, int);
  int i, ii, l, tt, t1, j, *xnew, k, ucat=0, nx;
  int *dth, death, **tm, nl, nr, *ivector(int ), **imatrix(int, int);
  void distinct_x(float *, int , int *, int *), sort2(int , float *, int *);
  void nrerror(const char *), free_imatrix(int **, int);

  min_mr=0.;
  dth=ivector(nxcat);
  xnew=ivector(n);
  for(k=0; k<nxcat; k++)
  {
    xu[k]=xcat[k]=0; xdist[k]=dth[k]=0;
  }

  for(i=0, t1=0; i<n; i++)
  {
    j=p[i];
    if(j < 0 || j >= n_obs)
    {
      Rprintf("wrong case %d\n", j);
      nrerror("exit from catg_best_split2");
    }
    nx=(int)x[j];
    xu[nx]=1;
    xdist[nx]++;
    if(yy[1][j]!=1.) dth[nx]++;
    ty[1][i]=yy[1][j]; ty[0][i]=yy[0][j];
    xnew[i]=(int)x[j];
  }

  for (i=0, ucat=0; i<nxcat; i++)
    if(xu[i]) xcat[ucat++]=i;

    if(ucat<=1)
    {
      free((char*)xnew);
      free((char*)dth);
      return min_mr;
    }
    for(i=t1=0; i<n; i++)
    {
      if(ty[1][i]==0.)
      { unceny[1+t1]=ty[0][i];
        indc[1+t1]=t1;
        t1++;
      }

    }

    if(t1<=1)
    {
      free((char*)xnew);
      free((char*)dth);
      return min_mr;
    }
    sort2(t1, unceny, indc);
    for(i=0; i<t1; i++) unceny[i]=unceny[i+1];
    distinct_x(unceny, t1, &tt, bin);
    if(tt<=1)
    {
      free((char*)xnew);
      free((char*)dth);
      return min_mr;
    }
    tm=imatrix(nxcat, tt);

    for(j=0; j<tt; j++)
    {
      for(i=0; i<nxcat; i++) tm[i][j]=0;
    }

    for(j=0; j<tt; j++)
    {
      tr[j]=tl[j]=0;
      for(ii=0; ii<n; ii++)
      {
        for(i=0; i<ucat; i++) if(xnew[ii]==xcat[i])
          tm[xcat[i]][j]+=(ty[0][ii]>=unceny[bin[j]]);
      }
    }

    for(j=0; j<tt; j++)
      for(i=0; i<ucat; i++)
        tr[j]+=tm[xcat[i]][j];

    nl=0; nr=n;
    death=0;

    /* ulim=(((unsigned)1)<<(ucat-1))-1; */

    for(i=1; i<=(int)tempvalue; i++)
    {
      l=i; k=0;
      while(!(l&1))
      {
        l>>=1; k++;
      }

      if (!(l&2))
      {
        for(ii=0; ii<tt; ii++) tl[ii]+=tm[xcat[k]][ii];
        death+=dth[xcat[k]];
        nl+=xdist[xcat[k]], nr-=xdist[xcat[k]];
      }

      else
      {
        for(ii=0; ii<tt; ii++) tl[ii]-=tm[xcat[k]][ii];
        death-=dth[xcat[k]];
        nl-=xdist[xcat[k]], nr+=xdist[xcat[k]];
      }
    }
    imp = logrank(death, tr, tl, bin, tt, t1);
    free((char*)xnew);
    free((char*)dth);
    free_imatrix(tm, nxcat);
    return imp;
}

/***Kaplan-Meier****/
void grow_large_tree3(float **xx, float **yy, int n, int *nnd, float *r, IntegerVector spv, NumericVector spvl, IntegerVector pt,
                      IntegerVector dt, int *npl, IntegerVector ncases, NumericVector median)
{
  void initialize_tree2(int, int *, IntegerVector, IntegerVector, NumericVector);
  void split_a_node3(int, float **, float **, int ***, int **, int *, float *, IntegerVector, NumericVector,IntegerVector, IntegerVector, int *, int **, IntegerVector, float **, NumericVector);
  void find_yvar(float *, int , float *, float *);
  int i, t;
  void map_distx(float **, int , int , int **, int **);
  float tmp;

  find_yvar(yy[1],n,&sumy1[0][1],&tmp);
  sumy1[0][1]=n-sumy1[0][1];
  initialize_tree2(n, npl, dt, spv, spvl);
  r[0]=0.;
  *nnd=1;
  map_distx(xx, n_meas, n, posi1[0], rank1);
  for(i=0; i<n; i++) casesInt[0][i]=i;
  ncases[0]=n;
  for(t=0, *nnd=1; t < *nnd; t++)
  {
    split_a_node3(t, xx, yy, posi1, rank1, nnd, r, spv, spvl, pt,
                  dt, npl, casesInt, ncases, sumy1, median);

  }
}

void split_a_node3(int t, float **xx, float **yy, int ***posi, int **rank, int *nnd, float *r, IntegerVector spv, NumericVector spvl,
                   IntegerVector pt, IntegerVector dt, int *npl, int **cases_in_t, IntegerVector ncases, float **sumy, NumericVector median)
{
  float min_mr, *dvector(int ), find_best_split3(int, float **, float **, int **, int, int *, float *, IntegerVector, NumericVector, int *, float **);
  void split_via_best(int , float **, int **, int ***, int , int *, IntegerVector , NumericVector ,
                      IntegerVector , IntegerVector , int *, int **, IntegerVector ), sort2(int , float *, int *);
  int n, yes, i, j, kk,  *ivector(int ), check_if_needed_split(float **, float **, int *, int , int , IntegerVector );
  int count_death=0;

  n=ncases[t];

  for(i=0; i<n; i++)
  {
    order[i]=i+1;
    j=cases_in_t[t][i];
    for(kk=0; kk<n_out; kk++)
    {
      if(n_class[kk]==1)
        medy[i+1]=yy[kk][j];
      else
        count_death+=yy[kk][j];
    }
  }
  death_catg[t]=count_death;
  sort2(n, medy, order);
  if(n%2==1)  median[t]=medy[n/2+1];
  else    median[t]=(medy[n/2]+medy[n/2+1])/2.0;

  yes = (*nnd < MAXNODES-1);
  if(yes) yes=check_if_needed_split(xx, yy, cases_in_t[t], n, t, pt);
  if(yes)
  {
    min_mr=find_best_split3(t, xx, yy, rank, n, nnd, r, spv, spvl, cases_in_t[t], sumy);
    if(min_mr > 0.)
    {
      split_via_best(t, xx, rank, posi, n, nnd, spv, spvl,
                     pt, dt, npl, cases_in_t, ncases);
      r[t]=min_mr;
    }
    else dt[t]=0, r[t]=0.;
  }
  else dt[t]=0, r[t]=0.;
  /*  free((char*)order);
  free((char*)medy);  */
}

float find_best_split3(int t, float **xx, float **yy, int ** rank, int n, int *nnd, float *r, IntegerVector spv,
                       NumericVector spvl, int *cases_in_t, float **sumy)
{
  int i;
  float min_mr, minv, maxv, r1, r2;
  float ord_best_split3(float *, float **, int *, int *, int, float *, float *, float *, float *), thresh, it;
  float catg_best_split3(float *, float **, int, int *, int, float *, float *, float *, float *);
  void find_min_max(float *, int , float *, float *);
  float log_rank_catg(float *, float **, int, int *, int, float ), log_rank_ord(float *, float **, int *, int *, int, float);

  min_mr=0.;
  for(i=0; i<n_meas; i++)
  {
    if(catg[i])
    {
      if(catg[i] > 1)
      {
        find_min_max(xx[i], n_obs, &minv, &maxv);
        catg[i]=(int)maxv+1;
        it=catg_best_split3(xx[i], yy, catg[i], cases_in_t, n, &r1, &r2, &thresh, sumy[t]);
      }
      else
        it=ord_best_split3(xx[i], yy, rank[i], cases_in_t, n, &r1, &r2, &thresh, sumy[t]);
      if(it > min_mr)
      {
        min_mr=it;
        spv[t]=i;
        spvl[t]=thresh;
        sumy[*nnd][1]=sums[0][1];
        sumy[*nnd+1][1]=sums[0][3];
      }
    }
  }
  if(catg[spv[t]]>1)
  { lrank[t]=log_rank_catg(xx[spv[t]], yy, catg[spv[t]], cases_in_t, n,
                           spvl[t]);
  }
  else
  {
    lrank[t]=log_rank_ord(xx[spv[t]], yy, rank[spv[t]], cases_in_t, n, spvl[t]);
  }
  return min_mr;
}

void prune_tree(IntegerVector dt, int *npl, float *lrank, int *nnd)
{
  int i, j, k=0, t1, t2;

  beforeprune = *nnd;
  for(i=*nnd-1; i>=0; i--)
  {
    if(dt[i] && lrank[i]<=4)
    {
      t1 = dt[i], t2 = dt[i]+1;
      if(dt[t1]==0 && dt[t2]==0)
      {
        dt[i]=0, k++;
        for(j=0; j<n_obs; j++)
        {
          if(npl[j]==t1 || npl[j]==t2) npl[j]=i;
        }
      }
    }
  }
  *nnd -= 2*k;
}

float logrank(int dth, int *ni, int *n1i, int *mi, int tt, int t1)
{
  void nrerror(const char *);
  float a1, a2;
  int m1i, i;

  a1 = (float)dth, a2=0.;
  for(i=0; i<tt; i++)
  {
    m1i= (i < tt-1) ? mi[i+1]-mi[i] : t1-mi[tt-1];
    if(ni[i] > 0)
      a1 -= 1.*m1i*n1i[i]/ni[i];
    if(ni[i] > 1)
      a2 += 1.*m1i*(ni[i]-m1i)*n1i[i]*(ni[i]-n1i[i])/(ni[i]*ni[i]*(ni[i]-1.));
  }
  if(a2 <= 0.)
  {
    if(a2 < -1.e-6) Rprintf("warning: negative variace %f\n", a2);
    return 0.;
  }
  return a1*a1/a2;
}

float km_dist(int *ni, int *n1i, int *mi, int *m2i, int tt)
{
  void nrerror(const char *);
  float a1, a2, kmd, tmp9;
  int m1i, i;

  if(2*ni[0] <= MINNODE || 2*(ni[0]-n1i[0]) <= MINNODE) return 0.;
  a1 = a2 = 1.;
  kmd=0.;
  for(i=0; i<tt-1; i++)
  {
    m1i= mi[i+1]-mi[i];
    if(m1i > 0 && n1i[i] > 0 && ni[i] > n1i[i])
    {
      if(m1i < m2i[i])
      {
        Rprintf("%d %d %f\n", m1i, m2i[i], unceny[mi[i]+1]);
        nrerror("error1 in km\n");
      }
      if(m1i - m2i[i] > ni[i]-n1i[i]) nrerror("error2 in km\n");
      if(m2i[i] > n1i[i]) nrerror("error3 in km\n");
      a1 *= 1. - m2i[i]*1./n1i[i];
      a2 *= 1. - (m1i - m2i[i])*1./(ni[i]-n1i[i]);
      tmp9 = (unceny[mi[i+1]]-unceny[mi[i]])*(a1-a2);
      if(tmp9 > 0.) kmd += tmp9;
      else kmd -= tmp9;
    }
  }
  return kmd;
}

float ord_best_split3(float *x, float **yy, int *rank, int *cases_in_t, int n, float *r1, float *r2, float *thresh, float *sumy)
{
  float min_mr, imp;
  float km_dist(int *, int *, int *, int *, int);
  int i, ii, tt, t1, j, jobs;
  int k, kk;
  void distinct_x(float *, int , int *, int *), sort2(int , float *, int *);
  void nrerror(const char *);

  for(i=0; i<n; i++)
  {
    j=cases_in_t[i];
    if(j < 0 || j >= n_obs)
    {
      Rprintf("wrong case %d\n", j);
      nrerror("exit from ord_best_split3");
    }
    k=rank[j];
    tmpx[k]=x[j];
    if(k < 0 || k >= n_obs)
    {
      Rprintf("rank %d out of bound for case %d\n", k, j);
      for(j=0; j<n_obs; j++) Rprintf("%d(%d, %d) ", rank[j], j, (int)x[j]);
      nrerror("exit from ord_best_split3");
    }
    for(kk=0; kk<2; kk++)
      ty[kk][k]=yy[kk][j];
  }
  distinct_x(tmpx, n, &k, disx);
  min_mr=0.;
  jobs=0;
  if(k>1)
  {
    for(j=t1=0; j<n; j++) if(ty[1][j] != 1.)
    {
      unceny[1+t1]=ty[0][j];
      indc[1+t1]=t1;
      t1++;
    }
    if(t1<=1) return min_mr;
    sort2(t1, unceny, indc);
    for(i=0; i<t1; i++) unceny[i] = unceny[i+1]; /*changed*/
  distinct_x(unceny, t1, &tt, bin);
  if(tt<=1) return min_mr;
  for(j=0; j<tt; j++)
  {
    tr[j]=tl[j]=indc[j]=0;
    for(ii=0; ii<n; ii++) tr[j] += (ty[0][ii] >= unceny[bin[j]]);
  }
  cntl1[1]=0., cntr1[1]=sumy[1];
  for(i=1; i<k; i++)
  {
    for(j=disx[i-1]; j<disx[i]; j++)
    {
      for(ii=0; ii<tt; ii++)
      {
        tl[ii] += (ty[0][j] >= unceny[bin[ii]]);
        indc[ii] += (ty[0][j] == unceny[bin[ii]]) && (ty[1][j]==0);
      }
      cntl1[1] += 1.-ty[1][j], cntr1[1] -= 1.- ty[1][j];
    }
    if(2*disx[i] > MINNODE && 2*(n-disx[i]) > MINNODE)
    {
      imp = km_dist(tr, tl, bin, indc, tt);
      if(imp > min_mr)
      {
        min_mr = imp;
        jobs = disx[i]-1;
        sums[0][1]=cntl1[1];
        sums[0][3]=cntr1[1];
      }
    }
  }
  }
  *thresh=tmpx[jobs];
  return min_mr;
}

float catg_best_split3(float *x, float **yy, int nxcat, int *p, int n, float *r1, float *r2, float *lcutc, float *sumy)
{
  float min_mr, imp;
  float km_dist(int *, int *, int *, int*, int);
  int i, ii, l, tt, t1, j, *xnew, k, ucat=0, ulim, nx, isave;
  int **tm, **indm, nl, nr, *ivector(int ), **imatrix(int, int);
  void distinct_x(float *, int , int *, int *), sort2(int, float *, int *);
  void nrerror(const char *);
  void free_imatrix(int **, int);

  min_mr=0.;
  isave=1;

  xnew=ivector(n);
  for(k=0; k<nxcat; k++)
  {
    xu[k]=xcat[k]=0; xdist[k]=0;
    ty[1][k]=0.;
  }
  cntl1[1]=0.;
  cntr1[1]=0.;
  for(i=0, t1=0; i<n; i++)
  {
    j=p[i];
    if(j < 0 || j >= n_obs)
    {
      Rprintf("wrong case %d\n", j);
      nrerror("exit from catg_best_split2");
    }
    nx=(int)x[j];
    xu[nx]=1;
    xdist[nx]++;
    ty[1][i]=yy[1][j]; ty[0][i]=yy[0][j];
    xnew[i]=(int)x[j];
  }

  for (i=0, ucat=0; i<nxcat; i++)
    if(xu[i]) xcat[ucat++]=i;

    if(ucat<=1)
    {
      free((char*)xnew);
      return min_mr;
    }
    for(i=t1=0; i<n; i++)
    {
      if(ty[1][i]==0.)
      { unceny[1+t1]=ty[0][i];
        indc[1+t1]=t1;
        t1++;
      }

    }
    if(t1<=1)
    {
      free((char*)xnew);
      return min_mr;
    }
    sort2(t1, unceny, indc);
    for(i=0; i<t1; i++) unceny[i]=unceny[i+1];
    distinct_x(unceny, t1, &tt, bin);
    if(tt<=1)
    {
      free((char*)xnew);
      return min_mr;
    }

    tm=imatrix(nxcat, tt);
    indm=imatrix(nxcat, tt);

    for(j=0; j<tt; j++)
    {
      for(i=0; i<nxcat; i++)
      {
        tm[i][j]=0; indm[i][j]=0;
      }
    }

    for(j=0; j<tt; j++)
    {
      tr[j]=tl[j]=indc[j]=0;
      for(ii=0; ii<n; ii++)
      {
        for(i=0; i<ucat; i++) if(xnew[ii]==xcat[i])
        {
          tm[xcat[i]][j]+=(ty[0][ii]>=unceny[bin[j]]);
          indm[xcat[i]][j]+=(ty[0][ii]==unceny[bin[j]])&&(ty[1][ii]==0);
        }
      }
    }

    for(j=0; j<tt; j++)
      for(i=0; i<ucat; i++)
        tr[j]+=tm[xcat[i]][j];

    nl=0; nr=n;

    ulim=(((unsigned)1)<<(ucat-1))-1;

    for(i=1; i<=ulim; i++)
    {
      l=i; k=0;
      while(!(l&1))
      {
        l>>=1; k++;
      }

      if (!(l&2))
      {
        for(ii=0; ii<tt; ii++)
        {
          tl[ii]+=tm[xcat[k]][ii];
          indc[ii]+=indm[xcat[k]][ii];
        }
        nl+=xdist[xcat[k]], nr-=xdist[xcat[k]];
        cntl1[1]+=ty[1][xcat[k]], cntr1[1] -= ty[1][xcat[k]];
      }

      else
      {
        for(ii=0; ii<tt; ii++)
        {
          tl[ii]-=tm[xcat[k]][ii];
          indc[ii]-=indm[xcat[k]][ii];
        }
        nl-=xdist[xcat[k]], nr+=xdist[xcat[k]];
        cntl1[1]-=ty[1][xcat[k]], cntr1[1] += ty[1][xcat[k]];
      }

      if(2*nl > MINNODE && 2*nr > MINNODE)
      {
        imp = km_dist(tr, tl, bin, indc, tt);
        if(imp > min_mr)
        {
          min_mr = imp;
          isave=i;
          sums[0][1]=cntl1[1], sums[0][3]=cntr1[1];
        }
      }
    }
    free((char*)xnew);
    free_imatrix(indm, nxcat);
    free_imatrix(tm, nxcat);
    *lcutc=(float)isave;
    return min_mr;
}

int *find_ori_catg(float *covar, int *k)
{
  void sort2(int, float *, int *);
  int *idistinct_x(float *, int, int *);
  float *ftmp;
  int *itmp, *distx, i;

  itmp=(int*)malloc((unsigned)(n_obs+1)*sizeof(int));
  ftmp=(float*)malloc((unsigned)(n_obs+1)*sizeof(float));
  /* itmp = ivector(n_obs+1); */
  /* ftmp = dvector(n_obs+1); */
  for(i=0; i<n_obs; i++)
  {
    itmp[i+1] = 1;
    ftmp[i+1] = covar[i];
  }
  sort2(n_obs, ftmp, itmp);
  for(i=0; i<n_obs; i++)
  {
    ftmp[i]=ftmp[i+1];
    itmp[i]=itmp[i+1];
  }
  free((char*)itmp);
  itmp=idistinct_x(ftmp, n_obs, k);
  distx=(int*)malloc((unsigned)(*k)*sizeof(int));

  /* distx = ivector(*k); */
  for(i=0; i<*k; i++)
    distx[i] = (int) ftmp[itmp[i]-1];
  free((char*)ftmp);
  return distx;
}

void nrerror(const char *error_text)
  /*numerical recipes standard error handler*/
{
  //fprintf(stderr, "Numerical Recipes run-time error ...\n");
  //fprintf(stderr, "%s\n", error_text);
  //fprintf(stderr, "...now exiting to system...\n");
  //fprintf(stderr, "Press Enter to end this program!");
  getchar();
  //exit(1);
}

float **matrix(int nr, int nc)
{
  int i;
  float **m;
  void nrerror(const char *);
  /* unsigned malloc(); */

  m=(float **) malloc((unsigned) nr*sizeof(float*));
  if (!m) nrerror("error in matrix");
  for(i=0;i<nr;i++)
  {
    m[i]=(float *) malloc((unsigned) nc*sizeof(float));
    if (!m[i]) nrerror("error in matrix");
  }
  return m;
}

void free_matrix(float **m, int nr)
{
  int i;

  for(i=nr-1;i>=0;i--) free((char*)m[i]);
  free((char*)m);
}

int **imatrix(int nr, int nc)
{
  int i, j, **m;
  void nrerror(const char *);

  m=(int **) malloc((unsigned) nr*sizeof(int*));
  if (!m) nrerror("error in imatrix");
  for(i=0;i<nr;i++)
  {
    m[i]=(int *) malloc((unsigned) nc*sizeof(int));

    if (!m[i])
      nrerror("error in imatrix");
    else
      for(j=0;j<nc;j++)
        m[i][j]=0;
  }
  return m;
}

void free_imatrix(int **m, int nr)
{
  int i;

  for(i=nr-1;i>=0;i--) free((char*)m[i]);
  free((char*)m);
}

float *dvector(int nh)
  /*allocate a float vector with range [nl, nh]*/
{
  float *v;
  int i;
  void nrerror(const char *);

  v=(float *)malloc((unsigned)nh*sizeof(float));
  if(!v)
    nrerror("allocation failure in dvector()");
  else
    for(i=0;i<nh;i++)
      v[i]=0.0;
  return v;
}

int *ivector(int nh)
  /*allocate a int vector with range [nl, nh]*/
{
  int *v,i;
  void nrerror(const char *);

  v=(int *)malloc((unsigned)nh*sizeof(int));
  if(!v)
    nrerror("allocation failure in ivector()");
  else
    for(i=0;i<nh;i++)
      v[i]=0;  /* initialize */
  return v;
}

float **matrix2(int nr, int *nc)
{
  int i;
  float **m;
  void nrerror(const char *);

  m=(float **) malloc((unsigned) nr*sizeof(float*));
  if (!m) nrerror("error in matrix");
  for(i=0;i<nr;i++)
  {
    m[i]=(float *) malloc((unsigned) (nc[i]+1)*sizeof(float));
    if (!m[i]) nrerror("error in dmatrix");
  }
  return m;
}

int *idistinct_x(float *tmpx, int n_cases, int *k)
{
  int i, j, *disx;

  disx=ivector(n_cases);
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


void input(DataFrame dataset, IntegerVector colname)
{
  int num = colname.size();
  NumericVector coltemp;
  float tmp, **matrix(int , int );
  int i, j, k, kk, ncol=0, *ivector(int ), **imatrix(int , int );
  int *cflag, sub_obs, minv0, maxv0, adjust=0, itmp;
  void clean_y();
  float **matrix2(int , int *), *dvector(int );
  void nrerror(const char *);
  void find_min_max(float *, int , float *, float *);
  float minv, maxv;
  int *find_ori_catg(float *,int *);

  n_out=0;
  n_meas=0;
  ncol=0;

  for(i=0;i<num;i++)
  {
    if(colname[i]==0||colname[i]==1||colname[i]==2)
    {
      ncol+=1;
      if(colname[i]==1||colname[i]==2)
        n_meas+=1;
    }
    else if(colname[i]==-1||colname[i]==-2)
      n_out+=1;
    else
    {
      Rprintf("abnormal exit");
      //exit(-1);
    }
  }

  if(n_meas<1)
  {
    Rprintf("ERROR: No covariate!");
    //exit(-1);
  }

  coltemp = dataset[0];
  sub_obs = coltemp.size();
  n_class = ivector(n_out);
  cflag = ivector(ncol+n_out);
  catg = ivector(2*n_meas);
  maxcatg = ivector(2*n_meas);
  for (i=0; i<ncol+n_out; i++)
    cflag[i]=colname[i];
  y=matrix(n_out, sub_obs);
  X=matrix(2*n_meas, sub_obs);
  indi=imatrix(n_meas, sub_obs);
  indj=imatrix(n_out, sub_obs);
  for(j=0,k=0,kk=0; j<ncol+n_out; j++)
  {
    coltemp = dataset[j];
    if(cflag[j]>0)
    {
      for(i=0; i<sub_obs; i++)
      {
        indi[k][i]=0;
        X[k][i]=coltemp[i];
      }
      catg[k]=(cflag[j] > 1) + 1;
      k++;
    }
    else if(cflag[j]<0)
    {
      for(i=0; i<sub_obs; i++)
      {
        if (cflag[j]== -1 && kk==0)
          adjust=1;
        indj[kk][i]=false;
        y[kk][i]=coltemp[i];
      }
      n_class[kk]=1+(cflag[j] == -1);
      kk++;
    }
  }


  if (adjust ==1)
  {
    for (i=0; i< sub_obs; i++)
    {
      tmp=y[0][i];
      y[0][i]=y[1][i];
      y[1][i]=tmp;
    }
    itmp=n_class[0];
    n_class[0]=n_class[1];
    n_class[1]=itmp;
  }

  for(i=0, j=0; i<sub_obs; i++)
  {
    if((indj[0][i]==false) && (indj[1][i]==false))
    {
      for(k=0; k<n_meas; k++)
      {
        X[k][j]= X[k][i];
        indi[k][j]=indi[k][i];
        if((indi[k][j]==1) && catg[k]==1)
          X[k][j]=X[k][j-1];
      }
      for(kk=0; kk<n_out; kk++)
        y[kk][j]=y[kk][i];
      j++;
    }
  }

  n_obs=j;
  for(i=0; i<n_obs; i++)
    y[1][i] = (1.0-y[1][i]);

  var_to_var = (int*) calloc (2*n_meas, sizeof(int));
  if((ori_catg = (int **)calloc(n_meas, sizeof(int*)))==NULL)
    nrerror("error calloc ori_catg\n");

  for(j=0, k=0; j<n_meas; j++)
  {
    var_to_var[j]=j;
    find_min_max(X[j], n_obs, &minv, &maxv);
    find_min_max((float *)indi[j], n_obs, (float *)&minv0, (float *)&maxv0);
    if(maxv <= minv)
      catg[j]=0;
    else if(catg[j] > 1)
    {
      catg[j]=(int)maxv+1;
      maxcatg[j]=catg[j];
      ori_catg[j] = find_ori_catg(X[j], &catg[j]);
    }
    else if(catg[j]==1 && maxv0==1)
    {
      var_to_var[n_meas+k]=j;
      catg[n_meas+k]=catg[j];
      for(i=0; i<n_obs; i++)
      {
        if(indi[j][i]==0)
        {
          X[j][i]=X[j][i];
          X[n_meas+k][i]=X[j][i];
        }
        else
        {
          X[j][i]=minv-1;
          X[n_meas+k][i]=maxv-1;
        }
      }
      k++;
    }
  }
  n_meas=n_meas+k;

  prior=matrix2(n_out,n_class);
  ssy=dvector(n_out);
  yprop=matrix2(n_out,n_class);
  posi1=(int***)malloc((unsigned)MAXNODES*sizeof(int**));
  if(!posi1)
    nrerror("error in allocate posi");
  for(i=0; i<MAXNODES; i++)
    posi1[i]=imatrix(n_meas, n_obs+1);
  casesInt=imatrix(MAXNODES, n_obs+1);
  rank1=imatrix(n_meas, n_obs+1);
  sumy1 = matrix(MAXNODES, 2*n_out);
  sums = matrix(n_out, 4);
  cntl1=dvector(2), cntr1=dvector(2);
  cntl2=matrix(2,2), cntr2=matrix(2,2);
  disx=ivector(n_obs), tl=ivector(n_obs), tr=ivector(n_obs);
  tmpx=dvector(n_obs+1);
  cp=matrix(2,2);
  c1=ivector(n_out), c2=ivector(n_out);
  ty=matrix(n_out, n_obs+1), xdist=ivector(100);
  xu=ivector(100), xcat=ivector(100), xdist=ivector(100), bin=ivector(n_obs+1);
  counts=imatrix(n_out, 200), ynum=imatrix(n_out, 200);
  ycnt=imatrix(n_out, 200), ycat=imatrix(n_out, 200);
  indc=ivector(n_obs+1), tmp12=matrix(n_out, 200);
  unceny=dvector(n_obs+1);
  clean_y();
}
