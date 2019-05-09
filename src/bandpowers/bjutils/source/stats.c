// statistics functions

#include "../include/bjutils.h"


void random_gauss(double *arr,int size,double mean,double stddev)
{
  int i;
  const gsl_rng_type *T=gsl_rng_default;
  gsl_rng_default_seed=1;                      
  gsl_rng *rnd=gsl_rng_alloc(T);

  for (i=0;i<size;i++) {
    arr[i]=mean+gsl_ran_gaussian(rnd,stddev);
  }
  gsl_rng_free(rnd);
  return;
}


void gaussian_correlated(double **res,int dim,int nsample,double *mean,gsl_matrix *cov)
{
  int i,j,k;
  double sum;
  gsl_vector *varu=gsl_vector_alloc(dim);
  gsl_matrix *corr=gsl_matrix_alloc(dim,dim);

  // initialise random sampling
  const gsl_rng_type *T=gsl_rng_default;
  gsl_rng_default_seed=1;                      
  gsl_rng *rnd=gsl_rng_alloc(T);

  // compute correlation matrix and decompose
  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      gsl_matrix_set(corr,i,j,gsl_matrix_get(cov,i,j)/sqrt(gsl_matrix_get(cov,i,i)*gsl_matrix_get(cov,j,j)));
    }
  }
  gsl_linalg_cholesky_decomp(corr);

  // sample
  for (i=0;i<nsample;i++) {
    for (j=0;j<dim;j++) {
      gsl_vector_set(varu,j,gsl_ran_gaussian(rnd,sqrt(gsl_matrix_get(cov,j,j))));  // uncorrelated variates
    }

    for (j=0;j<dim;j++) {
      sum=0.0;
      for (k=0;k<=j;k++) {
	sum+=gsl_matrix_get(corr,j,k)*gsl_vector_get(varu,k);  // correlate
      }
      res[j][i]=mean[j]+sum;
    }
  }

  // clean up
  gsl_rng_free(rnd);
  gsl_matrix_free(corr);
  gsl_vector_free(varu);
  return;
}


// same as gaussian_correlated but provide random sampler externally
void gaussian_correlated_withrnd(double **res,int dim,int nsample,double *mean,gsl_matrix *cov,gsl_rng *rnd)
{
  int i,j,k;
  double sum;
  gsl_vector *varu=gsl_vector_alloc(dim);
  gsl_matrix *corr=gsl_matrix_alloc(dim,dim);

  // compute correlation matrix and decompose
  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      gsl_matrix_set(corr,i,j,gsl_matrix_get(cov,i,j)/sqrt(gsl_matrix_get(cov,i,i)*gsl_matrix_get(cov,j,j)));
    }
  }
  gsl_linalg_cholesky_decomp(corr);

  // sample
  for (i=0;i<nsample;i++) {
    for (j=0;j<dim;j++) {
      gsl_vector_set(varu,j,gsl_ran_gaussian(rnd,sqrt(gsl_matrix_get(cov,j,j))));  // uncorrelated variates
    }

    for (j=0;j<dim;j++) {
      sum=0.0;
      for (k=0;k<=j;k++) {
	sum+=gsl_matrix_get(corr,j,k)*gsl_vector_get(varu,k);  // correlate
      }
      res[j][i]=mean[j]+sum;
    }
  }

  // clean up
  gsl_matrix_free(corr);
  gsl_vector_free(varu);
  return;
}


gsl_matrix *sample_wishart(gsl_matrix *mean,int dof)
{
  int i,j,k,l;
  double v,sum;
  const int dim=mean->size1;
  gsl_matrix *sample=gsl_matrix_calloc(dim,dim);
  gsl_matrix *decomp=gsl_matrix_calloc(dim,dim);
  gsl_matrix *res=gsl_matrix_calloc(dim,dim);

  static int flag=0;
  static gsl_rng *r1,*r2;
  if (!flag) {
    r1=gsl_rng_alloc(gsl_rng_default);
    r2=gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r1,1);
    gsl_rng_set(r2,2);
  }

  for(i=0;i<dim;i++) {
    sum=0.0;
    for(k=0;k<i;k++) {
      sum+=gsl_matrix_get(res,k,i)*gsl_matrix_get(res,k,i);
    }
    v=gsl_ran_chisq(r2,dof-i);

    gsl_matrix_set(sample,i,i,sum+v);
    for(j=i+1;j<dim;j++) {
      sum=0.0;
      for(k=0;k<i;k++) {
	sum+=gsl_matrix_get(res,k,i)*gsl_matrix_get(res,k,j);
      }
      gsl_matrix_set(res,i,j,gsl_ran_gaussian(r1,1.)); //Gaussian samples in res
      gsl_matrix_set(sample,i,j,gsl_matrix_get(res,i,j)*sqrt(v)+sum);
      gsl_matrix_set(sample,j,i,gsl_matrix_get(sample,i,j));  //symmetric
    }
  }
  gsl_matrix_scale(sample,1./(1.*dof));  //return unbiased estimate of mean

  gsl_matrix_memcpy(decomp,mean);
  gsl_linalg_cholesky_decomp(decomp);
  for(i=0;i<dim;i++) {
    for(j=i;j<dim;j++) {  // resulting matrix is symmetric
      sum=0.0;
      for(k=0;k<=i;k++) {  // only lower diagonal
	for(l=0;l<=j;l++) {  // only upper diagonal
	  sum+=gsl_matrix_get(decomp,i,k)*gsl_matrix_get(sample,k,l)*gsl_matrix_get(decomp,l,j);
	}
      }
      gsl_matrix_set(res,i,j,sum);
      gsl_matrix_set(res,j,i,sum);
    }                     
  }

  flag=1;
  gsl_matrix_free(sample);
  gsl_matrix_free(decomp);
  return(res);
}



// currently same outbinning of smooth curve as input binning
void smoothingsplines(int n,double *x,double *yarray,double *warray,int order,int nbreak, double *fit)
{
  int i,j;
  gsl_vector *c, *w,*B, *y;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;
  gsl_bspline_workspace *bw;
  double chisq,dof,bj,yerr; 

  // allocate a cubic bspline workspace
  bw=gsl_bspline_alloc((size_t)order,nbreak);
  size_t ncoeffs=gsl_bspline_ncoeffs(bw);
  B=gsl_vector_alloc(ncoeffs);
    
  // further allocations 
  y=gsl_vector_alloc(n);
  X=gsl_matrix_alloc(n,ncoeffs);
  c=gsl_vector_alloc(ncoeffs);
  w=gsl_vector_alloc(n);
  cov=gsl_matrix_alloc(ncoeffs,ncoeffs);
  mw=gsl_multifit_linear_alloc(n,ncoeffs);

  for (i=0;i<n;i++) {
    gsl_vector_set(y,i,yarray[i]);
    gsl_vector_set(w,i,warray[i]);
  }

  // use uniform breakpoints on [XMIN,XMAX] 
  double xmin=x[0];
  double xmax=x[n-1];
  gsl_bspline_knots_uniform(xmin,xmax,bw);

  //construct the fit matrix X 
  for (i=0;i<n;i++) {
    gsl_bspline_eval(x[i],B,bw);      // compute B_j(xi) for all j 
    for (j=0;j<ncoeffs;j++) {       // fill in row i of X 
      bj=gsl_vector_get(B,j);
      gsl_matrix_set(X,i,j,bj);
    }
  }
     
  // perform the fit 
  gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);     
  dof=n-ncoeffs;
  fprintf(stderr,"Spline fit: chisq/dof = %.3f\n",chisq/dof);

  // obtain result
  for (i=0;i<n;i++) {     // can be different from input binning
    gsl_bspline_eval(x[i],B,bw);
    gsl_multifit_linear_est(B,c,cov,&fit[i],&yerr);
  }

  // clean up
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  return;
}


// boxcar smoothing
void boxcar(double *sample_in,double *sample_out,int size,int width)
{
  int i,j;
  double norm=2.*width+1.;

  if (width==0) return;  //nothing to be done
  if ((width<0)||(width>size/2)) {
    printf("Error in routine boxcar: width parameter should be between %i and %i\n",0,size/2);
    exit(-1);
  }

  for (i=0;i<size;i++) { 
    if ((i-width>=0)&&(i+width<size)) {
      sample_out[i]=0.0;
      for (j=i-width;j<=i+width;j++) { 
	sample_out[i]+=sample_in[j];
      }
      sample_out[i]/=norm;
    }
    else {
      sample_out[i]=sample_in[i];
    }
  }
  return;
}


// to makehisto
void bisect_makehisto(double *a,int dim,double x,int *erg)
{
  int l=0,u=dim-1,m;

  if ((x<a[l]) || (x>a[u])) {
    printf("error: index out of range\n");
    exit(-1);
  }
  while (u-l>1) {
    m=(u+l)/2;
    if (x>a[m]) l=m;
    else u=m;
  }
  *erg=l;
}

void makehisto(double *sample,int size,double *x,double *hist,int *nbin,double minfix,double maxfix)
{
  int i,j,ix,NBIN,ex;
  double min=1.e10,max=-1.e10,scale=1.,width;

  // define binning
  if (nbin[0]==0) {         // automatic binning
    for (i=0;i<size;i++) {
      if (min>sample[i]) min=sample[i];
      if (max<sample[i]) max=sample[i];
    }
    scale=pow(10,floor(log10(fabs(max-min))));
    min=scale*floor(min/scale);
    max=scale*ceil(max/scale);
    NBIN=(int)floor(sqrt(size));
    nbin[0]=NBIN;
  }
  else {
    NBIN=nbin[0];
    min=minfix;
    max=maxfix;
  }
  width=(max-min)/(1.*NBIN);
  printf("makehisto: Use %i bins in [%g:%g]... \n",NBIN,min,max);

  // create histogram
  double *xbord=calloc(NBIN+1,sizeof(double));

  for (i=0;i<NBIN;i++) { 
    x[i]=min+(i+.5)*width; 
    xbord[i]=min+i*width;
  }
  xbord[NBIN]=max;

  ex=0;
  for (j=0;j<size;j++) {
    if (sample[j]<min) ex++;
    else if (sample[j]>max) ex++;
    else {
      bisect_makehisto(xbord,NBIN+1,sample[j],&ix);
      hist[ix]+=1.0;
    }
  }
  printf("makehisto: %i samples not in range [%g:%g]... \n",ex,min,max);

  // normalise histogram
  for (i=0;i<NBIN;i++) {
    hist[i]/=1.*(size-ex);
  }

  free(xbord);
  return;
}


void makehisto_multidim(double **sample,int dimension,int size,double **x,double *hist,int *nhistobin,double *minfix,double *maxfix,int *exclude)
{
  int i,j,ix[dimension],NBIN[dimension],ex,nbinceil,flag,dimtot,dimpart,index=0;
  double min[dimension],max[dimension],scale[dimension],width[dimension];

  // define binning
  for (j=0;j<dimension;j++) {
    min[j]=1.e10;
    max[j]=-1.e10;
    if (nhistobin[j]==0) {         // automatic binning
      for (i=0;i<size;i++) {
	if (min[j]>sample[j][i]) min[j]=sample[j][i];
	if (max[j]<sample[j][i]) max[j]=sample[j][i];
      }
      scale[j]=pow(10,floor(log10(fabs(max[j]-min[j]))));
      min[j]=scale[j]*floor(min[j]/scale[j]);
      max[j]=scale[j]*ceil(max[j]/scale[j]);
      NBIN[j]=(int)floor(pow(size,1./(2.*dimension)));
      nhistobin[j]=NBIN[j];
    }
    else {
      NBIN[j]=nhistobin[j];
      min[j]=minfix[j];
      max[j]=maxfix[j];
    }
    width[j]=(max[j]-min[j])/(1.*NBIN[j]);
    printf("makehisto: Use %i bins in [%g:%g] in dimension %i\n",NBIN[j],min[j],max[j],j+1);
  }

  // create histogram
  nbinceil=1;
  dimtot=1;
  for (j=0;j<dimension;j++) {
    if (NBIN[j]>nbinceil) nbinceil=NBIN[j];
    dimtot*=NBIN[j];
  }
  double **xbord=calloc(dimension,sizeof(double *));
  for (j=0;j<dimension;j++) {
    xbord[j]=calloc(nbinceil+1,sizeof(double));
  }

  for (j=0;j<dimension;j++) {
    for (i=0;i<NBIN[j];i++) { 
      x[j][i]=min[j]+(i+.5)*width[j]; 
      xbord[j][i]=min[j]+i*width[j];
    }
    xbord[j][NBIN[j]]=max[j];
  }

  ex=0;
  for (i=0;i<size;i++) {
    flag=0;
    for (j=0;j<dimension;j++) {
      if ((sample[j][i]<min[j])||(sample[j][i]>max[j])) flag=1;
    }
    if (!flag) {
      index=0;
      dimpart=dimtot;
      for (j=0;j<dimension;j++) {
	bisect_makehisto(xbord[j],NBIN[j]+1,sample[j][i],&ix[j]);
	dimpart/=NBIN[j];
	index+=ix[j]*dimpart;
      }
      hist[index]+=1.0;
    }
    else ex++;
  }
  printf("makehisto: %i samples not in histogram range ... \n",ex);

  // normalise histogram
  for (i=0;i<dimtot;i++) {
    hist[i]/=1.*(size-ex);
  }
  *exclude=ex;

  for (j=0;j<dimension;j++) {
    free(xbord[j]);
  }
  free(xbord);
  return;
}



void sample_histogram(gsl_histogram *histo,int nsample,double *samples)
{
  int i;
  double uran;

  const gsl_rng_type *T=gsl_rng_default;
  gsl_rng *r=gsl_rng_alloc(T);
  gsl_rng_env_setup();

  gsl_histogram_pdf *phisto=gsl_histogram_pdf_alloc(histo->n);
  gsl_histogram_pdf_init(phisto,histo);

  for (i=0;i<nsample;i++) {
    uran=gsl_rng_uniform(r);
    samples[i]=gsl_histogram_pdf_sample(phisto,uran);
  }
  gsl_histogram_pdf_free(phisto);
  gsl_rng_free(r);
  return;
}


void sample_histogram_2D(gsl_histogram2d *histo,int nsample,double **samples)
{
  int i;
  double uran1,uran2;

  const gsl_rng_type *T=gsl_rng_default;
  gsl_rng *r=gsl_rng_alloc(T);
  gsl_rng_env_setup();

  gsl_histogram2d_pdf *phisto=gsl_histogram2d_pdf_alloc(histo->nx,histo->ny);
  gsl_histogram2d_pdf_init(phisto,histo);

  for (i=0;i<nsample;i++) {
    uran1=gsl_rng_uniform(r);
    uran2=gsl_rng_uniform(r);
    gsl_histogram2d_pdf_sample(phisto,uran1,uran2,&samples[0][i],&samples[1][i]);
  }
  gsl_histogram2d_pdf_free(phisto);
  gsl_rng_free(r);
  return;
}

void sample_histogram_3D(gsl_histogram2d **histo,double *r1,int nbin1,int nsample,double **samples)
{
  int i,j,k;
  double uran1,uran2,uran3,norm=0.0;
  double *marg=calloc(nbin1,sizeof(double));
  double *margbound=calloc(nbin1+1,sizeof(double));
  int nbin2=histo[0]->nx;
  int nbin3=histo[0]->ny;

  // Create marginal cdf
  for (i=0;i<nbin1;i++) {
    for (j=0;j<nbin2;j++) {
      for (k=0;k<nbin3;k++) {
	marg[i]+=gsl_histogram2d_get(histo[i],j,k);
      }
    }
    norm+=marg[i];
  }

  margbound[0]=0.0;
  for (i=0;i<nbin1;i++) {
    marg[i]/=norm;  // normalise to unity
    margbound[i+1]=margbound[i]+marg[i];
  }

  // Sample from histogram
  const gsl_rng_type *T=gsl_rng_default;
  gsl_rng *r=gsl_rng_alloc(T);
  gsl_rng_env_setup();

  gsl_histogram2d_pdf **phisto=calloc(nbin1,sizeof(gsl_histogram2d_pdf *));
  for (i=0;i<nbin1;i++) {
    phisto[i]=gsl_histogram2d_pdf_alloc(nbin2,nbin3);
    gsl_histogram2d_pdf_init(phisto[i],histo[i]);
  }

  for (k=0;k<nsample;k++) {
    uran1=gsl_rng_uniform(r);
    for (i=0;i<nbin1;i++) {
      if ((uran1>=margbound[i])&&(uran1<margbound[i+1])) break;
    }
    samples[0][k]=r1[i]+(uran1-margbound[i])/(margbound[i+1]-margbound[i])*(r1[i+1]-r1[i]);
    uran2=gsl_rng_uniform(r);
    uran3=gsl_rng_uniform(r);
    gsl_histogram2d_pdf_sample(phisto[i],uran2,uran3,&samples[1][k],&samples[2][k]);
  }

  for (i=0;i<nbin1;i++) {
    gsl_histogram2d_pdf_free(phisto[i]);
  }
  free(phisto);
  gsl_rng_free(r);
  free(marg);
  free(margbound);
  return;
}

void read_histogram_3D(FILE *dat,gsl_histogram2d **histo,double *r1,int nbin1,int nbin2,int nbin3)
{
  int i,j,k;
  double *r2=calloc(nbin2+1,sizeof(double));
  double *r3=calloc(nbin3+1,sizeof(double));
  double ***w=calloc(nbin1+1,sizeof(double *));
  for (i=0;i<nbin1+1;i++) {
    w[i]=calloc(nbin2+1,sizeof(double *));
    for (j=0;j<nbin2+1;j++) {
      w[i][j]=calloc(nbin3+1,sizeof(double));
    }
  }

  for (i=0;i<nbin1;i++) {
    for (j=0;j<nbin2;j++) {
      for (k=0;k<nbin3;k++) {
	fscanf(dat,"%lf  %lf",&r1[i],&r1[i+1]);
	fscanf(dat,"%lf  %lf",&r2[j],&r2[j+1]);
	fscanf(dat,"%lf  %lf",&r3[k],&r3[k+1]);
	fscanf(dat,"%lf",&w[i][j][k]);
      }
    }
  }
 
  for (i=0;i<nbin1;i++) {
    gsl_histogram2d_set_ranges(histo[i],r2,nbin2+1,r3,nbin3+1);
    for (j=0;j<nbin2;j++) {
      for (k=0;k<nbin3;k++) {
	gsl_histogram2d_accumulate(histo[i],0.5*(r2[j]+r2[j+1]),0.5*(r3[k]+r3[k+1]),w[i][j][k]);
      }
    }
  }

  for (i=0;i<nbin1+1;i++) {
    for (j=0;j<nbin2+1;j++) {
      free(w[i][j]);
    }
    free(w[i]);
  }
  free(w);
  free(r2);
  free(r3);
  return;
}



// to get_confidence_intervals
int compare_descend(double *a,double *b)
{
  if (*a > *b ) return -1;
  else if (*a < *b ) return 1;
  else return 0;
}


void get_confidence_intervals(double *x,double *hist,int nbin,double *levels,int nlevel,double **res)  //res[nlevel][2]
{
  int i,j,maxindex=0,lowindex=0.0,higindex=0.0;
  double sum=0.0,max=-1.e10,fraclow=0.0,frachig=0.0;
  double *hist_sorted=calloc(nbin,sizeof(double));
  int flag[nlevel];
  double bound[nlevel];

  // normalise histogram & find maximum
  for (i=0;i<nbin;i++) {
    sum+=hist[i];
  }
  for (i=0;i<nbin;i++) {
    hist[i]/=sum;
    hist_sorted[i]=hist[i];  //duplicate histogram
    if (hist[i]>max) {
      max=hist[i];
      maxindex=i;
    }
  }

  // find levels
  gsl_heapsort(hist_sorted,nbin,sizeof(double),(gsl_comparison_fn_t)compare_descend);

  sum=0.0;
  for (j=0;j<nlevel;j++) { 
    flag[j]=0;
  }

  for (i=0;i<nbin;i++) { 
    sum+=hist_sorted[i];
    for (j=0;j<nlevel;j++) { 
      if ((sum >= levels[j]) && (!flag[j])) {
	bound[j]=hist_sorted[i];
	flag[j]=1;
      }
    }
  }

  // get intervals
  for (j=0;j<nlevel;j++) {
    for (i=maxindex;i>=0;i--) { 
      if (hist[i]<bound[j]) {
	lowindex=i;
	fraclow=(bound[j]-hist[i])/(hist[i+1]-hist[i]);
	break;
      }
    }
    res[j][0]=x[maxindex]-(x[lowindex]+fraclow*(x[lowindex+1]-x[lowindex]));
    for (i=maxindex;i<nbin;i++) { 
      if (hist[i]<bound[j]) {
	higindex=i-1;
	frachig=(bound[j]-hist[i-1])/(hist[i]-hist[i-1]);
	break;
      }
    }
    res[j][1]=(x[higindex]+frachig*(x[higindex+1]-x[higindex]))-x[maxindex];
  }

  // clean up
  free(hist_sorted);
  return;
}


// get confidence levels
void get_confidence_levels(double *hist,int size,double *levels,int nlevel,double *res)  //res[nlevel]
{
  int i,j;
  double sum=0.0;
  double *hist_sorted=calloc(size,sizeof(double));
  int flag[nlevel];
  for (j=0;j<nlevel;j++) { 
    flag[j]=0;
  }

  // normalise
  for (i=0;i<size;i++) { 
    sum+=hist[i];
  }
  for (i=0;i<size;i++) { 
    hist[i]/=sum;
    hist_sorted[i]=hist[i];  //duplicate histogram
  }

  // find levels
  gsl_heapsort(hist_sorted,size,sizeof(double),(gsl_comparison_fn_t)compare_descend);
  sum=0.0;
  for (i=0;i<size;i++) {                
    sum+=hist_sorted[i];
    for (j=0;j<nlevel;j++) { 
      if ((sum >= levels[j]) && (!flag[j])) {
	res[j]=hist_sorted[i];
	flag[j]=1;
      }
    }
  }

  // clean up
  free(hist_sorted);
  return;
}


// compute all marginal distributions of 3D gridded probability
void marginalise_grid_3D(bj_prob_grid_3D *p)
{
  int i,j,k;
  double sum,check;

  // normalise 3D probability
  sum=0.0;
  for (i=0;i<p->nx*p->ny*p->nz;i++) {  
    sum+=p->prob[i];
  }
  for (i=0;i<p->nx*p->ny*p->nz;i++) {  
    p->prob[i]/=sum;
  }

  // 2D marginalisation - xy
  check=0.0;
  for (i=0;i<p->nx;i++) {                
    for (j=0;j<p->ny;j++) {
      sum=0.0;
      for (k=0;k<p->nz;k++) {                
	sum+=p->prob[k+p->nz*(j+p->ny*i)];
      }
      p->xy[j+p->ny*i]=sum;
      check+=sum;
    }
  }
  if (fabs(check-1.) > 1.e-5) printf("marginalisation (2D) for x-y incorrect: %g\n",check);

  // 2D marginalisation - xz
  check=0.0;
  for (i=0;i<p->nx;i++) {                
    for (j=0;j<p->nz;j++) {
      sum=0.0;
      for (k=0;k<p->ny;k++) {                
	sum+=p->prob[j+p->nz*(k+p->ny*i)];
      }
      p->xz[j+p->nz*i]=sum;
      check+=sum;
    }
  }
  if (fabs(check-1.) > 1.e-5) printf("marginalisation (2D) for x-z incorrect: %g\n",check);

  // 2D marginalisation - yz
  check=0.0;
  for (i=0;i<p->ny;i++) {                
    for (j=0;j<p->nz;j++) {
      sum=0.0;
      for (k=0;k<p->nx;k++) {                
	sum+=p->prob[j+p->nz*(i+p->ny*k)];
      }
      p->yz[j+p->nz*i]=sum;
      check+=sum;
    }
  }
  if (fabs(check-1.) > 1.e-5) printf("marginalisation (2D) for y-z incorrect: %g\n",check);

  // 1D marginalisation - x
  check=0.0;
  for (i=0;i<p->nx;i++) {  
    sum=0.0;
    for (j=0;j<p->ny;j++) {
      sum+=p->xy[j+p->ny*i];
    }
    p->x[i]=sum;
    check+=sum;
  }
  if (fabs(check-1.) > 1.e-5) printf("marginalisation (1D) for x incorrect: %g\n",check);

  // 1D marginalisation - y
  check=0.0;
  for (i=0;i<p->ny;i++) {  
    sum=0.0;
    for (j=0;j<p->nx;j++) {
      sum+=p->xy[i+p->ny*j];
    }
    p->y[i]=sum;
    check+=sum;
  }
  if (fabs(check-1.) > 1.e-5) printf("marginalisation (1D) for y incorrect: %g\n",check);

  // 1D marginalisation - z
  check=0.0;
  for (i=0;i<p->nz;i++) {  
    sum=0.0;
    for (j=0;j<p->ny;j++) {
      sum+=p->yz[i+p->nz*j];
    }
    p->z[i]=sum;
    check+=sum;
  }
  if (fabs(check-1.) > 1.e-5) printf("marginalisation (1D) for z incorrect: %g\n",check);

  return;
}


// get bootstrap resampling of data vector (same size)
// orig/resampled: size1 -> dimension of data vector; size2 -> # samples
void get_bootstrap_resampling(gsl_rng *rnd,gsl_matrix *orig,gsl_matrix *resampled)
{
  int i;
  const int nr=(int)orig->size2;
  int *indices=calloc(nr,sizeof(int));
  int *sequence=calloc(nr,sizeof(int));
  gsl_vector_view column;

  for (i=0;i<nr;i++) {
   sequence[i]=i;
  }
  gsl_ran_sample(rnd,indices,nr,sequence,nr,sizeof(int));

  for (i=0;i<nr;i++) {
    column=gsl_matrix_column(orig,indices[i]);      // extracts column indices[i] of orig
    gsl_matrix_set_col(resampled,i,&column.vector); // adds columns to resampled
  }

  free(indices);
  free(sequence);
  return;
}

