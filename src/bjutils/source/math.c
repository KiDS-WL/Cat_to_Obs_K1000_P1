// mathematical routines

#include "../include/bjutils.h"


double integrate_qag(double (*func)(double, void*),void *par,double low,double up,double *error)
{
  const int lim=2000;
  const int key=5;
  const double eabs=0.0;
  const double erel=1.e-4;
  double res,err;

  gsl_integration_workspace *work=gsl_integration_workspace_alloc(lim);

  gsl_function F;
  F.function=func;
  F.params=par;

  gsl_integration_qag(&F,low,up,eabs,erel,lim,key,work,&res,&err);

  gsl_integration_workspace_free(work);
  if (error!=NULL) *error=err; 
  return(res);
}


double findroot(double (*func)(double, void*),void *par,double low,double hig,double tol,double *error)
{
  #include <gsl/gsl_errno.h>
  int status,iter=0;
  double root,rnext;
  const int max_iter=1000;    // maximum no. of iterations

  gsl_function F;
  F.function=func;
  F.params=par;
     
  const gsl_root_fsolver_type *T=gsl_root_fsolver_brent;
  gsl_root_fsolver *solv=gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(solv,&F,low,hig);
     
  do {
    iter++;
    status=gsl_root_fsolver_iterate(solv);
    root=gsl_root_fsolver_root(solv);
    low=gsl_root_fsolver_x_lower(solv);
    hig=gsl_root_fsolver_x_upper(solv);
    status=gsl_root_test_interval(low,hig,0.0,tol);
    if (status==GSL_SUCCESS) printf ("Root finder converged after %i iterations\n",iter);  
  } while (status==GSL_CONTINUE && iter<max_iter);
     
  gsl_root_fsolver_iterate(solv);
  rnext=gsl_root_fsolver_root(solv);
  *error=fabs(root-rnext);
  gsl_root_fsolver_free(solv);
  return(root);
}





// index search via bisection (log N process)
// array, dimension of array, value to be found, returned: index of array lower than x
void bisect_interpol(double *a,int dim,double x,int *erg)
{
  int l=0,u=dim-1,m;

  if ((x<a[l]) || (x>a[u])) {
    printf("Error: index out of range: %g not in %g-%g\n",x,a[l],a[u]);
    exit(-1);
  }
  while (u-l>1) {
    m=(u+l)/2;
    if (x>a[m]) l=m;
    else u=m;
  }
  *erg=l;
}

// 1D linear interpolation + extrapolation
double interpol_linear_extra_1D(double *ax,int NX,double *z,double x)
{
  if (NX==1) return(z[0]);  // assume z constant if only one data point
  if (x<=ax[0]) return((z[1]-z[0])/(ax[1]-ax[0])*(x-ax[0])+z[0]);
  else if (x>=ax[NX-1]) return((z[NX-1]-z[NX-2])/(ax[NX-1]-ax[NX-2])*(x-ax[NX-1])+z[NX-1]);
  else {
    int ix;
    bisect_interpol(ax,NX,x,&ix);
    double t=(x-ax[ix])/(ax[ix+1]-ax[ix]);
    return((1-t)*z[ix]+t*z[ix+1]);
  }
}


// 2D linear interpolation + extrapolation
double interpol_linear_extra_2D(double *ax,int NX,double *ay,int NY,double **z,double x,double y)
{
  if ((x<=ax[0])&&(y>ay[0])&&(y<ay[NY-1])) {
    double int1=interpol_linear_extra_1D(ay,NY,z[0],y);
    double int2=interpol_linear_extra_1D(ay,NY,z[1],y);
    return(int1+(int2-int1)/(ax[1]-ax[0])*(x-ax[0]));
  }
  else if ((x>=ax[NX-1])&&(y>ay[0])&&(y<ay[NY-1])) {
    double int1=interpol_linear_extra_1D(ay,NY,z[NX-2],y);
    double int2=interpol_linear_extra_1D(ay,NY,z[NX-1],y);
    return(int1+(int2-int1)/(ax[NX-1]-ax[NX-2])*(x-ax[NX-2]));
  }
  else if ((y<=ay[0])&&(x>ax[0])&&(x<ax[NX-1])) {
    int i;
    double tmp1[NX],tmp2[NX];
    for(i=0;i<NX;i++) {
      tmp1[i]=z[i][0];
      tmp2[i]=z[i][1];
    }
    double int1=interpol_linear_extra_1D(ax,NX,tmp1,x);
    double int2=interpol_linear_extra_1D(ax,NX,tmp2,x);
    return(int1+(int2-int1)/(ay[1]-ay[0])*(y-ay[0]));
  }
  else if ((y>=ay[NY-1])&&(x>ax[0])&&(x<ax[NX-1])) {
    int i;
    double tmp1[NX],tmp2[NX];
    for(i=0;i<NX;i++) {
      tmp1[i]=z[i][NY-2];
      tmp2[i]=z[i][NY-1];
    }
    double int1=interpol_linear_extra_1D(ax,NX,tmp1,x);
    double int2=interpol_linear_extra_1D(ax,NX,tmp2,x);
    return(int1+(int2-int1)/(ay[NY-1]-ay[NY-2])*(y-ay[NY-2]));
  }
  else {  // produces error via 'bisect' if both dimensions outside range
    int ix,iy;
    bisect_interpol(ax,NX,x,&ix);
    bisect_interpol(ay,NY,y,&iy);
    double t=(x-ax[ix])/(ax[ix+1]-ax[ix]);
    double u=(y-ay[iy])/(ay[iy+1]-ay[iy]);
    return((1-t)*(1-u)*z[ix][iy]+t*(1-u)*z[ix+1][iy]+(1-t)*u*z[ix][iy+1]+t*u*z[ix+1][iy+1]);
  }
}


// 3D linear interpolation + extrapolation
double interpol_linear_extra_3D(double *ax,int NX,double *ay,int NY,double *az,int NZ,double ***f,double x,double y,double z)
{
  if ((x<=ax[0])&&(y>ay[0])&&(y<ay[NY-1])&&(z>az[0])&&(z<az[NZ-1])) {   // x out of range, low
    double int1=interpol_linear_extra_2D(ay,NY,az,NZ,f[0],y,z);
    double int2=interpol_linear_extra_2D(ay,NY,az,NZ,f[1],y,z);
    return(int1+(int2-int1)/(ax[1]-ax[0])*(x-ax[0]));
  }
  else if ((x>=ax[NX-1])&&(y>ay[0])&&(y<ay[NY-1])&&(z>az[0])&&(z<az[NZ-1])) {   // x out of range, high
    double int1=interpol_linear_extra_2D(ay,NY,az,NZ,f[NX-2],y,z);
    double int2=interpol_linear_extra_2D(ay,NY,az,NZ,f[NX-1],y,z);
    return(int1+(int2-int1)/(ax[NX-1]-ax[NX-2])*(x-ax[NX-2]));
  }
  else if ((y<=ay[0])&&(x>ax[0])&&(x<ax[NX-1])&&(z>az[0])&&(z<az[NZ-1])) {   // y out of range, low
    int i,j;
    double **tmp1=malloc(NX*sizeof(double *));
    double **tmp2=malloc(NX*sizeof(double *));
    for(i=0;i<NX;i++) {
      tmp1[i]=malloc(NZ*sizeof(double));
      tmp2[i]=malloc(NZ*sizeof(double));
    }
    for(i=0;i<NX;i++) {
      for(j=0;j<NZ;j++) {
	tmp1[i][j]=f[i][0][j];
	tmp2[i][j]=f[i][1][j];
      }
    }
    double int1=interpol_linear_extra_2D(ax,NX,az,NZ,tmp1,x,z);
    double int2=interpol_linear_extra_2D(ax,NX,az,NZ,tmp2,x,z);
    for(i=0;i<NX;i++) {
      free(tmp1[i]);
      free(tmp2[i]);
    }
    free(tmp1);
    free(tmp2);
    return(int1+(int2-int1)/(ay[1]-ay[0])*(y-ay[0]));
  }
  else if ((y>=ay[NY-1])&&(x>ax[0])&&(x<ax[NX-1])&&(z>az[0])&&(z<az[NZ-1])) {   // y out of range, high
    int i,j;
    double **tmp1=malloc(NX*sizeof(double *));
    double **tmp2=malloc(NX*sizeof(double *));
    for(i=0;i<NX;i++) {
      tmp1[i]=malloc(NZ*sizeof(double));
      tmp2[i]=malloc(NZ*sizeof(double));
    }
    for(i=0;i<NX;i++) {
      for(j=0;j<NZ;j++) {
	tmp1[i][j]=f[i][NY-2][j];
	tmp2[i][j]=f[i][NY-1][j];
      }
    }
    double int1=interpol_linear_extra_2D(ax,NX,az,NZ,tmp1,x,z);
    double int2=interpol_linear_extra_2D(ax,NX,az,NZ,tmp2,x,z);
    for(i=0;i<NX;i++) {
      free(tmp1[i]);
      free(tmp2[i]);
    }
    free(tmp1);
    free(tmp2);
    return(int1+(int2-int1)/(ay[NY-1]-ay[NY-2])*(y-ay[NY-2]));
  }
  else if ((z<=az[0])&&(x>ax[0])&&(x<ax[NX-1])&&(y>ay[0])&&(y<ay[NY-1])) {   // z out of range, low
    int i,j;
    double **tmp1=malloc(NX*sizeof(double *));
    double **tmp2=malloc(NX*sizeof(double *));
    for(i=0;i<NX;i++) {
      tmp1[i]=malloc(NY*sizeof(double));
      tmp2[i]=malloc(NY*sizeof(double));
    }
    for(i=0;i<NX;i++) {
      for(j=0;j<NY;j++) {
	tmp1[i][j]=f[i][j][0];
	tmp2[i][j]=f[i][j][1];
      }
    }
    double int1=interpol_linear_extra_2D(ax,NX,ay,NY,tmp1,x,y);
    double int2=interpol_linear_extra_2D(ax,NX,ay,NY,tmp2,x,y);
    for(i=0;i<NX;i++) {
      free(tmp1[i]);
      free(tmp2[i]);
    }
    free(tmp1);
    free(tmp2);
    return(int1+(int2-int1)/(az[1]-az[0])*(z-az[0]));
  }
  else if ((z>=az[NZ-1])&&(x>ax[0])&&(x<ax[NX-1])&&(y>ay[0])&&(y<ay[NY-1])) {   // z out of range, high
    int i,j;
    double **tmp1=malloc(NX*sizeof(double *));
    double **tmp2=malloc(NX*sizeof(double *));
    for(i=0;i<NX;i++) {
      tmp1[i]=malloc(NY*sizeof(double));
      tmp2[i]=malloc(NY*sizeof(double));
    }
    for(i=0;i<NX;i++) {
      for(j=0;j<NY;j++) {
	tmp1[i][j]=f[i][j][NZ-2];
	tmp2[i][j]=f[i][j][NZ-1];
      }
    }
    double int1=interpol_linear_extra_2D(ax,NX,ay,NY,tmp1,x,y);
    double int2=interpol_linear_extra_2D(ax,NX,ay,NY,tmp2,x,y);
    for(i=0;i<NX;i++) {
      free(tmp1[i]);
      free(tmp2[i]);
    }
    free(tmp1);
    free(tmp2);
    return(int1+(int2-int1)/(az[NZ-1]-az[NZ-2])*(z-az[NZ-2]));
  }
  else {  // produces error via 'bisect' if more than 1 dimension outside range
    int ix,iy,iz;
    double f1a,f1b,f2a,f2b;

    bisect_interpol(ax,NX,x,&ix);
    bisect_interpol(ay,NY,y,&iy);
    bisect_interpol(az,NZ,z,&iz);

    double q1=(x-ax[ix])/(ax[ix+1]-ax[ix]);
    double q2=(y-ay[iy])/(ay[iy+1]-ay[iy]);
    double q3=(z-az[iz])/(az[iz+1]-az[iz]);
    double p1=1-q1;
    double p2=1-q2;

    f1a=p1*f[ix][iy][iz]+q1*f[ix+1][iy][iz];
    f1b=p1*f[ix][iy+1][iz]+q1*f[ix+1][iy+1][iz];
    f2a=p2*f1a+q2*f1b;

    f1a=p1*f[ix][iy][iz+1]+q1*f[ix+1][iy][iz+1];
    f1b=p1*f[ix][iy+1][iz+1]+q1*f[ix+1][iy+1][iz+1];
    f2b=p2*f1a+q2*f1b;

    return((1-q3)*f2a+q3*f2b);
  }
}



// 4D linear interpolation + extrapolation
double interpol_linear_extra_4D(double *ax,int NX,double *ay,int NY,double *az,int NZ,double *aw,int NW,double ****f,double x,double y,double z,double w)
{
  int ix;
  double f1,f2,t;

  if (((x<ax[0])||(x>ax[NX-1]))&&((y<ay[0])||(y>ay[NY-1])||(z<az[0])||(z>az[NZ-1])||(w<aw[0])||(w>aw[NW-1]))) {
    printf("Error in routine 'interpol_linear_extra_4D': attempt extrapolation in more than one dimension.\n");
    exit(-1);
  }

  if (x<=ax[0]) {
    f1=interpol_linear_extra_3D(ay,NY,az,NZ,aw,NW,f[0],y,z,w);
    f2=interpol_linear_extra_3D(ay,NY,az,NZ,aw,NW,f[1],y,z,w);
    return((f2-f1)/(ax[1]-ax[0])*(x-ax[0])+f1);
  }
  else if (x>=ax[NX-1]) {
    f1=interpol_linear_extra_3D(ay,NY,az,NZ,aw,NW,f[NX-2],y,z,w);
    f2=interpol_linear_extra_3D(ay,NY,az,NZ,aw,NW,f[NX-1],y,z,w);
    return((f2-f1)/(ax[NX-1]-ax[NX-2])*(x-ax[NX-1])+f2);
  }
  else {
    bisect_interpol(ax,NX,x,&ix);

    f1=interpol_linear_extra_3D(ay,NY,az,NZ,aw,NW,f[ix],y,z,w);
    f2=interpol_linear_extra_3D(ay,NY,az,NZ,aw,NW,f[ix+1],y,z,w);

    t=(x-ax[ix])/(ax[ix+1]-ax[ix]);
    return((1-t)*f1+t*f2);
  }
}




void mattimesmat_gen(gsl_matrix *a,gsl_matrix *b,gsl_matrix *erg,int dim1,int dim2,int dim3)
{
  int i,j,k;
  double sum;
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim3;j++) {
      sum=0.0;
      for (k=0;k<dim2;k++) {
        sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(erg,i,j,sum);
    }
  }
}



void mattimesmat(gsl_matrix *a,gsl_matrix *b,gsl_matrix *erg,int N)
{
  int i,j,k;
  double sum;
  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      sum=0.0;
      for (k=0;k<N;k++) {
        sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(erg,i,j,sum);
    }
  }
}



void mattimesvec_gen(gsl_matrix *a,gsl_vector *b,gsl_vector *erg,int dim1,int dim2)
{
  int i,j;
  double sum;
  for (i=0;i<dim1;i++) {
    sum=0.0;
    for (j=0;j<dim2;j++) {
      sum+=gsl_matrix_get(a,i,j)*gsl_vector_get(b,j);
    }
    gsl_vector_set(erg,i,sum);
  }
}

void mattimesvec(gsl_matrix *a,gsl_vector *b,gsl_vector *erg,int N)
{
  int i,j;
  double sum;
  for (i=0;i<N;i++) {
    sum=0.0;
    for (j=0;j<N;j++) {
      sum+=gsl_matrix_get(a,i,j)*gsl_vector_get(b,j);
    }
    gsl_vector_set(erg,i,sum);
  }
}


int matrixinvers(gsl_matrix *m,gsl_matrix *inv,int DIM)
{
  int i,j,k,flag=0;
  double elem,sum;
  const int DOCHECK=1;

  gsl_matrix *v=gsl_matrix_alloc(DIM,DIM);
  gsl_matrix *id=gsl_matrix_alloc(DIM,DIM);
  gsl_matrix *u=gsl_matrix_alloc(DIM,DIM);
  gsl_vector *work=gsl_vector_alloc(DIM);
  gsl_vector *s=gsl_vector_alloc(DIM);


  gsl_matrix_memcpy(u,m);
  gsl_linalg_SV_decomp(u,v,s,work);

  for (i=0;i<DIM;i++) {
    elem=gsl_vector_get(s,i);
    if (fabs(elem)==0) {
      printf("Error: zero singular value!\n");
      exit(-1);
    }
    if (fabs(elem)<1.e-35) {
      printf("Possible singular value: %g\n",elem);
      flag=1;
    }
    //printf("eigenvalue: %g\n",elem);
  }

  for (i=0;i<DIM;i++) {
    for (j=0;j<DIM;j++) {
      sum=0.0;
      for (k=0;k<DIM;k++) {
	sum+=gsl_matrix_get(v,i,k)/gsl_vector_get(s,k)*gsl_matrix_get(u,j,k);
      }
      gsl_matrix_set(inv,i,j,sum);
    }
  }

  if (DOCHECK) {
    mattimesmat(m,inv,v,DIM);           //v reused
    gsl_matrix_set_identity(id);
    gsl_matrix_sub(v,id);                //result in v
    for (i=0;i<DIM;i++) {
      for (j=0;j<DIM;j++) {
	elem=gsl_matrix_get(v,i,j);
	if (fabs(elem)>1.e-2) {
	  printf("Problematic matrix inversion: %g\n",elem);
	  flag=2;
	  break;
	}
      }
    }
  }

  gsl_matrix_free(u);
  gsl_matrix_free(v);
  gsl_matrix_free(id);
  gsl_vector_free(work);
  gsl_vector_free(s);
  return(flag);
}


gsl_matrix *pseudoinvers(gsl_matrix *m,double thresh,double *eigen)
{
  int i,j,k;
  double sum;

  const int dim=(int)m->size1;
  if (m->size1!=m->size2) {
    printf("Error in 'pseudoinvers': non-square matrix - not implemented yet.\n");
    exit(-1);
  }

  gsl_matrix *res=gsl_matrix_alloc(dim,dim);
  gsl_matrix *v=gsl_matrix_alloc(dim,dim);
  gsl_matrix *cov=gsl_matrix_alloc(dim,dim);
  gsl_vector *work=gsl_vector_alloc(dim);
  gsl_vector *s=gsl_vector_alloc(dim);
  gsl_vector *invs=gsl_vector_alloc(dim);

  gsl_matrix_memcpy(cov,m);
  gsl_linalg_SV_decomp(cov,v,s,work);

  for (i=0;i<dim;i++) {
    eigen[i]=gsl_vector_get(s,i);
    if (eigen[i]<=thresh) gsl_vector_set(invs,i,0.0);
    else gsl_vector_set(invs,i,1./eigen[i]);
  }

  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      sum=0.0;
      for (k=0;k<dim;k++) {
	sum+=gsl_matrix_get(v,i,k)*gsl_vector_get(invs,k)*gsl_matrix_get(cov,j,k);
      }
      gsl_matrix_set(res,i,j,sum);
    }
  }

  gsl_matrix_free(cov);
  gsl_matrix_free(v);
  gsl_vector_free(work);
  gsl_vector_free(s);
  gsl_vector_free(invs);
  return(res);
}


int covinvers(gsl_matrix *m,gsl_matrix *inv,int DIM)
{
  #include<gsl/gsl_errno.h>
  gsl_set_error_handler_off();
  const int DOCHECK=1;

  int i,flag=0;
  gsl_matrix *check=gsl_matrix_alloc(DIM,DIM);
  gsl_vector *b=gsl_vector_alloc(DIM);

  gsl_matrix_memcpy(check,m);

  flag=gsl_linalg_cholesky_decomp(check);

  if (!flag) {
    for (i=0;i<DIM;i++) {
      gsl_vector_set_basis(b,i);
      gsl_linalg_cholesky_svx(check,b);
      gsl_matrix_set_col(inv,i,b);
    }

    if (DOCHECK) {
      gsl_matrix *id=gsl_matrix_alloc(DIM,DIM);
      gsl_matrix_set_identity(id);
      int j;
      double elem;

      mattimesmat(m,inv,check,DIM);           //check reused
      gsl_matrix_sub(check,id);               //result in check
      for (i=0;i<DIM;i++) {
	for (j=0;j<DIM;j++) {
	  elem=gsl_matrix_get(check,i,j);
	  if (fabs(elem)>1.e-6) {
	    printf("Problematic matrix inversion: %g\n",elem);
	    flag=2;
	    break;
	  }
	}
      }
      gsl_matrix_free(id);
    }
  }

  gsl_matrix_free(check);
  gsl_vector_free(b);
  return(flag);
}

double matrix_trace(gsl_matrix *a,int N)
{
  int i;
  double sum=0.0;
  for (i=0;i<N;i++) {
    sum+=gsl_matrix_get(a,i,i);
  }
  return(sum);
}

double matrix_det(gsl_matrix *a,int N)
{
  int sign;
  gsl_permutation *perm=gsl_permutation_alloc(N);
  gsl_matrix *LU=gsl_matrix_alloc(N,N);
  gsl_matrix_memcpy(LU,a);
  gsl_linalg_LU_decomp(LU,perm,&sign);
  gsl_permutation_free(perm);
  return(gsl_linalg_LU_det(LU,sign));
}


void matrix_transpose(gsl_matrix *a,gsl_matrix *t,int dim1,int dim2)
{
  int i,j;
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      gsl_matrix_set(t,j,i,gsl_matrix_get(a,i,j));
    }
  }
  return;
}

void matrix_print(gsl_matrix *a,char *file)
{
  int i,j;
  FILE *dat;
  if ((dat=fopen(file,"w"))==NULL) {
    printf("Error in routine 'matrix_print': could not open file %s\n",file);
    exit(-1);
  }
  for (i=0;i<a->size1;i++) {
    for (j=0;j<a->size2;j++) {
      fprintf(dat,"%15.10g\t",gsl_matrix_get(a,i,j));
    }
    fprintf(dat,"\n");
  }
  fclose(dat);
  return;
}


void vector_print(gsl_vector *v,char *file)
{
  int i;
  FILE *dat;
  if ((dat=fopen(file,"w"))==NULL) {
    printf("Error in routine 'vector_print': could not open file %s\n",file);
    exit(-1);
  }
  for (i=0;i<v->size;i++) {
      fprintf(dat,"%15.10g\n",gsl_vector_get(v,i));
  }
  fclose(dat);
  return;
}


void solve_lse(gsl_matrix *m,gsl_vector *in,gsl_vector *res)
{
  if (m->size1!=res->size) {
    printf("Error in routine 'solve_lse': matrix row dimension does not match dimension of solution vector: %zu - %zu\n",m->size1,res->size);
    exit(-1);
  }
  if (m->size2!=in->size) {
    printf("Error in routine 'solve_lse': matrix column dimension does not match dimension of input vector: %zu - %zu\n",m->size2,in->size);
    exit(-1);
  }

  gsl_matrix *v=gsl_matrix_alloc(m->size2,m->size2);
  gsl_matrix *u=gsl_matrix_alloc(m->size1,m->size2);
  gsl_vector *work=gsl_vector_alloc(m->size2);
  gsl_vector *s=gsl_vector_alloc(m->size2);

  gsl_matrix_memcpy(u,m);
  gsl_linalg_SV_decomp(u,v,s,work);
  gsl_linalg_SV_solve(u,v,s,in,res);

  gsl_matrix_free(u);
  gsl_matrix_free(v);
  gsl_vector_free(work);
  gsl_vector_free(s);
  return;
}



// inverse of symmetric matrix via Cholesky decomposition
int inverse_symm(gsl_matrix *in,gsl_matrix *out)
{
  #include<gsl/gsl_errno.h>
  gsl_set_error_handler_off();

  int i,flag=0;
  const int dim=(int)in->size1;
  gsl_matrix *check=gsl_matrix_alloc(dim,dim);
  gsl_vector *b=gsl_vector_alloc(dim);
  gsl_matrix_memcpy(check,in);  // do not overwrite 'in'

  flag=gsl_linalg_cholesky_decomp(check);

  if (!flag) {
    for (i=0;i<dim;i++) {
      gsl_vector_set_basis(b,i);
      gsl_linalg_cholesky_svx(check,b);
      gsl_matrix_set_col(out,i,b);
    }
  }
  gsl_matrix_free(check);
  gsl_vector_free(b);
  return(flag);
}


int fisher_ellipse(gsl_matrix *fisher,int index1,int index2,double *semi_major,double *semi_minor,double *angle)
{
  int flag;

  // setup
  if (index1==index2) return(2);  // can only make ellipse for two different indices
  const int dim=fisher->size1;
  gsl_matrix *inv=gsl_matrix_alloc(dim,dim);
  gsl_matrix *sub=gsl_matrix_alloc(2,2);
  gsl_matrix *v=gsl_matrix_alloc(2,2);
  gsl_vector *work=gsl_vector_alloc(2);
  gsl_vector *s=gsl_vector_alloc(2);

  // invert Fisher matrix
  flag=inverse_symm(fisher,inv);

  // get ellipse plotting parameters
  gsl_matrix_set(sub,0,0,gsl_matrix_get(inv,index1,index1));
  gsl_matrix_set(sub,0,1,gsl_matrix_get(inv,index1,index2));
  gsl_matrix_set(sub,1,0,gsl_matrix_get(inv,index2,index1));
  gsl_matrix_set(sub,1,1,gsl_matrix_get(inv,index2,index2));

  gsl_linalg_SV_decomp(sub,v,s,work);
  *semi_major=sqrt(gsl_vector_get(s,0));
  *semi_minor=sqrt(gsl_vector_get(s,1));
  *angle=atan(gsl_matrix_get(v,1,0)/gsl_matrix_get(v,0,0)); 

  // clean up
  gsl_vector_free(work);
  gsl_matrix_free(v);
  gsl_matrix_free(inv);
  gsl_matrix_free(sub);
  gsl_vector_free(s);
  return(flag);
}



// provides KDE in 1D with Gaussians, using Silverman's rule if width<=0 on call
void kde_1D(double *dp,int nd,int nstep,double step_min,double step_max,double h,double *steps,double *func)
{
  int i,j;
  double width=1.;
  const double gauss_norm=1./sqrt(2.*PI);

  // set binning
  for (i=0;i<nstep;i++) {
    steps[i]=step_min+i*(step_max-step_min)/(1.*nstep-1.);
  }

  // get width of Gaussian
  if (h<=0.0) {  // apply Silverman's rule
    width=1.06*gsl_stats_sd(dp,1,nd)*pow(nd,-0.2);
  }
  else width=h;

  // perform KDE
  for (i=0;i<nstep;i++) {
    func[i]=0.0;
    for (j=0;j<nd;j++) {
      func[i]+=exp(-0.5*(steps[i]-dp[j])*(steps[i]-dp[j])/(width*width))/width;
    }
    func[i]*=gauss_norm;  // normalised to nd if integrated
  }

  return;
}
