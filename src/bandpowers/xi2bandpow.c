// xi2bandpow.c
// 31.08.2018
// converts xi_+/-, xi_g+/x, and xi_gg to bandpowers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>
#include "bjutils.h"

#define xi_pm_balance 0.5 // default weight for xi_+/- in P_E

void xi_to_derived2pt();
double K_bandpow_ee();
double K_bandpow_ne();
double K_bandpow_nn();
double apod_upper();
double apod_lower();



int main(int argc, char *argv[])
{
  int i,j;
  double dr,xip,xim,scale_terms,kernel;
  FILE *dat;
  char name[400],path[200],ident[200],outident[200];

  if (argc!=14) {
    printf("Syntax: %s\n1: <working directory>\n2: <input file identifier>\n3: <output file identifier>\n4: <number of input angular bins>\n5: <min input separation to use in conversion [arcmin] (xi_+ in case of ee)>\n6: <max input separation to use in conversion [arcmin] (xi_+ in case of ee)>\n7: <min input separation to use in conversion [arcmin] (xi_- in case of ee; otherwise unused)>\n8: <max input separation to use in conversion [arcmin] (xi_- in case of ee; otherwise unused)>\n9: <number of output angular bins>\n10: <min output angular frequency>\n11: <max output angular frequency>\n12: <correlation type (1: ee; 2: ne; 3: gg)>\n13: <log width of apodisation window [total width of apodised range is tmax/tmin=exp(width) in arcmin; <0 for no apodisation]>\n\n",argv[0]);
    exit(-1);
  }
  sprintf(path,"%s",argv[1]);
  sprintf(ident,"%s",argv[2]);
  sprintf(outident,"%s",argv[3]);
  const int NR=atoi(argv[4]); 
  const double angmin1=atof(argv[5]); 
  const double angmax1=atof(argv[6]); 
  const double angmin2=atof(argv[7]); 
  const double angmax2=atof(argv[8]); 
  const int NOUT=atoi(argv[9]); 
  const double ellmin=atof(argv[10]); 
  const double ellmax=atof(argv[11]); 
  const int corrtype=atoi(argv[12]); 
  const double logwidth=atof(argv[13]); 


  // allocations
  if (corrtype>3) {
    printf("Error: incorrect input for correlation type!\n");
    exit(-1);
  }

  double *tin=bj_alloc(NR);
  double *ell_centre=bj_alloc(NOUT);
  double *ell_bound=bj_alloc(NOUT+1);
  double *weight=bj_alloc(NR);
  double *weight_xim=bj_alloc(NR);
  double *pfrac=bj_alloc(NOUT);

  gsl_vector *xi=gsl_vector_calloc(2*NR);
  gsl_vector *derived2pt=gsl_vector_calloc(2*NOUT);
  gsl_matrix *trans=gsl_matrix_calloc(2*NOUT,2*NR);
  gsl_matrix *xierr=gsl_matrix_calloc(2*NR,2*NR);
  gsl_matrix *derived2pterr=gsl_matrix_calloc(2*NOUT,2*NOUT);

  double logbinwidth=log(ellmax/ellmin)/(1.*NOUT);  // output binning
  for (i=0;i<NOUT;i++) {
    ell_centre[i]=ellmin*exp((i+0.5)*logbinwidth);   // use log centre of bin
    ell_bound[i]=ellmin*exp(i*logbinwidth);
  }
  ell_bound[NOUT]=ellmin*exp(NOUT*logbinwidth);


  // read input data
  sprintf(name,"%s/xi2bandpow_input_%s.dat",path,ident);
  dat=bj_fileopen("r",name);
  for (i=0;i<NR;i++) {
    if (corrtype==3) {
      fscanf(dat,"%lf  %lf",&tin[i],&xip);  //theta (bin centres; in arcmin), w_gg
      xim=0.0;   // use only first half of vector
    }
    else {
      fscanf(dat,"%lf  %lf  %lf",&tin[i],&xip,&xim);  //theta (bin centres), xi_+, xi_- for ee; theta (bin centres), xi_g+, xi_gx for ne
    }

    gsl_vector_set(xi,i,xip);
    gsl_vector_set(xi,NR+i,xim);
    gsl_matrix_set(xierr,i,i,xip*xip);       // treat xi as noise, uncorrelated so diagonal
    gsl_matrix_set(xierr,NR+i,NR+i,xim*xim);

    if (logwidth<0) {   // no apodisation
      if ((tin[i]>angmin1)&&(tin[i]<angmax1)) {
	weight[i]=1.0;
      }
      else {
	weight[i]=0.0;
      }
    }
    else {
      weight[i]=apod_lower(tin[i],angmin1,logwidth)*apod_upper(tin[i],angmax1,logwidth);
    }

    if (corrtype==1) { //ee
      if (logwidth<0) {   // no apodisation
	if ((tin[i]>angmin2)&&(tin[i]<angmax2)) {
	  weight_xim[i]=1.0;
	}
	else {
	  weight_xim[i]=0.0;
	}
      }
      else {
	weight_xim[i]=apod_lower(tin[i],angmin2,logwidth)*apod_upper(tin[i],angmax2,logwidth);
      }
    }

    tin[i]*=PI/10800.;   // convert from arcmin to radian
  }  // end of i loop
  fclose(dat);

  if ((tin[NR-1]<angmax1*PI/10800.*exp(logwidth/2.))||(corrtype==1 && tin[NR-1]<angmax2*PI/10800.*exp(logwidth/2.))) {
    printf("Error: correlation function not available out to max. range requested: %g < %g ; %g\n",tin[NR-1]/PI*10800.,angmax1*exp(logwidth/2.),angmax2*PI/10800.*exp(logwidth/2.));
    exit(-1);
  }
  if ((tin[0]>angmin1*PI/10800.*exp(-logwidth/2.))||(corrtype==1 && tin[0]>angmin2*PI/10800.*exp(-logwidth/2.))) {
    printf("Error: correlation function not available out to min. range requested: %g > %g ; %g\n",tin[0]/PI*10800.,angmin1*exp(-logwidth/2.),angmin2*PI/10800.*exp(-logwidth/2.));
    exit(-1);
  }

  sprintf(name,"%s/xi2bandpow_pmweights_%s.dat",path,ident);
  if ((dat=fopen(name,"r"))==NULL) {
    printf("No file %s - use pmweight of %g instead.\n",name,xi_pm_balance);
    for (i=0;i<NOUT;i++) {  
      pfrac[i]=xi_pm_balance;
    }
  }
  else {
    printf("Using pmweights from file %s\n",name);
    for (i=0;i<NOUT;i++) { 
      fscanf(dat,"%lf",&pfrac[i]);  // xi_+/- weight in P_E; one value per output ell bin
      printf("K_+ for Bin %i: %g\n",i+1,pfrac[i]);
    }
    fclose(dat);
  }


  // construct transformation matrix
  logbinwidth=log(tin[NR-1]/tin[0])/(1.*NR-1.);
  for (j=0;j<NR;j++) {
    dr=tin[j]*(exp(0.5*logbinwidth)-exp(-0.5*logbinwidth));  // input bin width
    for (i=0;i<NOUT;i++) {         
      scale_terms=2.*PI*dr/(tin[j]*log(ell_bound[i+1]/ell_bound[i]));

      kernel=0.0;
      if (corrtype==1) {
	kernel=K_bandpow_ee(tin[j]*ell_bound[i+1],0)-K_bandpow_ee(tin[j]*ell_bound[i],0);
	gsl_matrix_set(trans,i,j,pfrac[i]*weight[j]*kernel*scale_terms);
	gsl_matrix_set(trans,NOUT+i,j,0.5*weight[j]*kernel*scale_terms);

	kernel=K_bandpow_ee(tin[j]*ell_bound[i+1],1)-K_bandpow_ee(tin[j]*ell_bound[i],1);
	gsl_matrix_set(trans,i,NR+j,(1.-pfrac[i])*weight_xim[j]*kernel*scale_terms);
	gsl_matrix_set(trans,NOUT+i,NR+j,-0.5*weight_xim[j]*kernel*scale_terms);
      }
      if (corrtype==2) {
	kernel=K_bandpow_ne(tin[j]*ell_bound[i+1])-K_bandpow_ne(tin[j]*ell_bound[i]);
	gsl_matrix_set(trans,i,j,weight[j]*kernel*scale_terms);
	gsl_matrix_set(trans,NOUT+i,NR+j,weight[j]*kernel*scale_terms);
      }
      if (corrtype==3) {
	kernel=K_bandpow_nn(tin[j]*ell_bound[i+1])-K_bandpow_nn(tin[j]*ell_bound[i]);
	gsl_matrix_set(trans,i,j,weight[j]*kernel*scale_terms);
      }
    }
  }


  // convert signals
  xi_to_derived2pt(xi,xierr,trans,derived2pt,derived2pterr,2*NR,2*NOUT);


  // output band powers
  sprintf(name,"%s/xi2bandpow_output_%s.dat",path,outident);
  dat=bj_fileopen("w",name);
  fprintf(dat,"#   ell         PEE             PEEerr          PBB                PBBerr\n");
  for (i=0;i<NOUT;i++) {
    if (corrtype==3) {
      fprintf(dat,"%15.10g\t%15.10g\t%15.10g\n",ell_centre[i],gsl_vector_get(derived2pt,i),sqrt(gsl_matrix_get(derived2pterr,i,i)));   // ell, bandpow, err
    }
    else {
      fprintf(dat,"%15.10g\t%15.10g\t%15.10g\t%15.10g\t%15.10g\n",ell_centre[i],gsl_vector_get(derived2pt,i),sqrt(gsl_matrix_get(derived2pterr,i,i)),gsl_vector_get(derived2pt,NOUT+i),sqrt(gsl_matrix_get(derived2pterr,NOUT+i,NOUT+i)));   // ell, E/+, err, B/x, err
    }
  }
  fclose(dat);

  sprintf(name,"%s/xi2bandpow_kernels_%s.dat",path,outident);  // write kernels
  dat=bj_fileopen("w",name);
  for (j=0;j<NR;j++) {
    fprintf(dat,"%15.10g\t",tin[j]/PI*10800.);
    for (i=0;i<NOUT;i++) {
      fprintf(dat,"%15.10g\t",gsl_matrix_get(trans,i,j));   // E-mode, xi_+
    }
    for (i=0;i<NOUT;i++) {
      fprintf(dat,"%15.10g\t",gsl_matrix_get(trans,i,NR+j));   // E-mode, xi_-
    }
    for (i=0;i<NOUT;i++) {
      fprintf(dat,"%15.10g\t",gsl_matrix_get(trans,NOUT+i,j));   // B-mode, xi_+
    }
    for (i=0;i<NOUT;i++) {
      fprintf(dat,"%15.10g\t",gsl_matrix_get(trans,NOUT+i,NR+j));   // B-mode, xi_-
    }
    fprintf(dat,"\n");
  }
  fclose(dat);


  // clean up
  gsl_vector_free(xi);
  gsl_vector_free(derived2pt);
  gsl_matrix_free(trans);
  gsl_matrix_free(xierr);
  gsl_matrix_free(derived2pterr);
  free(tin);
  free(ell_centre);
  free(ell_bound);
  free(weight);
  free(weight_xim);
  free(pfrac);
  return 0;
}




void xi_to_derived2pt(gsl_vector *xi,gsl_matrix *xierr,gsl_matrix *trans,gsl_vector *res,gsl_matrix *reserr,int nin,int nout)
{
  gsl_matrix *tmp=gsl_matrix_calloc(nout,nin);
  gsl_matrix *transt=gsl_matrix_calloc(nin,nout);

  // signal conversion
  mattimesvec_gen(trans,xi,res,nout,nin); 

  // error conversion
  mattimesmat_gen(trans,xierr,tmp,nout,nin,nin);
  matrix_transpose(trans,transt,nout,nin);
  mattimesmat_gen(tmp,transt,reserr,nout,nin,nout);

  gsl_matrix_free(tmp);
  gsl_matrix_free(transt);
  return;
}



double K_bandpow_ee(double x,int pm)
{
  double res=0.0;

  if (pm==0) {
    res=x*gsl_sf_bessel_Jn(1,x);
  }
  else if (pm==1) {
    res=(x-8./x)*gsl_sf_bessel_Jn(1,x)-8.*gsl_sf_bessel_Jn(2,x);
  }
  else {
    printf("Error in routine K_bandpow_ee: invalid option for 'pm'.\n");
    exit(-1);
  }

  return(res);
}


double K_bandpow_ne(double x)
{
  return((-1.)*x*gsl_sf_bessel_Jn(1,x)-2.*gsl_sf_bessel_Jn(0,x));
}


double K_bandpow_nn(double x)
{
  return(x*gsl_sf_bessel_Jn(1,x));
}


double apod_upper(double x,double scale,double logwidth)
{
  double res=0.0;
  double logx=log(x);
  double logxmin=log(scale)-logwidth/2.;
  double logxmax=log(scale)+logwidth/2.;

  if (logx<=logxmin) {
    res=1.;
  }
  else if (logx>logxmax) {
    res=0.;
  }
  else {
    res=cos(PI/2.*(logx-logxmin)/(logxmax-logxmin));
    res=res*res;  // cos^2
  }
  return(res);
}


double apod_lower(double x,double scale,double logwidth)
{
  double res=0.0;
  double logx=log(x);
  double logxmin=log(scale)-logwidth/2.;
  double logxmax=log(scale)+logwidth/2.;

  if (logx<=logxmin) {
    res=0.;
  }
  else if (logx>logxmax) {
    res=1.;
  }
  else {
    res=cos(PI/2.*(logx-logxmax)/(logxmax-logxmin));
    res=res*res;  // cos^2
  }
  return(res);
}
