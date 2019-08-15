// basic functions

#include "../include/bjutils.h"


FILE *bj_fileopen(char *rw,char *name)
{
  FILE *dat;
  if ((dat=fopen(name,rw))==NULL) {
    printf("Could not open file %s\n",name);
    exit(-1);
  }
  return dat;
}


// returns no. of columns in file
int bj_get_file_columns(char *name)
{
  FILE *dat;
  int ncol,read,cursor=0;
  char dummy[MAX_CHAR_FILE];
  char *line = calloc(MAX_CHAR_FILE,sizeof(char));

  dat=bj_fileopen("r",name);
  ncol=0;
  fgets(line,MAX_CHAR_FILE,dat);
  while (sscanf(line+cursor, "%s%n", dummy, &read)==1) {
    cursor+=read;
    ncol++;
  }
  free(line);
  fclose(dat);
  return(ncol);
}


// returns no. of rows in file
int bj_get_file_rows(char *name)
{
  FILE *dat;
  int nrow=0;
  dat=bj_fileopen("r",name);
  while ((fscanf(dat,"%*[^\n]"), fscanf(dat,"%*c"))!=EOF) nrow++; 
  fclose(dat);
  return(nrow);
}


double *bj_alloc(int dim)
{
  double *arr=calloc(dim,sizeof(double));
  return arr;
}


double **bj_alloc_2D(int dim1,int dim2)
{
  int i;
  double **arr=calloc(dim1,sizeof(double *));
  for (i=0;i<dim1;i++) {
    arr[i]=calloc(dim2,sizeof(double));
  }
  return arr;
}


double ***bj_alloc_3D(int dim1,int dim2,int dim3)
{
  int i,j;
  double ***arr=calloc(dim1,sizeof(double *));
  for (i=0;i<dim1;i++) {
    arr[i]=calloc(dim2,sizeof(double *));
    for (j=0;j<dim2;j++) {
      arr[i][j]=calloc(dim3,sizeof(double));
    }
  }
  return arr;
}


double ****bj_alloc_4D(int dim1,int dim2,int dim3,int dim4)
{
  int i,j,k;
  double ****arr=calloc(dim1,sizeof(double *));
  for (i=0;i<dim1;i++) {
    arr[i]=calloc(dim2,sizeof(double *));
    for (j=0;j<dim2;j++) {
      arr[i][j]=calloc(dim3,sizeof(double *));
      for (k=0;k<dim3;k++) {
	arr[i][j][k]=calloc(dim4,sizeof(double));
      }
    }
  }
  return arr;
}


void bj_free_2D(double **arr,int dim1,int dim2)
{
  int i;
  for (i=0;i<dim1;i++) {
    free(arr[i]);
  }
  free(arr);
  return;
}


void bj_free_3D(double ***arr,int dim1,int dim2,int dim3)
{
  int i,j;
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(arr[i][j]);
    }
    free(arr[i]);
  }
  free(arr);
  return;
}


void bj_free_4D(double ****arr,int dim1,int dim2,int dim3,int dim4)
{
  int i,j,k;
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
	free(arr[i][j][k]);
      }
      free(arr[i][j]);
    }
    free(arr[i]);
  }
  free(arr);
  return;
}

