/*************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - utilities
  Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1    first release
              May 08 2014  v.0.2    functions costheta and cosomega updated
              Jun 28 2019  v.0.3.0  adding functions for vector and matrix manipulation
              Mar 21 2020  v.0.3.1  adding functions areSameVector and areSameMatrix (for tests)
              May 19 2020  v.0.3.2  no changes
**************************************************************************************************/ 

#include "bp.h"

// this function allocates memory for a vector (1-dim array of double)
double* allocateVector(size_t n)
{
   return (double*) calloc(n,sizeof(double));
};

// this function copies a vector into another
void copyVector(size_t n,double *source,double *dest)
{
   int i;
   for (i = 0; i < n; i++)  dest[i] = source[i];
};

// this function computes the difference between two vectors of the same length
void differenceVector(size_t n,double *a,double *b,double *c)
{
   int i;
   for (i = 0; i < n; i++)  c[i] = a[i] - b[i];
};

// this function computes the norm of a given vector
double normVector(size_t n,double *v)
{
   int i;
   double sqnorm = 0.0;
   for (i = 0; i < n; i++)  sqnorm = sqnorm + v[i]*v[i];
   return sqrt(sqnorm);
};

// this function verifies whether two vectors contain the same sequence of values
bool areSameVector(size_t n,double *v1,double *v2)
{
   int i = 0;
   bool same = true;
   while (same && (i < n))
   {
      if (v1[i] != v2[i])  same = false;
      i++;
   }
   return same;
};

// this function computes the cross product between two 3D vectors
void crossProdVector(double *v1,double *v2,double *res)
{
   res[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]); 
   res[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
   res[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
};

// this function prints a vector
void printVector(size_t n,double *v)
{
   int i;
   for (i = 0; i < n; i++)  printf(" %20.17lf",v[i]);
   printf("\n");
};

// this function frees a vector
double* freeVector(double *v)
{
   free(v);
   return NULL;
};

// this function allocates memory for a matrix (2D array of double)
double** allocateMatrix(size_t n,size_t m)
{
   int i;
   double **A = (double**)calloc(n,sizeof(double));
   for (i = 0; i < n; i++)  A[i] = (double*)calloc(m,sizeof(double));
   return A;
};

// this function copies a matrix into another
void copyMatrix(size_t n,size_t m,double **source,double **dest)
{
   int i,j;
   for (i = 0; i < n; i++)  for (j = 0; j < m; j++)  dest[i][j] = source[i][j];
};

// this function copies and centers a matrix
void copyCenterMatrix(size_t n,size_t m,double **source,double **dest)
{
   int i,j;
   double sum;

   for (i = 0; i < n; i++)
   {
      sum = 0.0;
      for (j = 0; j < m; j++)  sum = sum + source[i][j];
      for (j = 0; j < m; j++)  dest[i][j] = source[i][j] - sum/m;
   };
};

// this function computes the difference between two matrices
void differenceMatrix(size_t n,size_t m,double **A,double **B,double **C)
{
   int i,j;
   for (i = 0; i < n; i++)  for (j = 0; j < m; j++)  C[i][j] = A[i][j] - B[i][j];
};

// this function verifies whether two matrices are identical
bool areSameMatrix(size_t n,size_t m,double **A,double **B)
{
   int i,j;
   bool same = true;

   i = 0;  j = 0;
   while (same && (i < n))
   {
      while (same & (j < m))
      {
         if (A[i][j] != B[i][j])  same = false;
         j++;
      };
      i++;
   };

   return same;
};

// this function computes the U matrix (stored column by column)
void UMatrix(int i3,int i2,int i1,int i,double **X,double *U)
{
   double nxaxis,nyaxis,nzaxis;
   double v1[3],v2[3];

   // x axis (first column)
   v1[0] = X[0][i1] - X[0][i2];  v1[1] = X[1][i1] - X[1][i2];  v1[2] = X[2][i1] - X[2][i2];
   v2[0] = X[0][i3] - X[0][i2];  v2[1] = X[1][i3] - X[1][i2];  v2[2] = X[2][i3] - X[2][i2];
   nxaxis = normVector(3,v1);
   U[0] = v1[0]/nxaxis;  U[1] = v1[1]/nxaxis;  U[2] = v1[2]/nxaxis;

   // z axis (third column)
   crossProdVector(v1,v2,&U[6]);  nzaxis = normVector(3,&U[6]);
   U[6] = U[6]/nzaxis;  U[7] = U[7]/nzaxis;  U[8] = U[8]/nzaxis;

   // y axis (second column)
   crossProdVector(&U[6],&U[0],&U[3]);  nyaxis = normVector(3,&U[3]);
   U[3] = U[3]/nyaxis;  U[4] = U[4]/nyaxis;  U[5] = U[5]/nyaxis;
};

// this function computs the coordinates of the current vertex i (using U matrix)
void genCoordinates(int i1,int i,double **X,double *U,double di1i,double ctheta,double stheta,double comega,double somega)
{
   double a[3];

   // computing vector a (depends on angles)
   a[0] = -di1i*ctheta;
   a[1] =  di1i*stheta*comega;
   a[2] =  di1i*stheta*somega;

   // generation of the coordinates
   X[0][i] = X[0][i1] + a[0]*U[0] + a[1]*U[3] + a[2]*U[6];
   X[1][i] = X[1][i1] + a[0]*U[1] + a[1]*U[4] + a[2]*U[7];
   X[2][i] = X[2][i1] + a[0]*U[2] + a[1]*U[5] + a[2]*U[8];
};

// this function prints a matrix
void printMatrix(size_t n,size_t m,double **a)
{
   int i,j;
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < m; j++)
      {
         printf(" %20.17lf",a[i][j]);
      };
      printf("\n");
   };
};

// this function frees a matrix
double** freeMatrix(size_t n,double **a)
{
   int i;
   for (i = 0; i < n; i++)  free(a[i]);
   free(a);
   return NULL;
};

