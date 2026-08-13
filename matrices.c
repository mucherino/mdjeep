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
	      Aug 13 2026  v.0.3.3  new function to evaluate the RMSD between two solutions
**************************************************************************************************/ 

#include "bp.h"

// this function allocates memory for a vector (1-dim array of double)
double* allocateVector(size_t n)
{
   return (double*)calloc(n,sizeof(double));
};

// this function copies a vector into another
void copyVector(size_t n,double *source,double *dest)
{
   for (size_t i = 0; i < n; i++)  dest[i] = source[i];
};

// this function computes the difference between two vectors of the same length
void differenceVector(size_t n,double *a,double *b,double *c)
{
   for (size_t i = 0; i < n; i++)  c[i] = a[i] - b[i];
};

// this function computes the norm of a given vector
double normVector(size_t n,double *v)
{
   double sqnorm = 0.0;
   for (size_t i = 0; i < n; i++)  sqnorm = sqnorm + v[i]*v[i];
   return sqrt(sqnorm);
};

// this function verifies whether two vectors contain the same sequence of values
bool areSameVector(size_t n,double *v1,double *v2)
{
   bool same = true;
   size_t i = 0;
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
   for (size_t i = 0; i < n; i++)  printf(" %20.17lf",v[i]);
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
   double **A = (double**)calloc(n,sizeof(double));
   for (size_t i = 0; i < n; i++)  A[i] = (double*)calloc(m,sizeof(double));
   return A;
};

// this function sets all matrix elements to zero (NEW)
void zeroMatrix(size_t n,size_t m,double **A)
{
   for (size_t i = 0; i < n; i++)  for (size_t j = 0; j < m; j++)  A[i][j] = 0.0;
};

// this function creates a new matrix with the squared elements of another one
void squaredElementsMatrix(size_t n,size_t m,double **A,double **SA)
{
   for (size_t i = 0; i < n; i++)  for (size_t j = 0; j < m; j++)  SA[i][j] = A[i][j]*A[i][j];
};

// this function copies a matrix into another
void copyMatrix(size_t n,size_t m,double **source,double **dest)
{
   for (size_t i = 0; i < n; i++)  for (size_t j = 0; j < m; j++)  dest[i][j] = source[i][j];
};

// this function copies and centers a matrix
void copyCenterMatrix(size_t n,size_t m,double **source,double **dest)
{
   double sum;
   if (m == 0)  return;
   for (size_t i = 0; i < n; i++)
   {
      sum = 0.0;
      for (size_t j = 0; j < m; j++)  sum = sum + source[i][j];
      for (size_t j = 0; j < m; j++)  dest[i][j] = source[i][j] - sum/m;
   };
};

// this function computes the difference between two matrices
void differenceMatrix(size_t n,size_t m,double **A,double **B,double **C)
{
   for (size_t i = 0; i < n; i++)  for (size_t j = 0; j < m; j++)  C[i][j] = A[i][j] - B[i][j];
};

// this function verifies whether two matrices are identical
bool areSameMatrix(size_t n,size_t m,double **A,double **B)
{
   bool same = true;
   size_t i,j;

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

// this function computes the sum of all matrix elements
double sumElementMatrix(size_t n,size_t m,double **A)
{
   double sum = 0.0;
   for (size_t i = 0; i < n; i++)  for (size_t j = 0; j < m; j++)  sum = sum + A[i][j];
   return sum;
};

// this function computes the determinant of a 3x3 matrix
double determinant3x3Matrix(double **A)
{
   return A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] -
          A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1];
};

// this function computes the determinant of a 4x4 matrix
double determinant4x4Matrix(double **A)
{
   return A[0][0] * (A[1][1]*A[2][2]*A[3][3] + A[1][2]*A[2][3]*A[3][1] + A[1][3]*A[2][1]*A[3][2] - A[1][3]*A[2][2]*A[3][1] - A[1][2]*A[2][1]*A[3][3] - A[1][1]*A[2][3]*A[3][2]) -
          A[0][1] * (A[1][0]*A[2][2]*A[3][3] + A[1][2]*A[2][3]*A[3][0] + A[1][3]*A[2][0]*A[3][2] - A[1][3]*A[2][2]*A[3][0] - A[1][2]*A[2][0]*A[3][3] - A[1][0]*A[2][3]*A[3][2]) +
          A[0][2] * (A[1][0]*A[2][1]*A[3][3] + A[1][1]*A[2][3]*A[3][0] + A[1][3]*A[2][0]*A[3][1] - A[1][3]*A[2][1]*A[3][0] - A[1][1]*A[2][0]*A[3][3] - A[1][0]*A[2][3]*A[3][1]) -
          A[0][3] * (A[1][0]*A[2][1]*A[3][2] + A[1][1]*A[2][2]*A[3][0] + A[1][2]*A[2][0]*A[3][1] - A[1][2]*A[2][1]*A[3][0] - A[1][1]*A[2][0]*A[3][2] - A[1][0]*A[2][2]*A[3][1]);
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

// this function computes the RMSD between two conformations
// D.L. Theobald, Acta Crystallographica Section A 61(4), 478-480, 2005.
double rmsd(size_t n,double **coordsA,double **coordsB,size_t maxit)
{
   unsigned short x = 0;
   unsigned short y = 1;
   unsigned short z = 2;
   unsigned short three = 3;
   unsigned short four = 4;

   // safeguard against dividing by 0
   if (n == 0)  return 0.0;

   // centering the coordinates
   double **centerA = allocateMatrix(three,n);
   copyCenterMatrix(three,n,coordsA,centerA);
   double **centerB = allocateMatrix(three,n);
   copyCenterMatrix(three,n,coordsB,centerB);

   // initializing the inner product matrix M
   double **M = allocateMatrix(three,three);
   zeroMatrix(three,three,M);

   // computing the inner products
   double GA = 0.0;
   double GB = 0.0;
   for (size_t k = 0; k < n; k++)
   {
      for (unsigned short i = 0; i < three; i++)
      {
         GA = GA + centerA[i][k]*centerA[i][k];
         GB = GB + centerB[i][k]*centerB[i][k];
         for (unsigned short j = 0; j < three; j++)
         {
            M[i][j] = M[i][j] + centerA[i][k]*centerB[j][k];
         };
      };
   };
   double lambda = 0.5*(GA + GB);  // initial guess for max eigenvalue

   // computing the matrix of squared elements
   double **SM = allocateMatrix(three,three);
   squaredElementsMatrix(three,three,M,SM);

   // computing the key matrix K
   double **K = allocateMatrix(four,four);
   K[0][0] = M[x][x] + M[y][y] + M[z][z];
   K[0][1] = M[y][z] - M[z][y];
   K[0][2] = M[z][x] - M[x][z];
   K[0][3] = M[x][y] - M[y][x];
   K[1][0] = M[y][z] - M[z][y];
   K[1][1] = M[x][x] - M[y][y] - M[z][z];
   K[1][2] = M[x][y] + M[y][x];
   K[1][3] = M[z][x] + M[x][z];
   K[2][0] = M[z][x] - M[x][z];
   K[2][1] = M[x][y] + M[y][x];
   K[2][2] = M[y][y] - M[z][z] - M[x][x];
   K[2][3] = M[y][z] + M[z][y];
   K[3][0] = M[x][y] - M[y][x];
   K[3][1] = M[z][x] + M[x][z];
   K[3][2] = M[y][z] + M[z][y];
   K[3][3] = M[z][z] - M[x][x] - M[y][y];

   // computing the coefficients C2 and C1 of the characteristic polynomial
   double C2 = -2.0*sumElementMatrix(three,three,SM);
   double C1 = -8.0*determinant3x3Matrix(M);
   double C0 = determinant4x4Matrix(K);

   // freeing temporary matrices
   freeMatrix(four,K);
   freeMatrix(three,SM);
   freeMatrix(three,M);
   freeMatrix(three,centerB);
   freeMatrix(three,centerA);

   // Newton-Raphson method
   size_t it = 0;
   double lambda0 = 0.0;
   do {
      lambda0 = lambda;
      double lambda2 = lambda0*lambda0;
      double poly = lambda2*lambda2 + C2*lambda2 + C1*lambda0 + C0;
      double dpoly = 4.0*lambda2*lambda0 + 2.0*C2*lambda0 + C1;
      if (fabs(dpoly) < 1e-15)  break;
      lambda = lambda0 - poly/dpoly;
      it++;
      if (it == maxit)  break;
   }
   while (fabs(lambda - lambda0) > fabs(1.e-9*lambda));

   // RMSD value
   double msd = (GA + GB - 2.0*lambda) / (double)n;
   if (msd < 0.0)  return 0.0;
   return sqrt(msd);
};

// this function prints a matrix
void printMatrix(size_t n,size_t m,double **a)
{
   for (size_t i = 0; i < n; i++)
   {
      for (size_t j = 0; j < m; j++)
      {
         printf(" %20.17lf",a[i][j]);
      };
      printf("\n");
   };
};

// this function frees a matrix
double** freeMatrix(size_t n,double **a)
{
   for (size_t i = 0; i < n; i++)  free(a[i]);
   free(a);
   return NULL;
};

