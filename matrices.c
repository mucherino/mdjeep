/****************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry
              management of matrices for the generation of atomic coordinates
  Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
              May 10 2014  v.0.2  implementation of 2nd method for coordinate computation
*****************************************************************************************/

#include "bp.h"

// Generation of atomic coordinates by using the method based on cumulative matrices (option -m1)
// ----------------------------------------------------------------------------------------------

// setting up cumulative matrix Q

void set_matrix(double b11,double b12,double b13,double b14,
                double b21,double b22,double b23,double b24,
                double b31,double b32,double b33,double b34,
                double b41,double b42,double b43,double b44,  // the last row is actually never used
                double **Q,int i)
{
   int k;
   double aux[12];

   // the matrix is stored row by row
   // more convenient because last row is not explicitly computed

   aux[0]  = b11;  aux[1]  = b12;  aux[2]  = b13;  aux[3]  = b14;		
   aux[4]  = b21;  aux[5]  = b22;  aux[6]  = b23;  aux[7]  = b24;
   aux[8]  = b31;  aux[9]  = b32;  aux[10] = b33;  aux[11] = b34;
// aux[12] = b41;  aux[13] = b42;  aux[14] = b43;  aux[15] = b44;

   if (i > 1)  
   {
	   // cumulating the new matrix with the previous ones
	   matrix_prod(Q[i-1],aux,Q[i]);
   }
   else
   {
	   // this is the first matrix: Q(1)
	   for (k = 0; k < 12; k++)  Q[i][k] = aux[k];
   };
};

// performing the product between two matrices
// (the last row is always (0 0 0 1), and it is never used)

void matrix_prod(double *Q1,double *Q2,double *Qout)
{
   Qout[0]  = Q1[0]*Q2[0] + Q1[1]*Q2[4] + Q1[2]*Q2[8];           // Q2[12] = 0  ==>  Q1[3]*Q2[12] = 0
   Qout[1]  = Q1[0]*Q2[1] + Q1[1]*Q2[5] + Q1[2]*Q2[9];           // Q2[13] = 0  ==>  Q1[3]*Q2[13] = 0
   Qout[2]  = Q1[0]*Q2[2] + Q1[1]*Q2[6] + Q1[2]*Q2[10];          // Q2[14] = 0  ==>  Q1[3]*Q2[14] = 0
   Qout[3]  = Q1[0]*Q2[3] + Q1[1]*Q2[7] + Q1[2]*Q2[11] + Q1[3];  // Q2[15] = 1  ==>  Q1[3]*Q2[15] = Q1[3]

   Qout[4]  = Q1[4]*Q2[0] + Q1[5]*Q2[4] + Q1[6]*Q2[8];           // Q2[12] = 0  ==>  Q1[7]*Q2[12] = 0
   Qout[5]  = Q1[4]*Q2[1] + Q1[5]*Q2[5] + Q1[6]*Q2[9];           // Q2[13] = 0  ==>  Q1[7]*Q2[13] = 0
   Qout[6]  = Q1[4]*Q2[2] + Q1[5]*Q2[6] + Q1[6]*Q2[10];	         // Q2[14] = 0  ==>  Q1[7]*Q2[14] = 0
   Qout[7]  = Q1[4]*Q2[3] + Q1[5]*Q2[7] + Q1[6]*Q2[11] + Q1[7];  // Q2[15] = 1  ==>  Q1[7]*Q2[15] = Q1[7]

   Qout[8]  = Q1[8]*Q2[0] + Q1[9]*Q2[4] + Q1[10]*Q2[8];          // Q2[12] = 0  ==>  Q1[11]*Q2[12] = 0
   Qout[9]  = Q1[8]*Q2[1] + Q1[9]*Q2[5] + Q1[10]*Q2[9];          // Q2[13] = 0  ==>  Q1[11]*Q2[13] = 0
   Qout[10] = Q1[8]*Q2[2] + Q1[9]*Q2[6] + Q1[10]*Q2[10];         // Q2[14] = 0  ==>  Q1[11]*Q2[14] = 0
   Qout[11] = Q1[8]*Q2[3] + Q1[9]*Q2[7] + Q1[10]*Q2[11] + Q1[11];// Q2[15] = 1  ==>  Q1[11]*Q2[15] = Q1[11]

// Qout[12] = 0.0;
// Qout[13] = 0.0;
// Qout[14] = 0.0;
// Qout[15] = 1.0;
};


// Generation of atomic coordinates by using the method based on the change of basis (option -m2)
// ----------------------------------------------------------------------------------------------

// setting up matrix U (stored column by column)

void gen_U(int i,SOLUTION *sol,double *U)
{
   int i1,i2,i3,ii;
   double nxaxis,nyaxis,nzaxis;
   double v1[3],v2[3];

   // reference atoms
   ii = i - 1;
   i3 = sol[ii].ref3;
   i2 = sol[ii].ref2;
   i1 = sol[ii].ref1;

   // x axis (first column)
   v1[0] = sol[i1].x - sol[i2].x;  v1[1] = sol[i1].y - sol[i2].y;  v1[2] = sol[i1].z - sol[i2].z;
   v2[0] = sol[i3].x - sol[i2].x;  v2[1] = sol[i3].y - sol[i2].y;  v2[2] = sol[i3].z - sol[i2].z;
   nxaxis = norm(v1);
   U[0] = v1[0]/nxaxis;  U[1] = v1[1]/nxaxis;  U[2] = v1[2]/nxaxis;

   // z axis (third column)
   cross_prod(v1,v2,&U[6]);  nzaxis = norm(&U[6]);
   U[6] = U[6]/nzaxis;  U[7] = U[7]/nzaxis;  U[8] = U[8]/nzaxis;

   // y axis (second column)
   cross_prod(&U[6],&U[0],&U[3]);  nyaxis = norm(&U[3]);
   U[3] = U[3]/nyaxis;  U[4] = U[4]/nyaxis;  U[5] = U[5]/nyaxis;
};

// computing the coordinates (using U)

void gen_coords(int i,SOLUTION *sol,double *U,double di1i,double ctheta,double stheta,double comega,double somega)
{
   int i1,ii;
   double a[3];

   // reference atoms
   ii = i - 1;  i1 = sol[ii].ref1;

   // computing vector a (depends on angles)
   a[0] = -di1i*ctheta;
   a[1] =  di1i*stheta*comega;
   a[2] =  di1i*stheta*somega;

   // generation of the coordinates
   sol[ii].x = sol[i1].x + a[0]*U[0] + a[1]*U[3] + a[2]*U[6];
   sol[ii].y = sol[i1].y + a[0]*U[1] + a[1]*U[4] + a[2]*U[7];
   sol[ii].z = sol[i1].z + a[0]*U[2] + a[1]*U[5] + a[2]*U[8];
};

