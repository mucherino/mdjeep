/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - 
              management of matrices B and Q
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    April 16 2010  v.0.1  first release
*****************************************************************************/

#include "bp.h"

// setting the cumulative matrix Q

void set_matrix(double b11,double b12,double b13,double b14,
                double b21,double b22,double b23,double b24,
                double b31,double b32,double b33,double b34,
                double b41,double b42,double b43,double b44,  // the last row is actually never used
                double **Q,int i)
{
   int k;
   double aux[12];

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


