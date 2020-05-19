/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - objective functions
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version
              Mar 21 2020  v.0.3.1  no changes
              May 19 2020  v.0.3.2  no changes
****************************************************************************************************/

#include "bp.h"

extern int K;  // the function "stress_gradient" uses the parameter K (space dimension);
               // in this version of MDjeep, K is always fixed to 3.

// Mean Distance Error (MDE)
// given a VERTEX array (n,v) and a realization X, this function computes the MDE value
// (eps is the tolerance to discriminate between exact and interval distances)
double compute_mde(int n,VERTEX *v,double **X,double eps)
{
   int i,j,m;
   REFERENCE *ref;
   double avg,dist;
   double value = 0.0;

   m = 0;
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         j = otherVertexId(ref);
         dist = distance(j,i,X);
         if (isExactDistance(ref,eps))
         {
            value = value + fabs(dist - lowerBound(ref))/lowerBound(ref);
         }
         else
         {
            avg = 0.5*(lowerBound(ref) + upperBound(ref));
            if (dist < lowerBound(ref))
            {
               value = value + fabs(dist - lowerBound(ref))/avg;
            }
            else if (dist > upperBound(ref))
            {
               value = value + fabs(dist - upperBound(ref))/avg;
            };
         };
         ref = ref->next;
         m++;
      };
   };
   if (m > 0)  value = value/n;

   return value;
};

// Largest Distance Error (LDE)
// given a VERTEX array (n,v) and a realization X, this function computes the LDE value
// (eps is the tolerance to discriminate between exact and interval distances)
double compute_lde(int n,VERTEX *v,double **X,double eps)
{
   int i,j;
   REFERENCE *ref;
   double diff,dist;
   double max = 0.0;

   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         j = otherVertexId(ref);
         dist = distance(j,i,X);
         if (isExactDistance(ref,eps))
         {
            diff = fabs(dist - lowerBound(ref));
            if (diff > max)  max = diff;
         }
         else
         {
            if (dist < lowerBound(ref))
            {
               diff = fabs(dist - lowerBound(ref));
               if (diff > max)  max = diff;
            }
            else if (dist > upperBound(ref))
            {
               diff = fabs(dist - upperBound(ref));
               if (diff > max)  max = diff;
            }
         };
         ref = ref->next;
      };
   };

   return max;
};

// STRESS function
// given a VERTEX array (n,v), a realization X, and vector y of selected distances from the intervals [lb,ub],
// this function computes the stress function [Glunt at al, "Molecular Conformations from Distance Matrices", 1993]
// (eps is the tolerance to discriminate between exact and interval distances)
double compute_stress(int n,VERTEX *v,double **X,double *y)
{
   int i,j,h;
   REFERENCE *ref;
   double term;
   double sigma = 0.0;

   h = 0;
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         j = otherVertexId(ref);
         term = distance(j,i,X) - y[h];
         term = term*term;
         sigma = sigma + term;
         ref = ref->next;
         h++;
      };
   };

   return sigma;
};

// this function computes the gradient of the stress function (see above)
// output arguments: the gradient wrt the variables X (gX), and the gradient wrt the variables y (gy)
// (the "memory" space needs to have at least size n)
void stress_gradient(int n,VERTEX *v,double **X,double *y,double **gX,double *gy,double *memory)
{
   int i,j,k,h;
   double tmp;
   REFERENCE *ref;

   // cleaning memory space
   for (i = 0; i < n; i++)  memory[i] = 0.0;

   // initialization for gX
   for (k = 0; k < K; k++)
   {
      for (i = 0; i < n; i++)
      {
         gX[k][i] = 0.0;
      }
   };

   // computation of gy and gX (compact form, all steps in one, except case u==v)
   h = 0;
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         j = otherVertexId(ref);
         tmp = distance(j,i,X);
         gy[h] = -2.0*(tmp - y[h]);
         if (tmp > 0.0)
         {
            tmp = -y[h]/tmp;
            memory[i] = memory[i] + tmp + 1.0;
            memory[j] = memory[j] + tmp + 1.0;
            tmp = -2.0*(1.0 + tmp);
            for (k = 0; k < K; k++)
            {
               gX[k][i] = gX[k][i] + tmp*X[k][j];
               gX[k][j] = gX[k][j] + tmp*X[k][i];
            };
         };
         ref = ref->next;
         h++;
      };
   };

   // completing the computation of gX (case u==v)
   for (k = 0; k < K; k++)
   {
      for (i = 0; i < n; i++)
      {
         gX[k][i] = gX[k][i] + 2.0*memory[i]*X[k][i];
      };
   };
};

