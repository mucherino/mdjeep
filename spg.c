/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - SPG
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version
****************************************************************************************************/

#include "bp.h"

// SPG parameters: the user cannot modify these parameters in version 0.3.0
extern int K;  // fixed to 3 in this version of MDjeep
double eta = 0.99;
double gam = 1.e-4;
double epsobj = 1.e-7;
double epsg = 1.e-8;
double epsalpha = 1.e-12;
double mumin = 1.e-12;
double mumax = 1.e+12;

// this function computes the scalar product between two pairs (X1,y1) and (X2,y2)
// where X* are matrices, and y* are vectors
double scalarProd(int n,double **X1,double **X2,int m,double *y1,double *y2)
{
   int i,j,k;
   double prod = 0.0;

   for (k = 0; k < K; k++)  for (i = 0; i < n; i++)  prod = prod + (X1[k][i] * X2[k][i]);
   for (j = 0; j < m; j++)  prod = prod + (y1[j]*y2[j]);

   return prod;
};

// this function computes the norm for pair (X,y)
double norm(int n,double **X,int m,double *y)
{
   return sqrt(scalarProd(n,X,X,m,y,y));
};

/* Spectral Projected Gradient (SPG)
 *
 *           input: the DGP instance (n,v), and the starting point X
 *          output: the found solution replaces the starting point in X
 *                  the stress function value in the found solution (obj, pointer)
 *                  the number of iterations (it, pointer)
 * returning value: the flag indicating the termination status (0 = normal, 
 *                                                              1 = direction norm too small,
 *                                                              2 = max number of iterations)
 * Additional memory and parameters in the SEARCH structure S; all memory needs to be pre-allocated.
 */
int spg(int n,VERTEX *v,double **X,SEARCH S,int *its,double *obj)
{
   int i,j,k;
   int m;
   int it,maxIt;
   REFERENCE *ref;
   double mu,alpha;
   double C,Q;
   double objval,newobjval;
   double scalprod;
   double dist;
   short flag = 0;

   // the max number of iterations depends on the problem size
   maxIt = 50 + 10*n;

   // fixing y variables at the centers of the given interval distances
   m = 0;
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         j = ref->otherId;
         dist = distance(j,i,X);
         S.y[m] = projection(dist,lowerBound(ref),upperBound(ref),gam);
         ref = ref->next;
         m++;
      };
   };

   // computing initial objective function and gradient values
   objval = compute_stress(n,v,X,S.y);
   stress_gradient(n,v,X,S.y,S.gX,S.gy,S.memory);
   C = objval;

   // running Spectral Projected Gradient Descent method
   it = 1;  Q = 1;  alpha = 1.0;
   while (maxIt > it && objval > epsobj && alpha > epsalpha)
   {
      // computing spectral parameter
      if (it == 1)
      {
         mu = 1.0;
      }
      else
      {
         differenceMatrix(K,n,S.gX,S.gXp,S.YX);  differenceVector(m,S.gy,S.gyp,S.Yy);
         differenceMatrix(K,n,X,S.Xp,S.ZX);  differenceVector(m,S.y,S.yp,S.Zy);
         mu = scalarProd(n,S.YX,S.ZX,m,S.Yy,S.Zy) / scalarProd(n,S.ZX,S.ZX,m,S.Zy,S.Zy);
         if (mu < mumin)  mu = mumin;
         if (mu > mumax)  mu = mumax;
      };

      // making a full step over the opposite direction of the gradient
      for (i = 0; i < n; i++)
      {
         for (k = 0; k < K; k++)
         {
            S.sX[k][i] = X[k][i] - S.gX[k][i]/mu;
         };
      };
      for (j = 0; j < m; j++)  S.sy[j] = S.y[j] - S.gy[j]/mu;

      // performing projection on the box constraints (x variables)
      for (i = 0; i < n; i++)
      {
         for (k = 0; k < K; k++)
         {
            S.sX[k][i] = projection(S.sX[k][i],S.lX[k][i]-S.be,S.uX[k][i]+S.be,gam);
         };
      };

      // performing projection on the box constraints (y variables)
      k = 0;
      for (i = 0; i < n; i++)
      {
         ref = v[i].ref;
         while (ref != NULL)
         {
            S.sy[k] = projection(S.sy[k],lowerBound(ref),upperBound(ref),gam);
            ref = ref->next;
            k++;
         };
      };

      // computing new descent direction D
      for (i = 0; i < n; i++)
      {
         for (k = 0; k < K; k++)   
         {
            S.DX[k][i] = S.sX[k][i] - X[k][i];
         };
      };
      for (j = 0; j < m; j++)  S.Dy[j] = S.sy[j] - S.y[j];

      if (norm(n,S.DX,m,S.Dy) < epsg)
      {
         flag = 1;
         break;
      };

      // performing nonmonotone line-search
      alpha = 2.0;
      copyMatrix(K,n,X,S.Xp); copyVector(m,S.y,S.yp);
      copyMatrix(K,n,S.gX,S.gXp);  copyVector(m,S.gy,S.gyp);
      scalprod = scalarProd(n,S.gX,S.DX,m,S.gy,S.Dy);
      do
      {
         alpha = 0.5*alpha;
         for (i = 0; i < n; i++)
         {
            for (k = 0; k < K; k++)  X[k][i] = S.Xp[k][i] + alpha*S.DX[k][i];
         };
         for (j = 0; j < m; j++)  S.y[j] = S.yp[j] + alpha*S.Dy[j];

         newobjval = compute_stress(n,v,X,S.y);
      }
      while (alpha > epsalpha && newobjval > C + gam*alpha*scalprod);

      if (alpha <= epsalpha)  scalprod = scalprod/(norm(n,S.gX,m,S.gy)*norm(n,S.DX,m,S.Dy));
      newobjval = compute_stress(n,v,X,S.y);

      // preparing for next iteration
      C = eta*Q*C;
      Q = eta*Q + 1.0;
      C = (C + newobjval)/Q;
      objval = newobjval;
      stress_gradient(n,v,X,S.y,S.gX,S.gy,S.memory);

      it++;
   };

   // updating output info
   *its = it;
   *obj = objval;
   if (it == maxIt)  flag = 2;

   return flag;
};

