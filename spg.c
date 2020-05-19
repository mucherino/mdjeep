/***********************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - SPG
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version
              Mar 21 2020  v.0.3.1  the variable S.be is not used directly in SPG to enlarge the bounds
              May 19 2020  v.0.3.2  parameters are now in the OPTION structure
                                    features to monitor and print added (SPG may be invoked as a main method)
************************************************************************************************************/

#include "bp.h"

// the dimension is fixed to 3 in this version of MDjeep
int K = 3;

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
int spg(int n,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info,int *its,double *obj)
{
   int i,j,k;
   int m;
   int it,maxIt;
   int ldigits;
   REFERENCE *ref;
   double mu,alpha;
   double C,Q;
   double objval,newobjval;
   double scalprod;
   double dist;
   short flag = 0;

   // if spg is refinement method, max number of iterations depends on the problem size
   if (info->refinement == 1)
      maxIt = 50 + 10*n;
   else
      maxIt = op.maxit;

   // computing y variables
   m = 0;
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         j = otherVertexId(ref);
         dist = distance(j,i,X);
         S.y[m] = projection(dist,lowerBound(ref),upperBound(ref),op.gam);
         ref = ref->next;
         m++;
      };
   };

   // computing initial objective function and gradient values
   objval = compute_stress(n,v,X,S.y);
   stress_gradient(n,v,X,S.y,S.gX,S.gy,S.memory);
   C = objval;

   // running Spectral Projected Gradient Descent method
   it = 1;  Q = 1.0;  alpha = 1.0;
   while (maxIt > it && objval > op.epsobj && alpha > op.epsalpha)
   {
      // monitor
      if (info->method == 1 && op.monitor)
      {
         ldigits = numberOfDigits(it);
         for (k = 0; k < info->ndigits + 9; k++)  fprintf(stderr,"\b");
         for (k = 0; k < info->ndigits - ldigits; k++)  fprintf(stderr," ");
         fprintf(stderr,"%d %8.2le",it,objval);
      };

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
         if (mu < op.mumin)  mu = op.mumin;
         if (mu > op.mumax)  mu = op.mumax;
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
            S.sX[k][i] = projection(S.sX[k][i],S.lX[k][i],S.uX[k][i],op.gam);
         };
      };

      // performing projection on the box constraints (y variables)
      k = 0;
      for (i = 0; i < n; i++)
      {
         ref = v[i].ref;
         while (ref != NULL)
         {
            S.sy[k] = projection(S.sy[k],lowerBound(ref),upperBound(ref),op.gam);
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

      if (norm(n,S.DX,m,S.Dy) < op.epsg)
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
      while (alpha > op.epsalpha && newobjval > C + op.gam*alpha*scalprod);

      if (alpha <= op.epsalpha)  scalprod = scalprod/(norm(n,S.gX,m,S.gy)*norm(n,S.DX,m,S.Dy));
      newobjval = compute_stress(n,v,X,S.y);

      // preparing for next iteration
      C = op.eta*Q*C;
      Q = op.eta*Q + 1.0;
      C = (C + newobjval)/Q;
      objval = newobjval;
      stress_gradient(n,v,X,S.y,S.gX,S.gy,S.memory);

      it++;
   };

   // printing (optional)
   if (info->method == 1)
   {
      if (op.print > 0)
      {
         if (op.format == 0)
            printfile(i,v,X,info->output,0);
         else
            printpdb(i,v,X,info->output,0);
      };
   };

   // updating output info
   *its = it;
   *obj = objval;
   if (it == maxIt)  flag = 2;

   return flag;
};

