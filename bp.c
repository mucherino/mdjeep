/*******************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - function BP
  Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
              May 10 2014  v.0.2  use of 2nd method for coordinate computation
********************************************************************************************/

#include "bp.h"

double epsomega = 1.e-4;

int bp(int i,int n,int m,PROBL *p,int **ind,SOLUTION *sol,OPTION op,PARTIALDE *plde,double **Q,INFORMATION *info)
{
   int i3,i2,i1,k;
   int nlde,pru;
   double U[9];
   double dist34,lde;
   double ctheta,stheta;
   double cosomega12,sinomega1,sinomega2;

   // is this the first atom?
   if (i == 1)
   {
      // the first three atoms can be positioned by using the available information
      sol[0].x = 0.0;  sol[0].y = 0.0;  sol[0].z = 0.0;
      k = ind[1][2];
      sol[1].x = -p[k].l;  sol[1].y = 0.0;  sol[1].z = 0.0;
      sol[2].x = -p[k].l;
      k = ind[2][3];
      ctheta = sol[2].costheta;  stheta = sqrt(1.0 - ctheta*ctheta);
      sol[2].x = sol[2].x + p[k].l*ctheta;  sol[2].y = p[k].l*stheta;  sol[2].z = 0.0;

      // the Q matrices are computed because necessary for obtaining next atoms (when option -m1)
      if (op.method == 1)
      {
         // 1st matrix - corresponding to sol[1]
         k = ind[1][2];
         set_matrix( -1, 0, 0, -p[k].l, 
                      0, 1, 0,       0,
                      0, 0,-1,       0,
                      0, 0, 0,       1,
                    Q,1);

         // 2nd matrix - corresponding to sol[2]
         k = ind[2][3];
         set_matrix( -ctheta, -stheta, 0, -p[k].l*ctheta,
                      stheta, -ctheta, 0,  p[k].l*stheta,
                           0,       0, 1,            0,
                           0,       0, 0,            1,
                    Q,2);
      };

      i = 4;  // branching starts at atom 4
   };

   // omega angles (torsion angles)
   cosomega12 = sol[i-1].cosomega;
   if (cosomega12 == -2.0)
   {
      i3 = sol[i-1].ref3;  i2 = sol[i-1].ref2;  i1 = sol[i-1].ref1;
      cosomega12 = cosomega_coords(i3+1,i2+1,i1+1,i,sol,p,ind);
   };
   sinomega1  =  sqrt(1.0 - cosomega12*cosomega12);
   sinomega2  = -sinomega1;

   // theta angles ("bond" angles)
   ctheta = sol[i-1].costheta;
   if (ctheta == -2.0)
   {
      i2 = sol[i-1].ref2;  i1 = sol[i-1].ref1;
      ctheta = costheta_coords(i2+1,i1+1,i,sol,p,ind);
   };
   stheta = sqrt(1.0 - ctheta*ctheta);

   // distance between current atom and first reference
   k = ind[1+sol[i-1].ref1][i];  dist34 = p[k].l;

   // computing atomic coordinates
   if (op.method == 1)
   {
      // mathod1
      set_matrix(            -ctheta,            -stheta,          0,           -dist34*ctheta,
                   stheta*cosomega12, -ctheta*cosomega12, -sinomega1, dist34*stheta*cosomega12,
                    stheta*sinomega1,  -ctheta*sinomega1, cosomega12,  dist34*stheta*sinomega1,
                                   0,                  0,          0,                        1,
                 Q,i-1);

      sol[i-1].x = Q[i-1][3];  sol[i-1].y = Q[i-1][7];  sol[i-1].z = Q[i-1][11];
   }
   else
   {
      // method2
      gen_U(i,sol,U);  // only once
      gen_coords(i,sol,U,dist34,ctheta,stheta,cosomega12,sinomega1);
   };

   // left-handed branch ("+",-1)
   sol[i-1].branch = -1;

   // performing the pruning test
   pru = pruningtest(i,p,ind,sol,op,info,plde);

   if (pru == 0)
   {
      if (i < n)  
      {
         // working on next atom
         bp(i+1,n,m,p,ind,sol,op,plde,Q,info);
      }
      else
      {
         // a complete solution is found
         (*info).nsols = (*info).nsols + 1;
         nlde = 0;  lde = 0.0;
         for (k = 0; k < n; k++)  lde = lde + plde[k].val;
         for (k = 0; k < n; k++)  nlde = nlde + plde[k].n;
         if (nlde != 0)  lde = lde / (double) nlde;

         // printing the found molecule (if required)
         if (op.print > 1)  printpdb(n,sol,lde,(*info).filename,(*info).nsols);

         // best solution found so far
         if (lde < (*info).best_lde)
         {
            (*info).best_lde = lde;
            (*info).best_sol = (*info).nsols;
            if (op.print == 1)  printpdb(n,sol,lde,(*info).filename,0);
         };
      };
   };

   // if only one solution is required, bp stops as soon as the first one is found
   if (op.allone == 1)  if ((*info).nsols > 0)  goto END;

   // the search stops when maxsols solutions are found
   if ((*info).nsols >= (*info).maxsols)  goto END;

   // if the sine of omega is zero, there is actually no branching
   if (fabs(sinomega1) < epsomega)  goto END;

   // skipping solutions that belong to the 2nd symmetric half of the tree (if required)
   if (op.symmetry == 1)  if (i == 4)  goto END;

   // computing atomic coordinates
   if (op.method == 1)
   {
      // method1
      set_matrix(            -ctheta,            -stheta,          0,           -dist34*ctheta,
                   stheta*cosomega12, -ctheta*cosomega12, -sinomega2, dist34*stheta*cosomega12,
                    stheta*sinomega2,  -ctheta*sinomega2, cosomega12,  dist34*stheta*sinomega2,
                                   0,                  0,          0,                        1,
                 Q,i-1);

      sol[i-1].x = Q[i-1][3];  sol[i-1].y = Q[i-1][7];  sol[i-1].z = Q[i-1][11];
   }
   else
   {
      // method2
      gen_coords(i,sol,U,dist34,ctheta,stheta,cosomega12,sinomega2);
   };

   // right-handed branch ("-",+1)
   sol[i-1].branch = +1;

   // performing the pruning test
   pru = pruningtest(i,p,ind,sol,op,info,plde);

   if (pru == 0)
   {
      if (i < n)
      {
         // working on next atom
         bp(i+1,n,m,p,ind,sol,op,plde,Q,info);
      }
      else
      {
         // a complete solution is found
         (*info).nsols = (*info).nsols + 1;
         nlde = 0;  lde = 0.0;
         for (k = 0; k < n; k++)  lde = lde + plde[k].val;
         for (k = 0; k < n; k++)  nlde = nlde + plde[k].n;
         if (nlde != 0)  lde = lde / (double) nlde;

         // printing the found molecule (if required)
         if (op.print > 1)  printpdb(n,sol,lde,(*info).filename,(*info).nsols);

         // best solution found so far
         if (lde < (*info).best_lde)
         {
            (*info).best_lde = lde;
            (*info).best_sol = (*info).nsols;
            if (op.print == 1)  printpdb(n,sol,lde,(*info).filename,0);
         };
      };
   };

END:return (*info).nsols;
};


