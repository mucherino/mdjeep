/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - function BP
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
*****************************************************************************/

#include "bp.h"

int bp(int i,int n,int m,PROBL *p,int **ind,SOLUTION *sol,OPTION op,PARTIALDE *plde,double **Q,INFORMATION *info)
{
   int k;
   int nlde,pru;
   double dist34,lde;
   double ctheta,stheta;
   double cosomega12,sinomega1,sinomega2;

   // is this the first atom?
   if (i == 1)
   {
      // the first three atoms can be positioned by using the available information
      // the matrices Q for these atoms are also computed because needed for the following atoms

      // first atom
      sol[0].x = 0.0;  sol[0].y = 0.0;  sol[0].z = 0.0;

      // second atom
      k = ind[1][2];
      sol[1].x = -p[k].l;  sol[1].y = 0.0;  sol[1].z = 0.0;
      set_matrix( -1, 0, 0, -p[k].l, 
                   0, 1, 0,       0,
                   0, 0,-1,       0,
                   0, 0, 0,       1,
                 Q,1);

      // third atom
      sol[2].x = -p[k].l;
      k = ind[2][3];
      ctheta = sol[2].costheta;
      stheta = sqrt(1.0 - ctheta*ctheta);

      sol[2].x = sol[2].x + p[k].l*ctheta;  sol[2].y = p[k].l*stheta;  sol[2].z = 0.0;
      set_matrix( -ctheta, -stheta, 0, -p[k].l*ctheta,
                   stheta, -ctheta, 0,  p[k].l*stheta,
                        0,       0, 1,            0,
                        0,       0, 0,            1,
                 Q,2);

      i = 4;  // branching starts at atom 4
   };

   // omega angles (torsion angles)
   cosomega12 =  sol[i-1].cosomega;
   sinomega1  =  sqrt(1.0 - cosomega12*cosomega12);
   sinomega2  = -sinomega1;

   // distances and angles needed for computing the matrix B
   k = ind[i-1][i];  dist34 = p[k].l;
   ctheta = sol[i-1].costheta;
   stheta = sqrt(1.0 - ctheta*ctheta);

   // computing the matrix Q (cosomega,+sinomega)
   set_matrix(            -ctheta,            -stheta,          0,           -dist34*ctheta,
                stheta*cosomega12, -ctheta*cosomega12, -sinomega1, dist34*stheta*cosomega12,
                 stheta*sinomega1,  -ctheta*sinomega1, cosomega12,  dist34*stheta*sinomega1,
                                0,                  0,          0,                        1,
              Q,i-1);

   // computing the atom coordinates
   sol[i-1].x = Q[i-1][3];  sol[i-1].y = Q[i-1][7];  sol[i-1].z = Q[i-1][11];  sol[i-1].branch = -1;

   // performing the pruning test
   pru = pruningtest(i,p,ind,sol,op,info,plde);

   if (pru == 0)
   {
      if (i < n)  
      {
         // working on the following atom
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

   // if the sine of omega is zero, there is actually no branching
   if (fabs(sinomega1) < 1.e-6)  goto END;

   // removing solutions that can be obtained by symmetry (if required)
   if (op.symmetry == 1)  if (i == 4)  goto END;

   // computing the matrix Q (cosomega,-sinomega)
   set_matrix(            -ctheta,            -stheta,          0,           -dist34*ctheta,
                stheta*cosomega12, -ctheta*cosomega12, -sinomega2, dist34*stheta*cosomega12,
                 stheta*sinomega2,  -ctheta*sinomega2, cosomega12,  dist34*stheta*sinomega2,
                                0,                  0,          0,                        1,
              Q,i-1);

   // computing the atom coordinates
   sol[i-1].x = Q[i-1][3];  sol[i-1].y = Q[i-1][7];  sol[i-1].z = Q[i-1][11];  sol[i-1].branch = 1;

   // performing the pruning test
   pru = pruningtest(i,p,ind,sol,op,info,plde);

   if (pru == 0)
   {
      if (i < n)
      {
         // working on the following atom
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


