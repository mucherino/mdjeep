/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - pruning test
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
*****************************************************************************/

#include "bp.h"

int pruningtest(int natoms,PROBL *p,int **ind,SOLUTION *sol,OPTION op,INFORMATION *info,PARTIALDE *plde)
{
   int i,j,k;
   int check;
   double diff,dist;

   j = natoms - 1;
   plde[j].n = 0;
   plde[j].val = 0.0;
   check = 0;

   for (i = 0; i < j && check == 0; i++)
   {
      k = ind[i+1][j+1];
      if (k != -1)
      {
         dist = distance(sol[i].x,sol[i].y,sol[i].z,sol[j].x,sol[j].y,sol[j].z);
         diff = dist - p[k].l;
         if (diff < 0.0)  diff = -diff;
         if (diff > op.eps)
         {
             check = 1;
         };

         if (check == 0 && p[k].l != 0.0)
         {
            plde[j].n = plde[j].n + 1;
            plde[j].val = plde[j].val + diff/p[k].l;
         };
      };
   };

   if (check == 1)
   {
      (*info).pruning = (*info).pruning + 1;
   };

   return check;
};


