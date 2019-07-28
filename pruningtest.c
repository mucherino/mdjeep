/*****************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - pruning devices
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    no changes w.r.t previous version
              Jul 28 2019  v.0.3.0  the function DDF is adapted to the new data structures
                                    LDE and MDE values are now computed in another functions (objfun.c)
******************************************************************************************************/

#include "bp.h"

// Direct Distance Feasibility pruning device
// it outputs the number of distances that are not satisfied
int DDF(int i,VERTEX *v,double **X,double eps)
{
   int count;
   double dist;
   REFERENCE *ref = v[i].ref;

   count = 0;
   while (ref != NULL)
   {
      dist = distance(ref->otherId,i,X);
      if (isExactDistance(ref,eps))
      {
         if (fabs(lowerBound(ref) - dist) > eps)  count++;
      }
      else
      {
         if (dist < lowerBound(ref) - eps)  count++;
         if (dist > upperBound(ref) + eps)  count++;
      }
      ref = ref->next;
   };
   return count;
};

