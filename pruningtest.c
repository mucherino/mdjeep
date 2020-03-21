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
              Mar 21 2020  v.0.3.1  function BoxDDF added
******************************************************************************************************/

#include "bp.h"

double boxeps = 0.1;

// Direct Distance Feasibility pruning device
// -> it outputs the number of distances that are not satisfied
int DDF(int id,VERTEX *v,double **X,double eps)
{
   int count;
   double dist;
   REFERENCE *ref = v[id].ref;

   count = 0;
   while (ref != NULL)
   {
      dist = distance(ref->otherId,id,X);
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

// Box Direct Distance Feasibility pruning device
// -> it works on the boxes used in the coarse-grained representation
// -> it outputs the number of distances that are not satisfied
int BoxDDF(int id,VERTEX *v,double **lX,double **uX,double eps)
{
   int i;
   int count;
   int otherId;
   double diff;
   double min,max;
   REFERENCE *ref = v[id].ref;

   count = 0;
   while (ref != NULL)
   {
      if (isIntervalDistance(ref,eps))
      {
         min = 0.0;  max = 0.0;
         otherId = ref->otherId;
         for (i = 0; i < 3; i++)
         {
            if (uX[i][otherId] < lX[i][id])  // [lX,uX](otherId) | [lX,uX](id)
            {
               diff = lX[i][id] - uX[i][otherId];
               min = min + diff*diff;
               diff = uX[i][id] - lX[i][otherId];
               max = max + diff*diff;
            }
            else if (uX[i][id] < lX[i][otherId])  // [lX,uX](id) | [lX,uX](otherId)
            {
               diff = lX[i][otherId] - uX[i][id];
               min = min + diff*diff;
               diff = uX[i][otherId] - lX[i][id];
               max = max + diff*diff;
            }
            else  // they intersect, min distance component is 0
            {
               if (lX[i][otherId] < lX[i][id])
                  diff = lX[i][otherId];
               else
                  diff = lX[i][id];
               if (uX[i][otherId] > uX[i][id])
                  diff = uX[i][otherId] - diff;
               else
                  diff = uX[i][id] - diff;
               max = max + diff*diff;
            };
         };

         min = sqrt(min);  max = sqrt(max);
         if (max < lowerBound(ref) - boxeps || min > upperBound(ref) + boxeps)  count++;
      };
      ref = ref->next;
   };

   return count;
};

