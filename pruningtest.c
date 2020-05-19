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
              May 19 2020  v 0.3.2  DDF and BoxDDF now output the partial error
                                    BoxDDF uses the function box_distance (distance.c)
******************************************************************************************************/

#include "bp.h"

// Direct Distance Feasibility pruning device
// -> id is the vertex id for which it is necessary to verify the reference distances
// -> v is the set of VERTEX structures (with size > id), and X is the current conformation
// -> DDF outputs the partial error on the entire set of reference distances related to the vertex id
//   (the partial error is computed as the sum of the MDE terms related to this id)
double DDF(int id,VERTEX *v,double **X)
{
   int n;
   double error,dist,diff;
   REFERENCE *ref = v[id].ref;

   // collecting distances and verifying error
   n = 0;  error = 0.0;
   while (ref != NULL)
   {
      n++;
      dist = distance(otherVertexId(ref),id,X);
      diff = lowerBound(ref) - dist;  if (diff > 0.0)  error = error + diff;  // only one of the two
      diff = dist - upperBound(ref);  if (diff > 0.0)  error = error + diff;  // per time can be true
      ref = ref->next;
   };

   // normalizing over the number of reference distances
   if (n != 0)  error = error/n;

   return error;
};

// Box Direct Distance Feasibility pruning device
// -> id is the vertex id whose box needs to be verified for feasibility
// -> v is the set of VERTEX structures (with size > id), and [lX,uX] is the set of boxes up to vertex id
// -> DDF outputs the partial error on the entire set of reference distances related to the vertex id
//   (the function box_distance is used to compute distances between pairs of boxes)
double BoxDDF(int id,VERTEX *v,double **lX,double **uX)
{
   int n;
   int otherId;
   double error,diff;
   double min,max;
   REFERENCE *ref = v[id].ref;

   // collecting distances and verifying error
   n = 0;  error = 0.0;
   while (ref != NULL)
   {
      otherId = otherVertexId(ref);
      min = box_distance(id,otherId,lX,uX,&max);
      diff = lowerBound(ref) - max;  if (diff > 0.0)  error = error + diff;  // only one of the two
      diff = min - upperBound(ref);  if (diff > 0.0)  error = error + diff;  // per time can be true
      ref = ref->next;
   };

   // normalizing over the number of reference distances
   if (n != 0)  error = error/n;

   return error;
};

