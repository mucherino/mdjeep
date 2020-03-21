/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - distance functions
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version for working with the newly introduced 
                                    data structures
              Mar 21 2020  v.0.3.1  adding numberOfExactDistances and rangeOfDistance
****************************************************************************************************/

#include "bp.h"

// this function computes the distance between two sets of coordinates in 3D
double pairwise_distance(double xA,double yA,double zA,double xB,double yB,double zB)
{
   double value;
   double dx,dy,dz;
   dx = xB - xA;
   dy = yB - yA;
   dz = zB - zA;
   value = (dx*dx) + (dy*dy) + (dz*dz);
   return sqrt(value);
};

// this function computes the distance between two columns of X (dimension is set to 3)
double distance(int i,int j,double **X)
{
   double value;
   double dx,dy,dz;
   dx = X[0][i] - X[0][j];
   dy = X[1][i] - X[1][j];
   dz = X[2][i] - X[2][j];
   value = (dx*dx) + (dy*dy) + (dz*dz);
   return sqrt(value);
}

// this function initializes a REFERENCE structure with the first distance
REFERENCE* initReference(int otherId,double lb,double ub)
{
   REFERENCE *ref = (REFERENCE*)calloc(1,sizeof(REFERENCE));
   ref->otherId = otherId;
   ref->lb = lb;
   ref->ub = ub;
   ref->next = NULL;
   return ref;
};

// this function adds a new distance to an already-initialized REFERENCE structure
// it returns the created REFERENCE
REFERENCE* addDistance(REFERENCE *ref,int otherId,double lb,double ub)
{
   while (ref->next != NULL)  ref = ref->next;
   ref->next = (REFERENCE*)calloc(1,sizeof(REFERENCE));
   ref->next->otherId = otherId;
   ref->next->lb = lb;
   ref->next->ub = ub;
   ref->next->next = NULL;
   return ref->next;
};

// given a *valid* REFERENCE, this function outputs the index of the other vertex
int otherVertexId(REFERENCE *ref)
{
   return ref->otherId;
};

// given a *valid* REFERENCE, this function outputs the distance lower bound
double lowerBound(REFERENCE* ref)
{
   return ref->lb;
};

// given a *valid* REFERENCE, this function outputs the distance upper bound
double upperBound(REFERENCE* ref)
{
   return ref->ub;
};

// given a REFERENCE, this function outputs the total number of referenced distances 
// (including the current one)
int numberOfDistances(REFERENCE *ref)
{
   int count = 0;
   if (ref != NULL)
   {
      count++;
      while (ref->next != NULL)
      {
         count++;
         ref = ref->next;
      };
   };
   return count;
};

// given a REFERENCE, this function outputs the total number of exact referenced distances (including the current one)
// the double value "eps" is the tolerance to distinguish between exact and interval distances
int numberOfExactDistances(REFERENCE *ref,double eps)
{
   int count = 0;
   if (ref != NULL)
   {
      if (isExactDistance(ref,eps))  count++;
      while (ref->next != NULL)
      {
         if (isExactDistance(ref->next,eps))  count++;
         ref = ref->next;
      };
   };
   return count;
};

// given a REFERENCE, this function outputs the range of the corresponding distance
// -> it corresponds to 0 if this is a exact reference
// -> the function outputs 0 also when the REFERENCE is NULL
double rangeOfDistance(REFERENCE *ref)
{
   if (ref != NULL)
      return upperBound(ref) - lowerBound(ref);
   return 0;
};

// this function verifies whether the given reference corresponds to an "exact" distance
// (with tolerance eps)
bool isExactDistance(REFERENCE *ref,double eps)
{
   if (ref != NULL)
      return upperBound(ref) - lowerBound(ref) <= eps;
   return false;
};

// given a REFERENCE, this function outputs the next REFERENCE related to an "exact" distance
// it gives NULL if the input REFERENCE is NULL or when no such a distance exists
REFERENCE* nextExactDistance(REFERENCE *current,double eps)
{
   REFERENCE *ref = current;
   if (ref != NULL)
   {
      while (ref->next != NULL)
      {
         if (upperBound(ref->next) - lowerBound(ref->next) <= eps)  return ref->next;
         ref = ref->next;
      };
      return ref->next;
   };
   return NULL;
};

// this function verifies whether the given reference corresponds to an "interval" distance
// (with tolerance eps)
bool isIntervalDistance(REFERENCE *ref,double eps)
{
   if (ref != NULL)
      return upperBound(ref) - lowerBound(ref) > eps;
   return false;
};

// given a REFERENCE, this function outputs the next REFERENCE related to an "interval" distance
// it gives NULL if the input REFERENCE is NULL or when no such a distance exists
REFERENCE* nextIntervalDistance(REFERENCE *current,double eps)
{
   REFERENCE *ref = current;
   if (ref != NULL)
   {
      while (ref->next != NULL)
      {
         if (upperBound(ref->next) - lowerBound(ref->next) > eps)  return ref->next;
         ref = ref->next;
      };
      return ref->next;
   };
   return NULL;
};

// given a REFERENCE, all its distances are printed on the screen (stdout)
void printDistances(REFERENCE *ref)
{
   while (ref != NULL)
   {
      printf("%3d) [%10.7lf,%10.7lf]\n",ref->otherId,ref->lb,ref->ub);
      ref = ref->next;
   };
};

// given a REFERENCE, this function frees its memory
REFERENCE* freeReference(REFERENCE *ref)
{
   while (ref != NULL)
   {
      if (ref->next == NULL)
      {
         free(ref);
         ref = NULL;
      }
      else
      {
         while (ref->next->next != NULL)  ref = ref->next;
         free(ref->next);
         ref->next = NULL;
      };
   };
   return NULL;
};

