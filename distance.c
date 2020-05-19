/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - distance functions
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version for working with the newly introduced 
                                    data structures
              Mar 21 2020  v.0.3.1  adding numberOfExactDistances and rangeOfDistance
              May 19 2020  v.0.3.2  adding box_distance and nextDistance
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
};

// this function computes the minimal and maximal distance between two boxes
// -> [lX,uX] is a list of boxes in dimension 3
// -> the returning value is the minimal distance between the two selected boxes (indices i and j)
// -> the maximal distance between the two boxes is returned via the double pointer m
// -> if the boxes are two singletons, then the min and max distance coincide
double box_distance(int i,int j,double **lX,double **uX,double *m)
{
   int k;
   double diff;
   double min,max;

   min = 0.0;  max = 0.0;
   for (k = 0; k < 3; k++)
   {
      if (uX[k][i] < lX[k][j])  // [lX,uX](i) | [lX,uX](j)
      {
         diff = lX[k][j] - uX[k][i];
         min = min + diff*diff;
         diff = uX[k][j] - lX[k][i];
         max = max + diff*diff;
      }
      else if (uX[k][j] < lX[k][i])  // [lX,uX](j) | [lX,uX](i)
      {
         diff = lX[k][i] - uX[k][j];
         min = min + diff*diff;
         diff = uX[k][i] - lX[k][j];
         max = max + diff*diff;
      }
      else  // they intersect: min distance component is 0
      {
         if (lX[k][i] < lX[k][j])
            diff = lX[k][i];
         else
            diff = lX[k][j];
         if (uX[k][i] > uX[k][j])
            diff = uX[k][i] - diff;
         else
            diff = uX[k][j] - diff;
         max = max + diff*diff;
      };
   };

   // ending
   (*m) = sqrt(max);
   return sqrt(min);
};

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
// -> it returns the created REFERENCE
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

// given a *valid* REFERENCE, this function outputs the index of vertex it makes reference to
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
// -> the double value "eps" is the tolerance to distinguish between exact and interval distances
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

// given a REFERENCE, this function outputs the total number of "precise" referenced distances (including the current one)
// -> the precision is evaluated by:
//    - the tolerance eps to distinguish between exact and interval distances
//    - the number ndigits of decimal digits used in the representation (min 0 / max 16)
int numberOfPreciseDistances(REFERENCE *ref,double eps,int ndigits)
{
   int count = 0;
   if (ref != NULL)
   {
      if (isExactDistance(ref,eps))  if (precisionOf(ref->lb) >= ndigits)  count++;
      while (ref->next != NULL)
      {
         if (isExactDistance(ref->next,eps))  if (precisionOf(ref->next->lb) >= ndigits)  count++;
         ref = ref->next;
      };
   };
   return count;
};

// given a REFERENCE, this function verifies whether all the corresponding distances are precise
// -> the precision is given through the number of digits used in the representation of the real numbers
bool onlyPreciseDistances(REFERENCE *ref,int ndigits)
{
   int i;
   double eps;

   // definition of eps 
   eps = 1.0;  for (i = 0; i < ndigits; i++)  eps = 0.1*eps;

   // performing the verification
   if (ref != NULL)
   {
      if (isIntervalDistance(ref,eps))  return false;
      if (precisionOf(ref->lb) < ndigits)  return false;
      while (ref->next != NULL)
      {
         if (isIntervalDistance(ref->next,eps))  return false;
         if (precisionOf(ref->next->lb) < ndigits)  return false;
         ref = ref->next;
      };
   };

   // all tests passed, all distances are precise
   return true;
};

// given a REFERENCE, this function outputs the range of the corresponding distance
// -> the function outputs 0 also when the REFERENCE is NULL
double rangeOfDistance(REFERENCE *ref)
{
   if (ref != NULL)
      return upperBound(ref) - lowerBound(ref);
   return 0;
};

// given a REFERENCE, this function outputs its next REFERENCE
// -> it gives NULL if the input REFERENCE is NULL or when no such a distance exists
REFERENCE* nextDistance(REFERENCE *current)
{
   REFERENCE* next = NULL;
   if (current != NULL)  next = current->next;
   return next;
};    

// this function verifies whether the input reference corresponds to an "exact" distance
// (with tolerance eps)
bool isExactDistance(REFERENCE *ref,double eps)
{
   if (ref != NULL)
      return upperBound(ref) - lowerBound(ref) <= eps;
   return false;
};

// given a REFERENCE, this function outputs the next REFERENCE related to an "exact" distance
// -> it gives NULL if the input REFERENCE is NULL or when no such a distance exists
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

// this function verifies whether the input reference corresponds to an "interval" distance
// (with tolerance eps)
bool isIntervalDistance(REFERENCE *ref,double eps)
{
   if (ref != NULL)
      return upperBound(ref) - lowerBound(ref) > eps;
   return false;
};

// given a REFERENCE, this function outputs the next REFERENCE related to an "interval" distance
// -> it gives NULL if the input REFERENCE is NULL or when no such a distance exists
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

