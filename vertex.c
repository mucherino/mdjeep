/***********************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - vertex functions
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version for managing the vertex structure
              Mar 21 2020  v.0.3.1  functions initialClique, isDDGP, isDMDGP, findSymmetries, 
                                    and findReferences added
              May 19 2020  v.0.3.2  functions nullTriplet, isNullTriplet, isValidTriplet, cloneTriplet added
                                    the function findReferences is replaced by the functions nextTripletRef,
                                    isExactClique, findReferencesExactCase and findReferencesIntervalCase
              Apr 13 2020  v.0.3.2  patch
************************************************************************************************************/

#include "bp.h"

extern double INFTY;

// this function initializes a VERTEX structure
void initVertex(VERTEX *v,int Id,int groupId,char *Name,char *Group)
{
   v->Id = Id;  
   v->groupId = groupId;
   v->Name = strdup(Name);
   v->Group = strdup(Group);
   v->ref = NULL;
};

// given a VERTEX, this function outputs its Id
int getVertexId(VERTEX v)
{
   return v.Id;
};

// given a VERTEX, this function outputs its groupId
int getVertexGroupId(VERTEX v)
{
   return v.groupId;
};

// given a VERTEX, this function makes a copy of its name
char* getVertexName(VERTEX v)
{
   return strdup(v.Name);
};

// given a VERTEX, this function makes a copy of its group name
char* getVertexGroupName(VERTEX v)
{
   return strdup(v.Group);
};

// given a VERTEX array and two *valid* indices i and j in the array,
// this function gives the REFERENCE containing the distance between vertices i and j
// (or NULL if such a distance does not exist)
REFERENCE* getReference(VERTEX* v,int i,int j)
{
   int I,J;
   REFERENCE *ref;

   if (i < j)
   {
      I = i;  J = j;
   }
   else
   {
      I = j;  J = i;
   };

   ref = v[J].ref;
   while (ref != NULL && otherVertexId(ref) != I)  ref = ref->next;

   return ref;
};

// given an array of VERTEX, this function gives the total number of distances that are found in the several 
// REFERENCE structures
int totalNumberOfDistances(int n,VERTEX *v)
{
   int i;
   int count = 0;
   for (i = 0; i < n; i++)
   {
      count = count + numberOfDistances(v[i].ref);
   };
   return count;
};

// given an array of VERTEX, this function gives the total number of *exact* distances that are found in the 
// several REFERENCE structures
int totalNumberOfExactDistances(int n,VERTEX *v,double eps)
{
   int i;
   int count = 0;
   for (i = 0; i < n; i++)
   {  
      count = count + numberOfExactDistances(v[i].ref,eps);
   };
   return count;
};

// given an array of VERTEX, this function gives the total number of *precise* distances that are found in the
// several REFERENCE structures
int totalNumberOfPreciseDistances(int n,VERTEX *v,int ndigits)
{
   int i,count;
   double eps;

   // definition of eps 
   eps = 1.0;  for (i = 0; i < ndigits; i++)  eps = 0.1*eps;

   // counting
   count = 0;
   for (i = 0; i < n; i++)
   {
      count = count + numberOfPreciseDistances(v[i].ref,eps,ndigits);
   };

   return count;
};

// given an array of VERTEX, this function copies all distances found in the REFERENCE structures in one vector
// -> if the lower and the upper bounds are different, their average is stored in the vector
// -> the pointer to vector is given as returning argument, memory for this vector is automatically allocated
double* getDistanceList(int n,VERTEX *v)
{
   int i,k,m;
   double *y;
   REFERENCE *ref;

   // memory allocation
   m = totalNumberOfDistances(n,v);
   y = allocateVector(m);

   // loading distances in vector
   k = 0;
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         y[k] = 0.5*(lowerBound(ref) + upperBound(ref));
         ref = ref->next;
         k++;
      };
   };

   return y;
};

// given an array of VERTEX, this function verifies whether its initial 3 vertices form a 3-clique of exact distances
// -> the double value "eps" represents the tolerance to distinguish between exact and interval distances
bool initialClique(int n,VERTEX *v,double eps)
{
   REFERENCE *ref;

   // if the number of vertices is smaller than 3, then we cannot form an initial clique
   if (n < 3)  return false;

   // checking whether the first and second vertices are connected by an exact distance
   ref = getReference(v,0,1);
   if (ref == NULL)  return false;
   if (isIntervalDistance(ref,eps))  return false;

   // checking whether the first 3 vertices form a 3-clique
   ref = getReference(v,0,2);
   if (ref == NULL)  return false;
   if (isIntervalDistance(ref,eps))  return false;
   ref = getReference(v,1,2);
   if (ref == NULL)  return false;
   if (isIntervalDistance(ref,eps))  return false;

   // all tests passed, the first 3 vertices of the instance form a clique of exact distances
   return true;
};

// given an array of VERTEX, this function verifies whether the discretization assumptions (DDGP instance) are satisfied
// -> the double value "eps" represents the tolerance to distinguish between exact and interval distances
// -> the boolean "clique" indicates whether it is already known that the first three instance vertices form a clique
// -> the returning value is an integer:  0 if the instance is a DDGP; 
//                                       >0 is the first vertex rank in the instance for which the assumption is not satisfied
int isDDGP(int n,VERTEX *v,double eps,bool clique)
{
   int i;

   // checking the initial clique
   if (!clique)  if (!initialClique(n,v,eps))  return 2;

   // counting the number of exact and interval distances for all other vertices
   for (i = 3; i < n; i++)
   {
      // all distances
      if (numberOfDistances(v[i].ref) < 3)  return i;

      // exact distances
      if (numberOfExactDistances(v[i].ref,eps) < 2)  return i;
   };

   // all tests passed, this is a DDGP instance
   return 0;
};

// given an array of VERTEX, this function verifies whether the discretization assumptions (DMDGP instance) are satisfied
// (DMDGP instances are DDGP instances for which the "consecutivity assumption" is satisfied)
// -> the double value "eps" represents the tolerance to distinguish between exact and interval distances
// -> if the boolean "ddgp" is set to true, then the DDGP verification is not performed and only the cons. assumption is verified
bool isDMDGP(int n,VERTEX *v,double eps,bool ddgp)
{
   int i;
   REFERENCE *r1,*r2,*r3;

   // checking the initial clique
   if (!ddgp)  if (!initialClique(n,v,eps))  return false;

   // checking whether the instance is a DDGP
   if (!ddgp)  if (isDDGP(n,v,eps,true) != 0)  return false;

   // checking the consecutivity assumption
   for (i = 3; i < n; i++)
   {
      r3 = getReference(v,i-3,i);  if (r3 == NULL)  return false;
      r2 = getReference(v,i-2,i);  if (r2 == NULL)  return false;
      r1 = getReference(v,i-1,i);  if (r1 == NULL)  return false;
      if (isIntervalDistance(r1,eps) || isIntervalDistance(r2,eps))  return false;
   };

   // all tests passed, this is a DMDGP instance (the consecutivity assumption is satisfied)
   return true;
};

// given an array of VERTEX, this function identifies the symmetric vertices of the instance
// -> the function gives a meaningful result only when the instance satisfies the consecutivity assumption
// -> the result is stored in the boolean array "sym", which needs to be pre-allocated (size n)
void findSymmetries(int n,VERTEX *v,bool *sym)
{
   int i,j;
   REFERENCE *ref;

   // initialization
   sym[0] = false;  sym[1] = false;  sym[2] = false;
   for (j = 3; j < n; j++)  sym[j] = true;

   // checking the set of vertices every distance covers
   for (i = 0; i < n; i++)
   {
      ref = v[i].ref;
      while (ref != NULL)
      {
         for (j = otherVertexId(ref) + 4; j <= i; j++)  sym[j] = false;
         ref = ref->next;
      };
   };
};

// this function returns a "null" triplet (all references are set to NULL)
triplet nullTriplet()
{
   triplet t;
   t.r1 = NULL;
   t.r2 = NULL;
   t.r3 = NULL;
   return t;
};

// given a triplet of REFERENCEs, this functions verifies whether it is a "null" triplet or not
bool isNullTriplet(triplet t)
{
   return t.r1 == NULL || t.r2 == NULL || t.r3 == NULL;
};

// given a triplet of REFERENCEs, this function verifies whether it is a valid triplet or not
// -> a valid triplet is not a null triplet and is at most related to one interval distance
// -> the double value "eps" is the tolerance to distinguish between exact and interval distances
bool isValidTriplet(triplet t,double eps)
{
   int count = 0;
   if (isNullTriplet(t))  return false;
   if (t.r1 == t.r2)  return false;
   if (t.r2 == t.r3)  return false;
   if (t.r3 == t.r1)  return false;
   if (isExactDistance(t.r1,eps))  count++;
   if (isExactDistance(t.r2,eps))  count++;
   if (isExactDistance(t.r3,eps))  count++;
   return count >= 2;
};

// this function clones a *valid* triplet of REFERENCEs and makes sure that the only interval
// distance is always placed in r3
triplet cloneTriplet(triplet t,double eps)
{
   triplet q;

   if (isIntervalDistance(t.r1,eps))
   {
      q.r3 = t.r1;  q.r2 = t.r2;  q.r1 = t.r3;
   }
   else if (isIntervalDistance(t.r2,eps))
   {
      q.r3 = t.r2;  q.r2 = t.r1;  q.r1 = t.r3;
   }
   else
   {
      q.r3 = t.r3;  q.r2 = t.r2;  q.r1 = t.r1;
   };

   return q;
};

// given a vertex reference (referring to a distance to all the following ones), and a triplet of REFERENCEs,
// this function finds the next valid triplet of REFERENCEs
// -> a starting triplet is searched if the input t is a triplet of NULL's
// -> the returning value is the next valid triplet (if any), or a triplet of NULL's otherwise
// -> the double value "eps" is the tolerance to distinguish between exact and interval distances
triplet nextTripletRef(REFERENCE *ref,triplet t,double eps)
{
   int i;
   triplet nt;

   // null triplet (in case the next triplet does not exist)
   nt.r1 = NULL;  nt.r2 = NULL;  nt.r3 = NULL;

   // starting triplet
   if (isNullTriplet(t))
   {
      t.r3 = ref;
      t.r2 = nextDistance(ref);
      t.r1 = nextDistance(t.r2);
      if (isNullTriplet(t))  return nt;  // no triplets exist
      if (isValidTriplet(t,eps))  return t;  // found already!
   };

   // searching the next triplet
   do // t.r3
   {
      do // t.r2
      {
         do // t.r1
         {
            t.r1 = nextDistance(t.r1);
            if (isValidTriplet(t,eps))  return t;
         }
         while (t.r1 != NULL);

         t.r2 = nextDistance(t.r2);
         t.r1 = nextDistance(t.r2);
         if (isValidTriplet(t,eps))  return t;
      }
      while (t.r2 != NULL);

      t.r3 = nextDistance(t.r3);
      t.r2 = nextDistance(t.r3);
      t.r1 = nextDistance(t.r2);
      if (isValidTriplet(t,eps))  return t;
   }
   while (t.r3 != NULL);

   // no valid triplet found
   return nt;
};

// given a vertex (via its rank) and its entire VERTEX array, and given a triplet of valid reference vertices for id,
// this function verifies whether the set of reference vertices from a clique of exact vertices or not
// -> if yes, the value of cosangle is the cosine of the angle formed by the triplet, it is unchanged otherwise
bool isExactClique(int id,VERTEX *v,triplet t,double eps,double *cosangle)
{
   int i,j,k;
   double dij,djk,dik;
   REFERENCE *ij,*jk,*ik;

   // reference ranks
   i = otherVertexId(t.r3);
   j = otherVertexId(t.r2);
   k = otherVertexId(t.r1);

   // is this a clique of exact distances?
   ij = getReference(v,i,j);
   if (ij == NULL)  return false;
   if (!isExactDistance(ij,eps))  return false;
   jk = getReference(v,j,k);
   if (jk == NULL)  return false;
   if (!isExactDistance(jk,eps))  return false;
   ik = getReference(v,i,k);
   if (ik == NULL)  return false;
   if (!isExactDistance(ik,eps))  return false;

   // computing the angle formed by the triplet (0 or 180 => flat triplet)
   dij = lowerBound(ij);  djk = lowerBound(jk);  dik = lowerBound(ik);
   (*cosangle) = dij*dij + djk*djk - dik*dik;
   (*cosangle) = (*cosangle) / (2.0*dij*djk);

   // the triplet is a clique of exact distances
   return true;
};

// given an array of VERTEX, this function creates a triplet structure containing its three reference 
// vertices (exact case)
// -> it is supposed that all reference vertices for id are exact distances
// -> the triplet that form an angle "far from" a multiple of pi/2
// -> the cosine of angle formed by the triplet is given in output
// -> if the function returns a null triplet, then the discretization assumptions may not be satisfied
triplet findReferencesExactCase(int id,VERTEX *v,double eps,double *cosine)
{
   double best,cosangle;
   triplet t,refs;

   // computing the initial triplet
   refs = nullTriplet();
   t = nextTripletRef(v[id].ref,nullTriplet(),eps);
   if (isNullTriplet(t))  return refs;

   // verifying all other triplets and computing the angle (its cosine) they form
   best = -2.0;
   do
   {
      if (isValidTriplet(t,eps))
      {
         cosangle = 0.0;
         if (isExactClique(id,v,t,eps,&cosangle))
         {
            if (isNullTriplet(refs))
            {
               best = cosangle;
               if (cosangle < 0.0)  best = -cosangle;
               refs = cloneTriplet(t,eps);
            }
            else
            {
               if (cosangle < 0.0)  cosangle = -cosangle;  // in [0,1]
               if (best > cosangle)  // we want it to be as close as possible to 0, so we minimize
               {
                  best = cosangle;
                  refs = cloneTriplet(t,eps);
               };
            };
         };
      };
      t = nextTripletRef(v[id].ref,t,eps);
   }
   while (!isNullTriplet(t));

   // returning the triplets (and its cosine)
   if (isExactClique(id,v,refs,eps,cosine))  return refs;
   return nullTriplet();
};

// given an array of VERTEX, this function creates a triplet structure containing its three reference 
// vertices (case where one distance may be represented by an interval)
// -> at least 2 references are related to exact distances
// -> the triplet containing the interval distance that is the smallest in range is selected
// -> if the function returns a null triplet, then the discretization assumptions may not be satisfied
triplet findReferencesIntervalCase(int id,VERTEX *v,double eps)
{
   double minrange,range;
   triplet t,refs;

   // computing the initial triplet
   refs = nullTriplet();
   t = nextTripletRef(v[id].ref,nullTriplet(),eps);
   if (isNullTriplet(t))  return refs;

   // verifying all other triplets and computing the range of the interval distance
   minrange = INFTY;
   do
   {
      if (isValidTriplet(t,eps))
      {
         if (isNullTriplet(refs))
         {
            refs = cloneTriplet(t,eps);
            minrange = maximum(rangeOfDistance(t.r1),rangeOfDistance(t.r2),rangeOfDistance(t.r3));
         }
         else
         {
            range = maximum(rangeOfDistance(t.r1),rangeOfDistance(t.r2),rangeOfDistance(t.r3));
            if (range < minrange)
            {
               refs = cloneTriplet(t,eps);
               minrange = range;
            };
        };
     };
     t = nextTripletRef(v[id].ref,t,eps);
  }
  while (!isNullTriplet(t));

  // returning the triplet
  return refs;
};

// given an array of VERTEX, this function prints all distance found in the several REFERENCE structures (stdout)
void printDistanceList(int n,VERTEX *v,bool symmetric)
{
   int i,j;
   REFERENCE *ref;

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         if (symmetric || i < j)
         {
            ref = getReference(v,i,j);
            if (ref != NULL)
            {
               fprintf(stderr," %3d %3d  [%10.6lf,%10.6lf]\n",i,j,lowerBound(ref),upperBound(ref));
            };
         };
      };
   };
};

// given a VERTEX, this function prints its main information (stdout)
void printVertex(VERTEX v)
{
   printf("[%d,%d,%s,%s] (%d distances)\n",v.Id,v.groupId,v.Name,v.Group,numberOfDistances(v.ref));
};

// this function properly frees the memory for a VERTEX array
VERTEX* freeVertex(int n,VERTEX *v)
{
   int i;
   for (i = 0; i < n; i++)  if (v[i].ref != NULL)  freeReference(v[i].ref);
   free(v);
   return NULL;
};

