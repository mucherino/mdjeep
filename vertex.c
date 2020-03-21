/********************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - vertex functions
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version for managing the vertex structure
              Mar 21 2020  v.0.3.1  functions initialClique, isDDGP, isDMDGP, findSymmetries, 
                                    and findReferences added
*********************************************************************************************************/

#include "bp.h"

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
   while (ref != NULL && ref->otherId != I)  ref = ref->next;

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

// given an array of VERTEX, this function creates a triplet structure containing its three reference vertices
// -> at least two selected reference vertices are exact distances
// -> when more than one triplet of exact distances is available, the one forming a theta angle "farther" from a
//    a multiple of pi is selected
// -> when more than one interval distance may be chosen for the triplet, then the one having the smallest range
//    is selected
triplet findReferences(int n,VERTEX *v,int id,double eps,bool onlyExact)
{
   int i,j,k;
   int m,mexact;
   double pimultiple,lesspimult;
   double dij,dik,djk;
   REFERENCE *ref;
   REFERENCE *rij,*rjk,*rik;
   REFERENCE **exact;
   REFERENCE **interval;
   triplet t;

   // initializing the triplet
   t.r1 = NULL;  t.r2 = NULL;  t.r3 = NULL;

   // checking the number of expected references
   ref = v[id].ref;
   m = numberOfDistances(ref);
   mexact = numberOfExactDistances(ref,eps);

   // are there enough references?
   if (m < 3)  return t;
   if (mexact < 2)  return t;

   // memory allocation
   exact = (REFERENCE**)calloc(mexact,sizeof(REFERENCE*));
   interval = (REFERENCE**)calloc(m-mexact,sizeof(REFERENCE*));

   // loading the references separated in the two allocated arrays
   // (the order is from the farthest to the nearest reference vertex)
   i = 0;  j = 0;
   while (ref != NULL)
   {
      if (isExactDistance(ref,eps))
      {
         exact[i] = ref;
         i++;
      }
      else
      {
         interval[j] = ref;
         j++;
      };
      ref = ref->next;
   };

   // selecting the best triplet of references
   if (mexact >= 3)
   {
      // initial option with exact distances (references nearest in terms of rank)
      t.r1 = exact[0];
      for (i = 0; i < mexact; i++)  if (t.r1->otherId < exact[i]->otherId)  t.r1 = exact[i];
      t.r2 = NULL;
      for (i = 0; i < mexact; i++)
      {
         if (t.r2 == NULL)
         {
            if (exact[i] != t.r1)  t.r2 = exact[i];
         }
         else
         {
            if (t.r2->otherId < exact[i]->otherId && exact[i] != t.r1)  t.r2 = exact[i];
         };
      };
      t.r3 = NULL;
      for (i = 0; i < mexact; i++)
      {
         if (t.r3 == NULL)
         {
            if (exact[i] != t.r1 && exact[i] != t.r2)  t.r3 = exact[i];
         }
         else
         {
            if (t.r1->otherId < exact[i]->otherId && exact[i] != t.r1 && exact[i] != t.r2)  t.r3 = exact[i];
         };
      };

      // exploring other possible options with exact distances
      if (!onlyExact && mexact > 3)
      {
         // we choose the triplet of exact references that form the theta angle
         // that is the farthest from a multiple of pi
         i = mexact - 3;  j = mexact - 2;  k = mexact - 1;
         lesspimult = 2.0;
         while (k >= 3)
         {
            // evaluating the alignment of the exact references
            rij = getReference(v,i,j);
            if (rij != NULL)  if (isExactDistance(rij,eps))
            {
               rjk = getReference(v,j,k);
               if (rjk != NULL)  if (isExactDistance(rjk,eps))
               {
                  rik = getReference(v,i,k);
                  if (rik != NULL)  if (isExactDistance(rik,eps))
                  {
                     dij = lowerBound(rij);  djk = lowerBound(rjk);  dik = lowerBound(rik);
                     if (dij != 0.0 && djk != 0.0 && dik != 0.0)
                     {
                        pimultiple = dij*dij + djk*djk - dik*dik;
                        pimultiple = pimultiple / (2.0*dij*djk);
                        if (pimultiple < 0.0)  pimultiple = -pimultiple;
                        if (fabs(pimultiple) < eps)
                           pimultiple = eps - fabs(pimultiple);
                        else if (fabs(pimultiple - 1.0) < eps)
                           pimultiple = eps - fabs(pimultiple - 1.0);
                        else
                           pimultiple = 0.0;
                        if (pimultiple < lesspimult)
                        {
                           lesspimult = pimultiple;
                           t.r3 = exact[i];
                           t.r2 = exact[j];
                           t.r1 = exact[k];
                        };
                     };
                  };
               };
            };

            // next possible triplet (in reverse order)
            i--;
            if (i == -1)
            {
               j--;
               i = j - 1;
            };
            if (j == 0)
            {
               k--;
               j = k - 1;
               i = j - 1;
            };
         };
      };
   }
   else if (!onlyExact)
   {
      // an interval is used for one of the distances
      // (we select the one with smallest range)
      k = 0;
      for (i = 1; i < m - mexact; i++)
      {
         if (rangeOfDistance(interval[k]) > rangeOfDistance(interval[i]))  k = i;
      };

      // defining the triplet
      t.r3 = interval[k];
      t.r2 = exact[1];
      t.r1 = exact[0];
   };

   // ending
   free(interval);
   free(exact);
   return t;
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

