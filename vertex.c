/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - vertex functions
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version for managing the vertex structure
****************************************************************************************************/

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

// given an array of VERTEX, this function copies all distances found in the REFERENCE structures in one vector
// (the pointer to vector is given as returning argument, memory for this vector will be allocated)
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

// given an array of VERTEX, this function prints all distance found in the several REFERENCE structures (stdout)
void printDistanceList(int n,VERTEX *v)
{
   int i,j;
   REFERENCE *ref;

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < n; j++)
      {
         ref = getReference(v,i,j);
         if (ref != NULL)
         {
            fprintf(stderr," %3d %3d  [%10.6lf,%10.6lf]\n",i,j,lowerBound(ref),upperBound(ref));
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

