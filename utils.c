/****************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - utilities
  Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 08 2014  v.0.2    functions costheta and cosomega updated
              Jul 28 2019  v.0.3.0  functions for computing cosines (omega,theta) from coords added
                                    functions for managing omega lists added
                                    functions projection and removExtension added
              Mar 21 2020  v.0.3.1  functions numberOfOmegaIntervals, expandBounds, numberOfDigits, 
                                    minimun and maximum added
*****************************************************************************************************/

#include "bp.h"

extern double INFTY;

/* functions to manage omega angle lists */

// this function initializes an omega list
omegaList initOmegaList(double l,double u)
{
   omegaList L;
   L.first = (Omega*)calloc(1,sizeof(Omega));
   if (l < u)
   {
      L.first->l = l;
      L.first->u = u;
   }
   else
   {
      L.first->l = u;
      L.first->u = l;
   };
   L.first->prev = NULL;
   L.first->next = NULL;
   return L;
};

// this function gives the first omega interval in the list L
Omega* firstOmegaInterval(omegaList L)
{
   return L.first;
};

// this function gives the last omega interval in the list L
Omega* lastOmegaInterval(omegaList L)
{
   Omega *current;
   if (L.first != NULL)
   {
      current = L.first;
      while (current->next != NULL)  current = current->next;
      return current;
   };
   return NULL;
};

// this function gives the lower bound of an omega interval
double omegaIntervalLowerBound(Omega *current)
{
   return current->l;
};

// this function gives the upper bound of an omega interval
double omegaIntervalUpperBound(Omega *current)
{
   return current->u;
};

// this function indicates whether another omega interval follows the current one
bool omegaIntervalHasNext(Omega *current)
{
   return current->next != NULL;
};

// this function gives the pointer to the next omega interval (NULL if it doesnt exist)
Omega* omegaIntervalNext(Omega *current)
{
   return current->next;
};

// this function indicates whether another omage interval precedes the current one
bool omegaIntervalHasPrev(Omega *current)
{
   return current->prev != NULL;
};

// this function gives the pointer to the previous omega interval (NULL if it doesnt exist)
Omega* omegaIntervalPrev(Omega *current)
{
   return current->prev;
};

// this function indicates whether another omega interval exists in the given direction
bool omegaIntervalHasNextAlongDirection(Omega *current,bool asNext)
{
   bool answer = false;
   if (asNext)
      answer = omegaIntervalHasNext(current);
   else
      answer = omegaIntervalHasPrev(current);
   return answer;
};

// this function gives the pointer to the previous or next omega interval 
// (it depends on the direction, it's NULL if it doesnt exist)
Omega* omegaIntervalNextAlongDirection(Omega *current,bool asNext)
{
   Omega *new;
   if (asNext)
      new = omegaIntervalNext(current);
   else
      new = omegaIntervalPrev(current);
   return new;
};

// this function "attaches" a new omega interval to a given omega interval
// -> it is supposed that the next omega interval is NULL (otherwise memory for next is deallocated)
void attachNewOmegaInterval(Omega *current,double l,double u)
{
   if (current != NULL)
   {
      if (current->next != NULL)  freeNextOmega(current);
      current->next = (Omega*)calloc(1,sizeof(Omega));
      if (l < u)
      {
         current->next->l = l;
         current->next->u = u;
      }
      else
      {
         current->next->l = u;
         current->next->u = l;
      };
      current->next->prev = current;
      current->next->next = NULL;
   };
};

// this function splits the omega interval list starting with "current" in subintervals if 
// their "arclength" (radius*angle) is larger than the given threshold "resolution"
void splitOmegaIntervals(Omega *current,double radius,double resolution)
{
   int i;
   int div;
   double range,l,u;
   double arclength,diff;
   Omega *next;

   if (current != NULL)
   {
      do
      {
         l = current->l;
         u = current->u;
         diff = u - l;
         arclength = radius*diff;
         div = 1;
         while (arclength/div > resolution)  div++;
         if (div > 1)
         {
            range = diff/div;
            current->u = current->l + range;
            next = current->next;
            current->next = NULL;
            for (i = 2; i <= div; i++)
            {
               attachNewOmegaInterval(current,current->u,current->u+range);
               current = current->next;
            };
            if (next != NULL)
            {
               current->next = next;
               current->next->prev = current;
            };
         };
         current = current->next;
      }
      while (current != NULL);
   };
};

// this function counts the number of omega interval from a given interval
int numberOfOmegaIntervals(Omega *current)
{
   int count = 0;
   if (current != NULL)
   {
      count = 1;
      while (current->next != NULL)
      {
         count++;
         current = current->next;
      };
   };
   return count;
};

// this function prints an interval list L (stdout)
void printOmegaList(omegaList L)
{
   int i;
   Omega *current = L.first;
   if (current == NULL)
   {
      printf("[empty Omega list]\n");
   }
   else
   {
      i = 0;
      while (current != NULL)
      {
         i++;
         printf("%3d) [%10.7lf,%10.7lf]\n",i,current->l,current->u);
         current = current->next;
      };
   };
};

// this function frees the next omega interval in the list L
void freeNextOmega(Omega *a)
{
   if (a->next != NULL)  freeNextOmega(a->next);
   free(a->next);
   a->next = NULL;
};

// this function frees the list L
omegaList freeOmegaList(omegaList L)
{
   if (L.first != NULL)
   {
      if (L.first->next != NULL)  freeNextOmega(L.first);
      free(L.first);
      L.first = NULL;
   };
   return L;
};

/* functions for computing the cosine and sine of angles */

// computing the cosine of theta by using the computed positions for the references
double costheta(int i,int j,int k,VERTEX *v,double **X)
{
   REFERENCE *r12,*r23,*r13;
   double d12,d23,d13;
   double val;

   r12 = getReference(v,i,j);
   if (r12 != NULL)  d12 = lowerBound(r12);  else  d12 = distance(i,j,X);

   r23 = getReference(v,j,k);
   if (r23 != NULL)  d23 = lowerBound(r23);  else  d23 = distance(j,k,X);

   r13 = getReference(v,i,k);
   if (r13 != NULL)  d13 = lowerBound(r13);  else  d13 = distance(i,k,X);

   val = d12*d12 + d23*d23 - d13*d13;
   val = val / (2.0*d12*d23);
   if (val < -1.0)  val = -1.0;
   if (val >  1.0)  val =  1.0;
   return val;
};

// computing the cosine of omega (with available distances when possible, or computed distances)
double cosomega(int i3,int i2,int i1,int i,VERTEX *v,double **X,double range)
{
   REFERENCE *r12,*r13,*r14,*r23,*r24,*r34;
   double d12,d13,d14,d23,d24,d34;
   double d12q,d13q,d14q,d23q,d24q,d34q;
   double a,b,c,e,f;
   double val;

   r12 = getReference(v,i3,i2);  r13 = getReference(v,i3,i1);  r23 = getReference(v,i2,i1);  // may not be available
   r14 = getReference(v,i3,i);   r24 = getReference(v,i2,i);   r34 = getReference(v,i1,i);   // known because of discretization

   if (r14 == NULL || r24 == NULL || r34 == NULL)
   {
      fprintf(stderr,"cosomega: internal error; it looks like the discretization assumptions are not satisfied\n");
      abort();
   };
   d14 = lowerBound(r14) + range*(upperBound(r14) - lowerBound(r14));  d14q = d14*d14;
   d24 = lowerBound(r24);  d24q = d24*d24;
   d34 = lowerBound(r34);  d34q = d34*d34;

   if (r12 != NULL)  d12 = lowerBound(r12);  else  d12 = distance(i3,i2,X);
   d12q = d12*d12;
   if (r13 != NULL)  d13 = lowerBound(r13);  else  d13 = distance(i3,i1,X);
   d13q = d13*d13;
   if (r23 != NULL)  d23 = lowerBound(r23);  else  d23 = distance(i2,i1,X);
   d23q = d23*d23;

   a = d12q + d24q - d14q;  a = a / (2.0*d12*d24);
   b = d24q + d23q - d34q;  b = b / (2.0*d24*d23);
   c = d12q + d23q - d13q;  c = c / (2.0*d12*d23);
   e = 1.0 - b*b;
   f = 1.0 - c*c;
   if (e < 0.0 || f < 0.0)
   {
      if (e < -0.1 && f < -0.1)
      {
         fprintf(stderr,"cosomega: internal error\n");
         abort();
      }
      else
      {
         if (e < 0.0)  e = 0.0;
         if (f < 0.0)  f = 0.0;
      };
   };
   e = sqrt(e);  f = sqrt(f);
   val = (a - b*c) / (e*f);
   if (val < -1.0)  val = -1.0;
   if (val >  1.0)  val =  1.0;

   return val;
};

/* other functions */

// expanding the bounds defining the 3D boxes (to help SPG to converge)
void expandBounds(int n,int *e,double **X,double **lX,double **uX,double factor)
{
   int i,j,k;
   int min,max;

   // verifying whether there are new not-yet expanded boxes
   min = e[0];  max = min;
   for (i = 1; i < n; i++)
   {
      if (min > e[i])  min = e[i];
      if (max < e[i])  max = e[i];
   };

   if (min != max)
   {
      // reversing the previous bound-expansion effect
      for (i = 0; i < n; i++)
      {
         for (k = 0; k < 3; k++)
         {
            j = e[i];
            while (j > 0 && lX[k][i] + (j+1)*factor < X[k][i] && X[k][i] < uX[k][i] - (j+1)*factor)
            {
               lX[k][i] = lX[k][i] + j*factor;
               uX[k][i] = uX[k][i] - j*factor;
               j--;
            };
         };
         e[i] = 0.0;
      };
   }
   else
   {
      // expanding the bounds
      for (i = 0; i < n; i++)
      {
         e[i]++;
         for (k = 0; k < 3; k++)
         {
            lX[k][i] = lX[k][i] - e[i]*factor;
            uX[k][i] = uX[k][i] + e[i]*factor;
         };
      };
   };
};

// projection of x over the real interval [a,b]: returns the projected x
// (with the use of the tolerance eps)
double projection(double x,double a,double b,double eps)
{
   if (a <= x && x <= b)
   {  
      return x;
   }
   else if (a > x)
   {  
      return a - eps;
   }
   else
   {  
      return b + eps;
   };
};

// counting the number of digits forming an integer number
int numberOfDigits(int integer)
{
   int nd = 0;
   while (integer > 0)
   {
      integer = integer/10;
      nd++;
   };
   return nd;
};

// removing the extension from the input file (if such an extension is present)
// it creates a new string with the file name without extension
char* removExtension(char *filename)
{
   int i,k;
   char* output = strdup(filename);

   k = strlen(output) - 1;
   while (k >= 0 && output[k] != '.' && output[k] != '/')  k--;
   if (output[k] == '.')
   {
      output[k] = '\0';
      for (i = k + 1; i < strlen(output); i++)  output[k] = '\0';
   };

   return output;
};

// finding the minimum value among three double values
double minimum(double a,double b,double c)
{
   double min = a;
   if (min > b)  min = b;
   if (min > c)  min = c;
   return min;
};

// finding the maximum value among three double values
double maximum(double a,double b,double c)
{
   double max = a;
   if (max < b)  max = b;
   if (max < c)  max = c;
   return max;
};

// BP usage - function invoked when too few arguments are passed to MDjeep
// -----------------------------------------------------------------------

void mdjeep_usage(void)
{
   fprintf(stderr,"mdjeep: too few arguments\n");
   fprintf(stderr,"        syntax: ./mdjeep [options] instance.nmr\n");
   fprintf(stderr," Options:\n");
   fprintf(stderr,"          -v | change input data format to previous versions (argument 0.1 or 0.2, with the same effect)\n");
   fprintf(stderr,"          -e | sets the tolerance epsilon (needs an extra argument, double, default is 0.001)\n");
   fprintf(stderr,"          -r | sets the resolution parameter (needs an extra argument, double, default is 1.0)\n");
   fprintf(stderr,"        -sym | only one symmetric half of the tree is explored (argument 1 or 2)\n");
   fprintf(stderr,"          -1 | the algorithm stops at the first solution\n");
   fprintf(stderr,"          -p | prints the best found solution in a text file\n");
   fprintf(stderr,"          -P | prints all found solutions (in the same text file)\n");
   fprintf(stderr,"             |  (when using -1, options -p and -P have the same effect)\n");
   fprintf(stderr,"          -f | specifies the output format (default is \"txt\", may be changed to \"pdb\")\n");
   fprintf(stderr,"     -consec | verifies whether the consecutivity assumption is satisfied\n");
   fprintf(stderr,"  -nomonitor | does no show the current layer number during the execution to improve performace\n");
};

