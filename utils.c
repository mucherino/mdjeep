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
              May 19 2020  v.0.3.2  functions for analyzing and reading input distance file added
                                    functions createBox and reCenterBounds added
                                    function expandBounds reimplemented
              Apr 13 2020  v.0.3.2  patch
*****************************************************************************************************/

#include "bp.h"

int errno;
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

// this function indicates whether another omega interval precedes the current one
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

// this function "attaches" a new omega interval to the omega list starting with the given omega interval
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

// this function splits the omega interval list starting with "current" in subintervals
// if their "arclength" (radius*angle) is larger than the given threshold "resolution"
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

// this function counts the number of omega intervals that are next from a given interval
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
// -> range is a real between 0.0 and 1.0 indicating the value of the distance in the interval distance between i3 and i
// -> when values of range equal to 0.0 or 1.0, the first (>=0.0 or <=1.0) value is chosen for which the angle is feasible
//   (so that the cosine can be actually computed)
// -> eps is the step between two consecutive range values, to be modified until the cosine can be computed
// -> the returning value is the cosine value on success, it is equal to -2.0 otherwise
double cosomega(int i3,int i2,int i1,int i,VERTEX *v,double **X,double range,double eps)
{
   REFERENCE *r12,*r13,*r14,*r23,*r24,*r34;
   double d12,d13,d14,d23,d24,d34;
   double d12q,d13q,d14q,d23q,d24q,d34q;
   double a,b,c,e,f;
   double r,val;

   // references
   r12 = getReference(v,i3,i2);  r13 = getReference(v,i3,i1);  r23 = getReference(v,i2,i1);  // may not be available
   r14 = getReference(v,i3,i);   r24 = getReference(v,i2,i);   r34 = getReference(v,i1,i);   // known because of discretization
   if (r14 == NULL || r24 == NULL || r34 == NULL)
   {
      fprintf(stderr,"cosomega: internal error; it looks like the discretization assumptions are not satisfied\n");
      abort();
   };

   // fixed discretization distances to vertex i
   d24 = lowerBound(r24);  d24q = d24*d24;
   d34 = lowerBound(r34);  d34q = d34*d34;

   // non-discretization distances
   if (r12 != NULL)  d12 = lowerBound(r12);  else  d12 = distance(i3,i2,X);
   d12q = d12*d12;
   if (r13 != NULL)  d13 = lowerBound(r13);  else  d13 = distance(i3,i1,X);
   d13q = d13*d13;
   if (r23 != NULL)  d23 = lowerBound(r23);  else  d23 = distance(i2,i1,X);
   d23q = d23*d23;

   // computing cosine of first feasible angle
   r = range;  e = -1.0;  f = -1.0;
   while (r >= 0.0 && r <= 1.0 && (e < 0.0 || f < 0.0))
   {
      // variable discretization distance (depends on value of range)
      d14 = lowerBound(r14) + r*(upperBound(r14) - lowerBound(r14));  d14q = d14*d14;

      // computing cosine
      a = d12q + d24q - d14q;  a = a / (2.0*d12*d24);
      b = d24q + d23q - d34q;  b = b / (2.0*d24*d23);
      c = d12q + d23q - d13q;  c = c / (2.0*d12*d23);
      e = 1.0 - b*b;
      f = 1.0 - c*c;

      // preparing for next loop?
      if (range == 0.0)  r = r + eps;
      if (range == 1.0)  r = r - eps;
      if (r == range)  break;
   };

   // perform correction if error is "relatively small"
   if (e < 0.0 && e > -1.0)  e = 0.0;
   if (f < 0.0 && f > -1.0)  f = 0.0;

   // stop if feasible angle was not found
   if (e < 0.0 || f < 0.0)  return -2.0;

   // final computation of cosine
   e = sqrt(e);  f = sqrt(f);
   val = (a - b*c) / (e*f);
   if (val < -1.0)  val = -1.0;
   if (val >  1.0)  val =  1.0;
   return val;
};

/* functions for "string" management */

// counting the number of digits forming a strictly positive integer number
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

// counting the number of digits forming the decimal part of a double number
int precisionOf(double real)
{
   int nd = 0;
   if (real < 0.0)  real = -real;
   while (real != floor(real))
   {
      real = 10.0*real;
      nd++;
   };
   return nd;
};

// function to verify whether the input char is the delimiter for the end of a string
bool isLastChar(char c)
{
   return c == '\0';
};

// function to verify whether the input char is a "separator"
// -> it is compared to the separator specified in sep, but also to blank chars and tabs
bool isSeparator(char c,char sep)
{
   if (c == sep)   return true;
   if (c == ' ')   return true;
   if (c == '\t')  return true;
   return false;
};

// function to move the char pointer to the next char that is not a blank char nor a tab
// -> the returning value is the pointer to the next char that is non-blank and not a tab
//   (NULL if it doesn't exit)
char* nextNonBlank(char *c)
{
   while (!isLastChar(c[0]) && (c[0] == ' ' || c[0] == '\t'))  c++;
   if (isLastChar(c[0]))  return NULL;
   return c;
};

// function to verify whether the input char is delimiter for "new line"
bool isNewLineDelimiter(char c)
{
   if (c == '\n')  return true;
   if (c == '\r')  return true;
   return false;
};

// function to remove the "ending" chars from an array of chars
// -> ending chars are blank chars, \n and \r
size_t removEndingChars(char *c)
{
   size_t l = strlen(c) - 1;
   while (l >= 0 && (c[l] == ' ' || isLastChar(c[l]) || isNewLineDelimiter(c[l])))
   {
      c[l] = '\0';
      l--;
   };
   return l;
};

// function which looks for the next colon (:) in a char string passing over blank chars (or tabs)
// -> the returning value is the pointer to the colon, if it is found before '\0' and any other char different from blank or tab
char* nextColon(char *c)
{
   c = nextNonBlank(c);
   if (c[0] == ':')  return c;
   return NULL;
};

// function to verify whether a string represents an integer number
// -> the string is supposed to be allocated and contain the ending '\0'
// -> the string cannot contain anything else other than the integer number
bool isInteger(char *c)
{
   char *pointer;

   // no digits
   if (isLastChar(c[0]))  return false;

   // first digit (sign?)
   if (c[0] == '-' || c[0] == '+')  c++;

   // first (actual) digit
   if (!isLastChar(c[0]))
   {
      if (isdigit(c[0]))
      {
         if (c[0] == '0' && !isLastChar(c[1]))  return false;  // avoiding integers starting with 0's
      }                                                       // (these are valid integers for strtol)
      else return false;
   };

   // all digits
   errno = 0;  // global error variable
   strtol(c,&pointer,10);
   if (errno != 0)  return false;  // it's not an integer
   if (pointer < c + strlen(c))  return false;  // the string does not only contain the integer

   return true;
};

// function to verify whether a string represents a real number
// -> the string is supposed to be allocated and contain the ending '\0'
// -> the string cannot contain anything else other than the real number
bool isReal(char *c)
{
   char *pointer;

   // no digits
   if (isLastChar(c[0]))  return false;

   // first digit (sign?)
   if (c[0] == '+' || c[0] == '-')  c++;

   // first (actual) digit
   if (!isLastChar(c[0]))
   {
      if (isdigit(c[0]))
      {
         if (c[0] == '0' && c[1] != '.')  return false;   // avoiding reals starting with 0's before actual value
      }
      else return false;
   };

   // all digits
   errno = 0;  // global error variable
   strtod(c,&pointer);
   if (errno != 0)  return false;  // it's not a real
   if (pointer < c + strlen(c))  return false;  // the string does not only contain the real number

   return true;
};

// removing the extension from the input file (if there is such an extension)
// -> it creates a new string with the file name without extension
//   (memory is allocated for the new string)
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

// detecting the types in a character string of words
// -> the array of char c needs to be a valid pointer (with allocated memory)
// -> "sep" is the separator (one character)
// -> the binary output is the format specified with this internal code: 01=integer,10=real,11=alpha,00=none/end)
//    0UL is returned when the char array is empty,
//                 or when the char array contains more than 4*sizeof(unsigned long) words
unsigned long detectTypes(char *c,char sep)
{
   int i,i0,k;
   int n,nwords;
   char tmp,*pointer;
   unsigned long ct;
   unsigned long type = 0UL;

   if (!isLastChar(c[0]))
   {
      // how many values on the line?
      i = 0;  nwords = 0;
      while (!isLastChar(c[i]) && isSeparator(c[i],sep))  i++;
      i0 = i;
      if (!isLastChar(c[i]))  nwords = 1;
      while (!isLastChar(c[i]))
      {
         if (isSeparator(c[i],sep))
         {
            i++;
            while (!isLastChar(c[i]) && isSeparator(c[i],sep))  i++;
            if (!isLastChar(c[i]))  nwords++;
         }
         else i++;
      };
      n = i + 1;

      // nothing to do if the line contains no words      
      if (nwords == 0)  return 0UL;

      // verifying the word types (with isInteger and isReal)
      i = i0;  k = 0;
      pointer = c + i;
      while (i < n)
      {
         if (isLastChar(c[i]) || isSeparator(c[i],sep))
         {
            if (!isLastChar(pointer[0]))
            {
               tmp = c[i];  c[i] = '\0';
               if (isInteger(pointer))
                  ct = 1;
               else if (isReal(pointer))
                  ct = 2;
               else
                  ct = 3;
               type = (type << 2) | ct;
               k = k + 2;  if (k > 8*sizeof(unsigned long))  return 0UL;
               c[i] = tmp;  i++;
               while (!isLastChar(c[i]) && isSeparator(c[i],sep))  i++;
               pointer = c + i;
            }
            else break;
         }
         else i++;
      };
   };

   return type;
};

/* other functions */

// creating the box with a predefined range
// -> i is the vertex rank for which the box needs to be created
// -> X[k][i] contains the coordinates of the vertex (k = 1,2,3)
// -> range is the specified interval range over the 3 dimensions
// -> the output are lX[k][i] and uX[k][i], for all i and k
void createBox(int i,double **X,double range,double **lX,double **uX)
{
   int k;

   range = 0.5*range;
   for (k = 0; k < 3; k++)
   {
      lX[k][i] = X[k][i] - range;
      uX[k][i] = X[k][i] + range;
   };
};

// expanding the bounds of a given vertex box as long as the added parts contain feasible positions
// -> i is the rank of the vertex whose box needs to be expanded
// -> v is the VERTEX array
// -> [lX,uX] is the list of boxes
// -> be is the "bound expanding" factor (ie, the added range at every expansion)
// -> eps is the tolerance error used for verify whether the added part to the box contains feasible positions
// -> in output, the list [lX,uX] is unchanged expect for the ith box, which is expanded 
//   (at least one expansion step is performed)
void expandBounds(int i,VERTEX *v,double **lX,double **uX,double be,double eps)
{
   int k,m;
   double *tmp,*d1,*d2,*d3,*d4;
   bool min,max,feasibility;
   REFERENCE *ref;

   // computing the number of reference distances
   m = numberOfDistances(v[i].ref);

   // memory allocation to keep the distances (avoids to compute them twice)
   d1 = allocateVector(m);  d2 = allocateVector(m);
   d3 = allocateVector(m);  d4 = allocateVector(m);

   // computing the box distance to reference vertices
   k = 0;  ref = v[i].ref;
   while (ref != NULL)
   {
      d1[k] = box_distance(i,otherVertexId(ref),lX,uX,d2+k);
      ref = ref->next;
      k++;
   };

   do // performing the expansions
   {
      // one (more) expansion
      lX[0][i] = lX[0][i] - be;  uX[0][i] = uX[0][i] + be;
      lX[1][i] = lX[1][i] - be;  uX[1][i] = uX[1][i] + be;
      lX[2][i] = lX[2][i] - be;  uX[2][i] = uX[2][i] + be;

      // recomputing the distance with the expanded box
      // and verifying feasibility of expanded part
      k = 0;  ref = v[i].ref;  feasibility = false;
      while (ref != NULL)
      {
         d3[k] = box_distance(i,otherVertexId(ref),lX,uX,d4+k);
         min = lowerBound(ref) >= d3[k] - eps && upperBound(ref) <= d1[k] + eps;
         max = lowerBound(ref) >= d2[k] - eps && upperBound(ref) <= d4[k] + eps;
         feasibility = feasibility || min || max;
         ref = ref->next;
         k++;
      };

      // swapping distances [d1,d2] and [d3,d4] to prepare next iteration
      tmp = d1;  d1 = d3;  d3 = tmp;
      tmp = d2;  d2 = d4;  d4 = tmp;
   }
   while (feasibility);

   // ending
   freeVector(d1);  freeVector(d2);
   freeVector(d3);  freeVector(d4);
};

// recentering the bounds defining the vertex boxes around new coodinates in X
// -> from a mathematical point of view, the recentering procedure can be described as follows:
//    1. every box is translated in the 3D space so that its center corresponds now to the position in X
//    2. the intersection between the old box and the new translated box is performed
//    3. the intersection is expanded as long as the added parts contain feasible positions
// -> v is a VERTEX array, with size n
// -> X is the corresponding set of 3D coordinates
// -> [lX,uX] is the list of vertex boxes, to be recentered
// -> be is the "bound expanding" factor
// -> eps is the tolerance error used for verify whether the added part to the box contains feasible positions
// -> in output, the list [lX,uX] contains the recentered boxes
void reCenterBounds(int n,VERTEX *v,double **X,double **lX,double **uX,double be,double eps)
{
   int i;
   double range;

   // recentering the box around every position in X
   for (i = 0; i < n; i++)
   {
      // x coordinates
      range = 0.5*(uX[0][i] - lX[0][i]);
      if (X[0][i] - range > lX[0][i])  lX[0][i] = X[0][i] - range;
      if (X[0][i] + range < uX[0][i])  uX[0][i] = X[0][i] + range;

      // y coordinates
      range = 0.5*(uX[1][i] - lX[1][i]);
      if (X[1][i] - range > lX[1][i])  lX[1][i] = X[1][i] - range;
      if (X[1][i] + range < uX[1][i])  uX[2][i] = X[2][i] + range;

      // z coordinates
      range = 0.5*(uX[2][i] - lX[2][i]);
      if (X[2][i] - range > lX[2][i])  lX[2][i] = X[2][i] - range;
      if (X[2][i] + range < uX[2][i])  uX[2][i] = X[2][i] + range;

      // re-expanding the bounds as long as feasible positions are added
      expandBounds(i,v,lX,uX,be,eps);
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

// finding the minimum of three double values
double minimum(double a,double b,double c)
{
   double min = a;
   if (min > b)  min = b;
   if (min > c)  min = c;
   return min;
};

// finding the maximum of three double values
double maximum(double a,double b,double c)
{
   double max = a;
   if (max < b)  max = b;
   if (max < c)  max = c;
   return max;
};

/* Help */

// BP usage - function invoked when too few arguments are passed to MDjeep
// -----------------------------------------------------------------------

void mdjeep_usage(void)
{
   fprintf(stderr,"mdjeep: too few arguments\n");
   fprintf(stderr,"        syntax: ./mdjeep [options] MDfile.mdf\n");
   fprintf(stderr," Options:\n");
   fprintf(stderr,"          -1 | the specified method stops at the first solution (always true for SPG)\n");
   fprintf(stderr,"          -l | specifies after how many solutions the method should stop (applies only to BP)\n");
   fprintf(stderr,"        -sym | only one symmetric half of the tree is explored (for BP, argument may be 1 or 2)\n");
   fprintf(stderr,"          -p | prints the best found solution in a text file\n");
   fprintf(stderr,"          -P | prints all found solutions (in the same text file)\n");
   fprintf(stderr,"             |  (when using -1, options -p and -P have the same effect)\n");
   fprintf(stderr,"          -f | specifies the output format (default is \"xyz\", may be changed to \"pdb\")\n");
   fprintf(stderr,"     -consec | verifies whether the consecutivity assumption is satisfied\n");
   fprintf(stderr,"  -nomonitor | does not show the current layer number during the execution to improve performance\n");
   fprintf(stderr,"          -r | obsolete, resolution parameter can now be specified in MDfile (method field)\n");
   fprintf(stderr,"          -e | obsolete, tolerance epsilon can now be specified in MDfile (method field)\n");
   fprintf(stderr,"          -v | obsolete, file formats can now be specified in MDfile (instance field)\n");
   fprintf(stderr," Please refer to the documentation for the MDfile syntax.\n");
};

