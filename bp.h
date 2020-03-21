/*******************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - header file
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    solution, option and info structures updated, new prototypes added
              Jun 28 2019  v.0.3.0  new and more efficient organization of data structures and functions
              Mar 21 2020  v.0.3.1  adding triplet structure and new function prototypes
********************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>

// Data Structures
// ---------------

// Omega angle interval (to be used in the linked list omegaList)
typedef struct omega Omega;
struct omega
{
   double l;      // omega angle lower bound
   double u;      // omega angle upper bound
   Omega *prev;   // pointer to previous angle
   Omega *next;   // pointer to next angle
};

// linked list of Omega intervals
typedef struct omegalist omegaList;
struct omegalist
{
   Omega *first;  // the first omega interval
};

// Reference distance
// (with reference vertex, distance interval [lb,ub], and the pointer to the next, if any)
typedef struct Reference REFERENCE;
struct Reference
{
   int otherId;      // id of the reference vertex
   double lb;        // distance lower bound
   double ub;        // distance upper bound
   REFERENCE *next;  // pointer to the next reference distance
};

// Triplet of references
// (for every vertex, the best triplet of reference vertices can be computed once and kept in memory)
typedef struct tpt triplet;
struct tpt
{
   REFERENCE *r1;
   REFERENCE *r2;
   REFERENCE *r3;
};

// Vertex (Id,groupId,Name,Group, and linked list of reference distances)
typedef struct Vertex VERTEX;
struct Vertex
{
   int Id;          // the vertex id
   int groupId;     // the vertex group id
   char *Name;      // char string containing the vertex name
   char *Group;     // char string containing the vertex group name
   REFERENCE *ref;  // pointer to the first reference distance
};

// SEARCH: collection of data and additional memory space necessary during the search (BP and SPG)
typedef struct Search SEARCH;
struct Search
{
   bool *sym;                    // boolean vector indicating whether a tree layer is symmetric or not
   triplet *refs;                // triplet vector containing the reference vertices for every vertex
   int *e;                       // counters of times boxes are expanded (to help SPG to converge)
// double *costheta,*cosomega;   // vectors of double with the cosine of theta and omega angles
   double **lX,**uX;             // bounds defining the boxes (necessary for SPG)
   double **pX;                  // additional memory for SPG
   double *y,*dL,*dU;            // additional memory for SPG
   double *gy,*sy,*yp,*gyp;      // additional memory for SPG
   double **gX,**sX,**Xp,**gXp;  // additional memory for SPG
   double **DX,**YX,**ZX;        // additional memory for SPG
   double *Dy,*Yy,*Zy;           // additional memory for SPG
   double *memory;               // additional memory for SPG
   double be;                    // bound expansion variable (for SPG)
   double pi;                    // pi
};

// options
typedef struct option OPTION;
struct option
{
   int version;   // is it necessary to load instances in the format of versions 0.1 and 0.2?
   int print;     // 0 = no print; 1 = print the best solution; >1 = print all solutions (default 0)
   int format;    // default format is "xyz" (0); it can be changed to "pdb" (1)
   int allone;    // all vs. one solution: 1 = only one solution (default 0)
   int symmetry;  // 0 = the entire tree is explored (default)
                  // 1 = only the first symmetric half of the tree is explored
                  // 2 = only the second symmetric half of the tree is explored
   bool monitor;  // if false, the small monitor indicating the currently explored layer is not printed
   double r;      // resolution parameter (default 1.0)
   double eps;    // tolerance epsilon (default 0.001)
};

// info
typedef struct information INFORMATION;
struct information
{
   bool exact;       // true if the instance has only exact distances
   int ndigits;      // number of digits forming the largest vertex rank
   int ncalls;       // number of BP calls
   int nspg;         // number of SPG calls
   int nspgok;       // number of successful SPG calls
   int nsols;        // number of solutions found by BP
   int maxsols;      // maximum number of solutions (set to 10000, the user cannot modify it)
   int pruning;      // number of times the pruning test pruned out tree branches
   int best_sol;     // integer label of best solution
   double best_mde;  // MDE function value in the best found solution
   double best_lde;  // LDE function value in the best found solution
   char *filename;   // name of input file
   char *output;     // name of output file
};

// Function prototypes
// -------------------

// bp.c
void bp(int i,int n,int m,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info);
void bp_exact(int i,int n,int m,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info);
void intHandler(int a);  // signal catcher

// distance.c
double pairwise_distance(double xA,double yA,double zA,double xB,double yB,double zB);
double distance(int i,int j,double **X);
REFERENCE* initReference(int otherId,double lb,double ub);
REFERENCE* addDistance(REFERENCE *ref,int otherId,double lb,double ub);
int otherVertexId(REFERENCE *ref);
double lowerBound(REFERENCE* ref);
double upperBound(REFERENCE* ref);
int numberOfDistances(REFERENCE *ref);
int numberOfExactDistances(REFERENCE *ref,double eps);
double rangeOfDistance(REFERENCE *ref);
bool isExactDistance(REFERENCE *ref,double eps);
REFERENCE* nextExactDistance(REFERENCE *current,double eps);
bool isIntervalDistance(REFERENCE *ref,double eps);
REFERENCE* nextIntervalDistance(REFERENCE *current,double eps);
void printDistances(REFERENCE *ref);
REFERENCE* freeReference(REFERENCE *ref);

// vertex.c
void initVertex(VERTEX *v,int Id,int groupId,char *Name,char *Group);
int getVertexId(VERTEX v);
int getVertexGroupId(VERTEX v);
char* getVertexName(VERTEX v);
char* getVertexGroupName(VERTEX v);
REFERENCE* getReference(VERTEX* v,int i,int j);
int totalNumberOfDistances(int n,VERTEX *v);
int totalNumberOfExactDistances(int n,VERTEX *v,double eps);
double* getDistanceList(int n,VERTEX *v);
bool initialClique(int n,VERTEX *v,double eps);
int isDDGP(int n,VERTEX *v,double eps,bool clique);
bool isDMDGP(int n,VERTEX *v,double eps,bool ddgp);
void findSymmetries(int n,VERTEX *v,bool *sym);
triplet findReferences(int n,VERTEX *v,int id,double eps,bool onlyExact);
void printDistanceList(int n,VERTEX *v,bool symmetric);
void printVertex(VERTEX v);
VERTEX* freeVertex(int n,VERTEX *v);

// matrices.c   
double* allocateVector(size_t n);
void copyVector(size_t n,double *source,double *dest);
void differenceVector(size_t n,double *a,double *b,double *c);
double normVector(size_t n,double *v);
bool areSameVector(size_t,double *v1,double *v2);
void crossProdVector(double *v1,double *v2,double *res);
void printVector(size_t n,double *v);
double* freeVector(double *v);
double** allocateMatrix(size_t n,size_t m);
void copyMatrix(size_t n,size_t m,double **source,double **dest);
void copyCenterMatrix(size_t n,size_t m,double **source,double **dest);
void differenceMatrix(size_t n,size_t m,double **A,double **B,double **C);
bool areSameMatrix(size_t n,size_t m,double **A,double **B);
void UMatrix(int i3,int i2,int i1,int i,double **X,double *U);
void genCoordinates(int i1,int i,double **X,double *U,double di1i,double ctheta,double stheta,double comega,double somega);
void printMatrix(size_t n,size_t m,double **a);
double** freeMatrix(size_t n,double **a);

// objfun.c
double compute_mde(int n,VERTEX *v,double **X,double eps);
double compute_lde(int n,VERTEX *v,double **X,double eps);
double compute_stress(int n,VERTEX *v,double **X,double *y);
void stress_gradient(int n,VERTEX *v,double **X,double *y,double **gX,double *gy,double *memory);

// pruningtest.c
int DDF(int id,VERTEX *v,double **X,double eps);
int BoxDDF(int id,VERTEX *v,double **lX,double **uX,double eps);

// spg.c
double scalarProd(int n,double **X1,double **X2,int m,double *y1,double *y2);
double norm(int n,double **X,int m,double *y);
int spg(int n,VERTEX *v,double **X,SEARCH S,int *its,double *obj);

// utils.c
omegaList initOmegaList(double l,double u);
Omega* firstOmegaInterval(omegaList L);
Omega* lastOmegaInterval(omegaList L);
double omegaIntervalLowerBound(Omega *current);
double omegaIntervalUpperBound(Omega *current);
bool omegaIntervalHasNext(Omega *current);
Omega* omegaIntervalNext(Omega *current);
bool omegaIntervalHasPrev(Omega *current);
Omega* omegaIntervalPrev(Omega *current);
bool omegaIntervalHasNextAlongDirection(Omega *current,bool asNext);
Omega* omegaIntervalNextAlongDirection(Omega *current,bool asNext);
void attachNewOmegaInterval(Omega *current,double l,double u);
void splitOmegaIntervals(Omega *current,double radius,double resolution);
int numberOfOmegaIntervals(Omega *current);
void printOmegaList(omegaList L);
void freeNextOmega(Omega *a);
omegaList freeOmegaList(omegaList L);
double costheta(int i,int j,int k,VERTEX *v,double **X);
double cosomega(int i3,int i2,int i1,int i,VERTEX *v,double **X,double range);
void expandBounds(int n,int *e,double **X,double **lX,double **uX,double factor);
double projection(double x,double a,double b,double eps);
int numberOfDigits(int integer);
char* removExtension(char *filename);
double minimum(double a,double b,double c);
double maximum(double a,double b,double c);
void mdjeep_usage(void);

// print.c
void printfile(int n,VERTEX *v,double **X,char *filename,int s);
void printpdb(int n,VERTEX *v,double **X,char *filename,int s);

// splitime.c
char* splitime(struct timeval start,struct timeval end);

