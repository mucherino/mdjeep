/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - header file
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// data structures
// ---------------

// probl - contains the problem - list of distances

struct probl
{
   int i;           // label first atom
   int j;           // label second atom
   char iatom[3];   // name of first atom (i)
   char jatom[3];   // name of second atom (j)
   char iamino[4];  // amino acid of first atom (i)
   char jamino[4];  // amino acid of second atom (j)
   double l;        // lower bound on the distance (i,j)
   double u;        // upper bound on the distance (i,j)
};

typedef struct probl PROBL;

// sol - contains a solution - list of coordinates and torsion angles

struct solution
{
   char atom[3];    // atom name
   char amino[4];   // amino acid name
   double costheta; // cosine of the angle among this atom and the two precedessors
   double cosomega; // cosine of the torsion angle among this atom and the three precedessors
   short branch;    // left vs right branch on the binary tree
   double x;        // x coordinate
   double y;        // y coordinate
   double z;        // z coordinate
};

typedef struct solution SOLUTION;

// option - list of options

struct option
{
   int print;         // 0 = no print; 1 = print the best solution; >1 = print all solutions (default 0)
   int allone;        // all vs. one solution: 1 = only one solution is required (default 0)
   int symmetry;      // 1 = symmetric solutions are not computed (default 0)
   int triangular;    // 1 = triangular inequalities are checked (default 0)
   double eps;        // tolerance epsilon (default 0.001)
};

typedef struct option OPTION;

// info - information on the execution of the BP algorithm

struct information
{
   int nsols;             // number of solutions found by BP
   int pruning;           // number of times the pruning test has been applied
   int best_sol;          // label of best solution
   double best_lde;       // LDE function value in the best found solution
   char filename[200];    // name of input file
};

typedef struct information INFORMATION;

// partial LDE - contains partial values for the LDE function

struct partialde
{
   int n;
   double val;
};

typedef struct partialde PARTIALDE;

// functions

int bp(int i,int n,int m,PROBL *p,int **ind,SOLUTION *sol,OPTION op,PARTIALDE *plde,double **Q,INFORMATION *info);  // bp.c

int pruningtest(int natoms,PROBL *p,int **ind,SOLUTION *sol,OPTION op,INFORMATION *info,PARTIALDE *plde);  // pruningtest.c

int printpdb(int n,SOLUTION *sol,double lde,char *filename,int s);  // printfile.c

void set_matrix(double b11,double b12,double b13,double b14,
                double b21,double b22,double b23,double b24,
                double b31,double b32,double b33,double b34,
                double b41,double b42,double b43,double b44,
                double **Q,int i);  // matrices.c

void matrix_prod(double *Q1,double *Q2,double *Qout);  // matrices.c

double costheta(int i,int j,int k,PROBL *p,int **ind);  // utils.c

double cosomega(int i,int j,int k,int h,PROBL *p,int **ind);  // utils.c

double distance(double x1,double y1,double z1,double x2,double y2,double z2);  // utils.c

void bp_usage(void);  // utils.c



