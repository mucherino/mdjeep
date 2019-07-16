/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - main program
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
*****************************************************************************/

#include "bp.h"

double INFTY = 1.e+30;

int main(int argc, char *argv[])
{
   int i,j,k;
   int ij,jk,ki;
   int n,m;
   int fidx;
   int **ind;
   double **Q;
   double time;
   double val1,val2;
   char ch1[4],ch2[4];
   char ch3[4],ch4[4];
   long t1,t2;
   PROBL *p;
   SOLUTION *sol;
   OPTION op;
   INFORMATION info;
   PARTIALDE *plde;
   FILE *input;

   // starting
   fprintf(stderr,"MD-jeep v0.1\n");
   fprintf(stderr,"GNU General Public License\n");
   fprintf(stderr,"Copyright (C) Mucherino, Liberti, Lavor, Maculan\n");

   // checking input arguments
   if (argc < 2)
   {
      bp_usage();
      return 1;
   };

   // default values for structures 'op' and 'info'
   op.print = 0;  op.allone = 0;  op.symmetry = 0;  op.triangular = 0;  op.eps = 0.001;
   info.nsols = 0;  info.pruning = 0;  info.best_sol = 0;  info.best_lde = INFTY;

   // checking input options
   fidx = 1;
   while (fidx < argc - 1) 
   { 
      if (strncmp(argv[fidx], "-e", 2) == 0)
      {
         if (fidx + 1 >= argc - 1) 
         {
            fprintf(stderr,"mdjeep: error: -e flag requires a floating point argument\n");
            return 1;
         };
         op.eps = atof(argv[fidx+1]);
         if (op.eps < 0 || op.eps > 100)
         {
            fprintf(stderr,"mdjeep: error: %g is not an acceptable tolerance\n",op.eps);
            return 1;
         };
         fidx = fidx + 2;
      }
      else if (strncmp(argv[fidx], "-p", 2) == 0)
      {
         op.print = 1;
         fidx++;
      }
      else if (strncmp(argv[fidx], "-P", 2) == 0)
      {
         op.print = 2;
         fidx++;
      }
      else if (strncmp(argv[fidx], "-1", 2) == 0)
      {
         op.allone = 1;
         fidx++;
      }
      else if (strncmp(argv[fidx], "-s", 2) == 0)
      {
         op.symmetry = 1;
         fidx++;
      }
      else if (strncmp(argv[fidx], "-t", 2) == 0)
      {
         op.triangular = 1;
         fidx++;
      }
      else
      {
         fprintf(stderr,"mdjeep: warning: '%s' is not an option for mdjeep\n",argv[fidx]);
         fidx++;
      };
   };

   // checking if there is an input file
   if (fidx >= argc) 
   {
      fprintf(stderr,"mdjeep: error: no input filename found on command line\n");
      return 1;
   };

   // welcome message
   fprintf(stderr,"mdjeep: tolerance epsilon = %g\n",op.eps);
   if (op.print == 1)  fprintf(stderr,"mdjeep: the best solution ");
   if (op.print == 2)  fprintf(stderr,"mdjeep: all solutions ");
   if (op.print != 0)  fprintf(stderr,"will be printed in PDB format\n");
   if (op.allone == 1)  fprintf(stderr,"mdjeep: only one solution is required by the user\n");
   if (op.symmetry == 1)  fprintf(stderr,"mdjeep: symmetric solutions are not computed\n");

   // opening input file
   strcpy(info.filename,argv[fidx]);
   input = fopen(argv[fidx],"r");
   if (input == NULL)
   {
      fprintf(stderr,"mdjeep: error while reading input file '%s'\n",argv[fidx]);
      return 1;
   };

   // computing the number of distances in the input file
   while (EOF != fscanf(input,"%d %d %lf %lf %s %s %s %s\n",&i,&j,&val1,&val2,ch1,ch2,ch3,ch4))  m++;
   rewind(input);

   // allocation of memory for the problem
   p = (PROBL*)calloc(m,sizeof(PROBL));

   // reading distances and computing the number of atoms
   n = 0;
   for (k = 0; k < m; k++)
   {
      fscanf(input,"%d %d %lf %lf %s %s %s %s\n",&p[k].i,&p[k].j,&p[k].l,&p[k].u,p[k].iatom,p[k].jatom,p[k].iamino,p[k].jamino);
      if (n < p[k].i)  n = p[k].i;
      if (n < p[k].j)  n = p[k].j;
   };

   // closing input file
   fclose(input);

   // checking lower and upper bounds
   for (k = 0; k < m; k++)
   {
      if (p[k].l != p[k].u)
      {
         fprintf(stderr,"mdjeep: error: lower and upper bounds must be the same in this version\n");
         fprintf(stderr,"               this feature will be implemented in future versions\n");
         return 1;
      };
   };

   // printing problem details
   fprintf(stderr,"mdjeep: file '%s' read: %d atoms / %d distances\n",argv[fidx],n,m);

   // defining the matrix of pointers 'ind'
   // (i,j) ---> position of distance (i,j) in p
   ind = (int**)calloc(n+1,sizeof(int*));
   for (i = 0; i < n+1; i++)  ind[i] = (int*)calloc(n+1,sizeof(int));
   for (i = 0; i < n+1; i++)  for (j = 0; j < n+1; j++)  ind[i][j] = -1;  // -1 indicates that the distance is unavailable
   for (k = 0; k < m; k++)  ind[p[k].i][p[k].j] = k;
   for (k = 0; k < m; k++)  ind[p[k].j][p[k].i] = k;

   // checking if the problem is discretizable
   for (i = 1; i < n-3; i++)
   {
      if (ind[i][i+1] == -1 || ind[i][i+2] == -1 || ind[i][i+3] == -1)
      {
         fprintf(stderr,"mdjeep: error: the problem in input is not discretizable\n");
         fprintf(stderr,"    distance(s) ");
         for (j = i+1; j <= i+3; j++)  if (ind[i][j] == -1)  fprintf(stderr,"(%d,%d) ",i,j);
         fprintf(stderr,"not available\n");
         fprintf(stderr,"    stopping here ... other necessary distances may also not be available\n");
         return 1;
      };
   };

   // checking the triangular inequalities (if required)
   if (op.triangular == 1)
   {
      for (i = 1; i <= n-2; i++)
      {
         for (j = i+1; j <= n-1; j++)
         {
            for (k = j+1; k <= n; k++)
            {
               ij = ind[i][j];  jk = ind[j][k];  ki = ind[k][i];
               if (ij != -1 && jk != -1 && ki != -1)
               {
                  if (p[ij].l >= p[jk].l + p[ki].l || p[jk].l >= p[ij].l + p[ki].l || p[ki].l >= p[ij].l + p[jk].l)
                  {
                     fprintf(stderr,"mdjeep: error: the input distances do not satisfy the triangular inequalities\n");
                     fprintf(stderr,"               triplet (%d,%d,%d): distances %g, %g, %g\n",i,j,k,p[ij].l,p[jk].l,p[ki].l);
                     fprintf(stderr,"        stopping here ... other triangular inequalities may be not satisfied\n");
                     return 1;
                  };
               };
            };
         };
      };
      fprintf(stderr,"mdjeep: all triangular inequelities are satisfied\n");
   };

   // allocation of memory for the solutions
   sol = (SOLUTION*)calloc(n,sizeof(SOLUTION));

   // defining atom and amino in sol
   for (k = 0; k < m; k++)
   {
      i = p[k].i - 1;
      strcpy(sol[i].atom,p[k].iatom);
      strcpy(sol[i].amino,p[k].iamino);
      j = p[k].j - 1;
      strcpy(sol[j].atom,p[k].jatom);
      strcpy(sol[j].amino,p[k].jamino);
   };

   // computing the cosine of theta
   for (i = 2; i < n; i++)
   {
      sol[i].costheta = costheta(i-1,i,i+1,p,ind);
   };

   // computing the cosine of the torsion angles
   for (i = 3; i < n; i++)
   {
      sol[i].cosomega = cosomega(i-2,i-1,i,i+1,p,ind);
      if (sol[i].cosomega < -1.0 || sol[i].cosomega > 1.0)
      {
         fprintf(stderr,"mdjeep: error while computing the torsion angle for the quadruplet (%d,%d,%d,%d)\n",i-2,i-1,i,i+1);
         return 1;
      };
   };

   // allocating memory for the matrices Q 
   // and the vector plde (partial LDE function value)
   plde = (PARTIALDE*)calloc(n,sizeof(PARTIALDE));
   Q = (double**)calloc(n,sizeof(double*));
   for (i = 0; i < n; i++)  Q[i] = (double*)calloc(12,sizeof(double));  // 4x4 matrices

   // calling BP
   t1 = clock();
   fprintf(stderr,"mdjeep: exploring the binary tree ... ");
   bp(1,n,m,p,ind,sol,op,plde,Q,&info);
   fprintf(stderr,"done!\n");
 
   // monitoring time
   t2 = clock();
   time = (double)(t2 - t1)/CLOCKS_PER_SEC; 

   // printing the results
   fprintf(stderr,"mdjeep: %d solutions found\n",info.nsols);
   fprintf(stderr,"mdjeep: %d branches were pruned\n",info.pruning);
   if (info.nsols > 0)  fprintf(stderr,"mdjeep: best LDE function (solution #%d): %g\n",info.best_sol,info.best_lde);
   fprintf(stderr,"mdjeep: CPU time (in seconds) = %g\n",time);

   // freeing memory
   free(Q);
   free(plde);
   free(sol);
   free(ind);
   free(p); 

   return 0;
};


