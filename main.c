/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - main program
  Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
              May 10 2014  v.0.2  method based on change of basis added for coordinate computation
                                  more robust check on input instances added (atom and amino labels)
****************************************************************************************************/

#include "bp.h"

double INFTY = 1.e+30;

int main(int argc, char *argv[])
{
   int i,j,k;
   int ij,jk,ki;
   int n,m;
   int fidx;
   int fly,enough;
   int consecutivity;
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
   fprintf(stderr,"MD-jeep v0.2\n");
   fprintf(stderr,"GNU General Public License\n");
   fprintf(stderr,"Copyright (C) Mucherino, Liberti, Goncalves, Lavor, Maculan\n");

   // checking input arguments
   if (argc < 2)
   {
      bp_usage();
      return 1;
   };

   // default values for structures 'op' and 'info'
   op.method = 1;  op.print = 0;  op.allone = 0;  op.symmetry = 0;  op.triangular = 0;  op.eps = 0.001;
   info.nsols = 0;  info.maxsols = 1000;  info.pruning = 0;  info.best_sol = 0;  info.best_lde = INFTY;

   // checking input options
   fidx = 1;
   while (fidx < argc - 1) 
   {
      if (!strncmp(argv[fidx], "-m1",3))
      {
         op.method = 1;
         fidx++;
      }
      else if (!strncmp(argv[fidx], "-m2",3))
      {
         op.method = 2;
         fidx++;
      }
      else if (!strncmp(argv[fidx], "-e", 2))
      {
         if (fidx + 1 >= argc - 1) 
         {
            fprintf(stderr,"mdjeep: error: -e flag requires a floating point argument\n");
            return 1;
         };
         op.eps = atof(argv[fidx+1]);
         if (op.eps < 0 || op.eps > 10)
         {
            fprintf(stderr,"mdjeep: error: %g is not an acceptable tolerance\n",op.eps);
            return 1;
         };
         fidx = fidx + 2;
      }
      else if (!strncmp(argv[fidx], "-p", 2))
      {
         op.print = 1;
         fidx++;
      }
      else if (!strncmp(argv[fidx], "-P", 2))
      {
         op.print = 2;
         fidx++;
      }
      else if (!strncmp(argv[fidx], "-1", 2))
      {
         op.allone = 1;
         fidx++;
      }
      else if (!strncmp(argv[fidx], "-s", 2))
      {
         op.symmetry = 1;
         fidx++;
      }
      else if (!strncmp(argv[fidx], "-t", 2))
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
   fprintf(stderr,"mdjeep: method used for coordinate computation: ");
   if (op.method == 1)  fprintf(stderr,"cumulative matrix based\n");
   if (op.method == 2)  fprintf(stderr,"change of basis based\n");
   if (op.symmetry == 1)  fprintf(stderr,"mdjeep: only one symmetric half of the tree is explored\n");

   // opening input file
   strcpy(info.filename,argv[fidx]);
   input = fopen(argv[fidx],"r");
   if (input == NULL)
   {
      fprintf(stderr,"mdjeep: error while reading input file '%s'\n",argv[fidx]);
      return 1;
   };

   // computing the number of distances in the input file
   m = 0;
   while (EOF != fscanf(input,"%d %d %lf %lf %s %s %s %s\n",&i,&j,&val1,&val2,ch1,ch2,ch3,ch4))  m++;
   rewind(input);

   // checking the number of available distances
   if (m <= 0)
   {
      fprintf(stderr,"mdjeep: error: the input file seems to be empty\n");
      return 1;
   };

   // allocation of memory for the problem
   p = (PROBL*)calloc(m,sizeof(PROBL));

   // reading distances and computing the number of atoms
   n = 0;
   for (k = 0; k < m; k++)
   {
      fscanf(input,"%d %d %lf %lf %s %s %s %s\n",&i,&j,&val1,&val2,p[k].iatom,p[k].jatom,p[k].iamino,p[k].jamino);
      p[k].i = i;  p[k].j = j;
      p[k].l = val1;  p[k].u = val2;
      if (n < p[k].i)  n = p[k].i;
      if (n < p[k].j)  n = p[k].j;
   };
   fclose(input);

   // checking the number of atoms in the instance
   if (n <= 0)
   {
      fprintf(stderr,"mdjeep: error: the input file seems to be empty\n");
      return 1;
   };

   // checking lower and upper bounds
   for (k = 0; k < m; k++)
   {
      if (p[k].l != p[k].u)
      {
         fprintf(stderr,"mdjeep: error: lower and upper bounds must be equal in this version\n");
         fprintf(stderr,"               interval distances will be considered in future versions\n");
         return 1;
      };
   };

   // printing problem details
   fprintf(stderr,"mdjeep: file '%s' read: %d atoms / %d distances\n",argv[fidx],n,m);

   // defining the matrix of pointers 'ind'
   // (i,j) ---> position of distance (i,j) in p
   ind = (int**)calloc(n+1,sizeof(int*));
   for (i = 0; i < n+1; i++)  ind[i] = (int*)calloc(n+1,sizeof(int));
   for (i = 0; i < n+1; i++)  for (j = 0; j < n+1; j++)  ind[i][j] = -1;  // -1 indicates that the distance is not available
   for (k = 0; k < m; k++)  ind[p[k].i][p[k].j] = k;
   for (k = 0; k < m; k++)  ind[p[k].j][p[k].i] = k;  // ind is a symmetric matrix

   // checking the triangular inequalities (if required)
   if (op.triangular == 1)
   {
      fprintf(stderr,"mdjeep: triangular inequalities check started ... ");
      if (n > 300)  fprintf(stderr,"it may take a while ...");
      fprintf(stderr,"\n");

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

   // checking whether the input instance is discretizable
   consecutivity = 1;
   for (i = 2; i <= n; i++)
   {
      if (i == 2)
      {
         if (ind[i-1][i] == -1)  consecutivity = 0;
      }
      else if (i == 3)
      {
         if (ind[i-2][i] == -1 || ind[i-1][i] == -1)  consecutivity = 0;
      }
      else if (consecutivity == 1)  
      { 
         if (ind[i-3][i] == -1 || ind[i-2][i] == -1 || ind[i-3][i] == -1)  consecutivity = 0;
      };
      k = 0;  for (j = 1; j < i; j++)  if (ind[j][i] != -1)  k++;

      enough = 3;  if (i <= 3)  enough = i - 1;
      if (k < enough)
      {
         fprintf(stderr,"mdjeep: error: the input instance is not discretizable\n");
         fprintf(stderr,"               %1d references for atom %d\n",k,i);
         fprintf(stderr,"        stopping here ... other necessary distances may be not available\n");
         return 1;
      };
   };

   // the input instance is discretizable: message
   fprintf(stderr,"mdjeep: the input instance is discretizable; the consecutivity assumption is ");
   if (!consecutivity)  fprintf(stderr,"not ");
   fprintf(stderr,"satisfied\n");

   // verifying consecutivity vs selected method for coordinate computation
   if (op.method == 1 && !consecutivity)
   {
      op.method = 2;
      fprintf(stderr,"mdjeep: forcing 'change of basis' based method (option -m2)\n");
   };

   // allocation of memory for the solutions
   sol = (SOLUTION*)calloc(n,sizeof(SOLUTION));

   // defining atom and amino in sol
   for (i = 0; i < n; i++)  strcpy(sol[i].atom,"X");
   for (i = 0; i < n; i++)  strcpy(sol[i].amino,"X");
   for (k = 0; k < m; k++)
   {
      i = p[k].i - 1;
      if (!strncmp(sol[i].atom,"X",1) && !strncmp(sol[i].amino,"X",1))
      {
         strcpy(sol[i].atom,p[k].iatom);
         strcpy(sol[i].amino,p[k].iamino);
      }
      else
      {
         if (strncmp(sol[i].atom,p[k].iatom,2) || strncmp(sol[i].amino,p[k].iamino,3))
         {
            fprintf(stderr,"mdjeep: input file: different labels (%s,%s)<>(%s,%s) for the same atom ... ",
                    sol[i].atom,sol[i].amino,p[k].iatom,p[k].iamino);
            fprintf(stderr,"using (%s,%s)\n",sol[i].atom,sol[i].amino);
         };
      };
    
      j = p[k].j - 1;
      if (!strncmp(sol[j].atom,"X",1) && !strncmp(sol[j].amino,"X",1))
      {
         strcpy(sol[j].atom,p[k].jatom);
         strcpy(sol[j].amino,p[k].jamino);
      }
      else
      {
         if (strncmp(sol[j].atom,p[k].jatom,2) || strncmp(sol[j].amino,p[k].jamino,3))
         {
            fprintf(stderr,"mdjeep: input file: different labels (%s,%s)<>(%s,%s) for the same atom ... ",
                    sol[j].atom,sol[j].amino,p[k].jatom,p[k].jamino);
            fprintf(stderr,"using (%s,%s)\n",sol[j].atom,sol[j].amino);
         };
      };
   };

   // defining the list of reference atoms
   i = 2;  sol[i-1].ref1 = 0;
   i = 3;  sol[i-1].ref1 = 1;  sol[i-1].ref2 = 0;
   for (i = 4; i <= n; i++)
   {
      j = i - 1;  while (ind[j][i] == -1)  j--;
      sol[i-1].ref1 = j - 1;
      j--;  while (ind[j][i] == -1)  j--;
      sol[i-1].ref2 = j - 1;
      j--;  while (ind[j][i] == -1)  j--;
      sol[i-1].ref3 = j - 1;

      // j should not reach 0 because of the discr. assumptions (verified above)
      if (j <= 0)
      {
         fprintf(stderr,"mdjeep: data inconsistency: stopping\n");
         return 1;
      };
   };

   // checking whether all cosines can be pre-computed
   fly = 0;

   // computing the cosine of theta
   // (when -2, the references don't form a clique: it will be computed later)
   for (i = 2; i < n; i++)
   {
      sol[i].costheta = costheta(1+sol[i].ref2,1+sol[i].ref1,i+1,p,ind);
      if (sol[i].costheta == -2.0)  fly = 1;
   };

   // computing the cosine of the torsion angles
   // (when -2, the references don't form a clique: it will be computed later)
   for (i = 3; i < n; i++)
   {
      sol[i].cosomega = cosomega(1+sol[i].ref3,1+sol[i].ref2,1+sol[i].ref1,i+1,p,ind);
      if (sol[i].cosomega == -2.0)  fly = 1;
      if (sol[i].cosomega != -2.0)  if (sol[i].cosomega < -1.0 || sol[i].cosomega > 1.0)
      {
         fprintf(stderr,"mdjeep: error while computing the torsion angle for the quadruplet (%d,%d,%d,%d)\n",
                 1+sol[i].ref3,1+sol[i].ref2,1+sol[i].ref1,i+1);
         return 1;
      };
   };

   // message about the cosines
   if (fly == 1)  fprintf(stderr,"mdjeep: not all cosines can be pre-computed; others will be computed on-the-fly\n");

   // allocating memory for vector plde (partial LDE function value)
   plde = (PARTIALDE*)calloc(n,sizeof(PARTIALDE));

   // if option -m1, allocating memory for cumulative matrices Q 
   if (op.method == 1)
   {
      Q = (double**)calloc(n,sizeof(double*));
      for (i = 0; i < n; i++)  Q[i] = (double*)calloc(12,sizeof(double));  // 4x4 matrices
   };

   // calling BP
   t1 = clock();
   fprintf(stderr,"mdjeep: exploring the binary tree ... ");
   bp(1,n,m,p,ind,sol,op,plde,Q,&info);
   fprintf(stderr,"done!\n");
 
   // monitoring time
   t2 = clock();
   time = (double)(t2 - t1)/CLOCKS_PER_SEC; 

   // printing the results
   fprintf(stderr,"mdjeep: %d solutions found ",info.nsols);
   if (info.nsols == info.maxsols)  fprintf(stderr,"(max %d)",info.maxsols);
   fprintf(stderr,"\n");
   fprintf(stderr,"mdjeep: %d branches were pruned\n",info.pruning);
   if (info.nsols > 0)  fprintf(stderr,"mdjeep: best LDE function (solution #%d): %g\n",info.best_sol,info.best_lde);
   fprintf(stderr,"mdjeep: CPU time (in seconds) = %g\n",time);

   // freeing memory
   if (op.method == 1)  free(Q);
   free(plde);
   free(sol);
   free(ind);
   free(p); 

   return 0;
};


