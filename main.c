/****************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - main program
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    method based on change of basis added for coordinate computation
                                    more robust check on input instances added
              Jul 28 2019  v.0.3.0  main adapted for BP with interval distances
                                    more efficient organization of the data structures
              Mar 21 2020  v.0.3.1  using new functions of "vertex" to verify the instance properties
                                    precomputing all triplets of reference vertices
*****************************************************************************************************/

#include "bp.h"

double INFTY = 1.e+30;

int main(int argc, char *argv[])
{
   int i,j;
   int I,J;
   int gi,gj;
   int n,n0,m;
   int eof;
   int fidx;
   bool clique;
   bool consecutivity;
   double **X;
   double val1,val2;
   char ch[4][20];
   char *timestring;
   VERTEX *v;
   SEARCH S;
   OPTION op;
   INFORMATION info;
   struct timeval t1,t2;
   FILE *input;

   // welcome message
   fprintf(stderr,"MD-jeep 0.3.1\n");

   // checking input arguments
   if (argc < 2)
   {
      mdjeep_usage();
      return 1;
   };

   // default values for structures 'op' and 'info'
   op.version = 3;  op.print = 0;  op.format = 0;  op.allone = 0;  op.symmetry = 0;  
   op.monitor = true;  op.r = 1.0;  op.eps = 0.001;
   info.exact = false; info.ncalls = 0;  info.nspg = 0;  info.nspgok = 0; info.nsols = 0;
   info.maxsols = 10000;  info.pruning = 0;  
   info.best_sol = 0;  info.best_mde = INFTY;  info.best_lde = INFTY;
   consecutivity = false;

   // checking input options
   fidx = 1;
   while (fidx < argc - 1) 
   {
      if (!strcmp(argv[fidx],"-nomonitor"))
      {
         op.monitor = false;
         fidx++;
      }
      else if (!strcmp(argv[fidx],"-v"))
      {
         if (fidx + 1 >= argc - 1)
         {
            fprintf(stderr,"mdjeep: error: -v flag requires a version number (eg. 0.1)\n");
            return 1;
         };
         if (!strcmp(argv[fidx+1],"0.1") || !strcmp(argv[fidx+1],"0.2"))
         {
            op.version = 0;
         };
         fidx = fidx + 2;
      }
      else if (!strcmp(argv[fidx],"-e"))
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
      else if (!strcmp(argv[fidx],"-r"))
      {
         if (fidx + 1 >= argc - 1)
         {
            fprintf(stderr,"mdjeep: error: -r flag requires a real argument\n");
            return 1;
         };
         op.r = atof(argv[fidx+1]);
         if (op.r <= 0 || op.r > 10.0)
         {
            fprintf(stderr,"mdjeep: error: %lf is not an acceptable resolution value\n",op.r);
            return 1;
         };
         fidx = fidx + 2;
      }
      else if (!strcmp(argv[fidx],"-1"))
      {
         op.allone = 1;
         fidx++;
      }
      else if (!strcmp(argv[fidx],"-sym"))
      {
         if (fidx + 1 >= argc - 1)
         {
            fprintf(stderr,"mdjeep: error: -sym flag requires an integer argument (1=first half of the tree; 2=second half)\n");
            return 1;
         };
         if (atoi(argv[fidx+1]) == 1)
         {
            op.symmetry = 1;
         }
         else
         {
            op.symmetry = 2;
         };
         fidx = fidx + 2;
      }
      else if (!strcmp(argv[fidx],"-p"))
      {
         op.print = 1;
         fidx++;
      }
      else if (!strcmp(argv[fidx],"-P"))
      {
         op.print = 2;
         fidx++;
      }
      else if (!strcmp(argv[fidx],"-f"))
      {
         if (fidx + 1 >= argc - 1)
         {
            fprintf(stderr,"mdjeep: error: -f flag requires a char string (either xyz or pdb)\n");
            return 1;
         };
         if (!strcmp(argv[fidx+1],"pdb") || !strcmp(argv[fidx+1],"PDB"))  op.format = 1;
         fidx = fidx + 2;
      }
      else if (!strcmp(argv[fidx],"-consec"))
      {
         consecutivity = true;
         fidx++;
      }
      else
      {
         fprintf(stderr,"mdjeep: error: unknown option (%s)\n",argv[fidx]);
         return 1;
      };
   };

   // checking if there is an input file
   if (fidx >= argc) 
   {
      fprintf(stderr,"mdjeep: error: no input filename found on command line\n");
      return 1;
   };

   // additional information is printed on the screen
   if (op.print == 1)  fprintf(stderr,"mdjeep: the best solution ");
   if (op.print == 2)  fprintf(stderr,"mdjeep: all solutions ");
   if (op.print != 0)
   {
      fprintf(stderr,"will be printed in ");
      if (op.format == 0)
         fprintf(stderr,"XYZ");
      else
         fprintf(stderr,"PDB");
      fprintf(stderr," format\n");
   };
   if (op.allone == 1)  fprintf(stderr,"mdjeep: only one solution is requested by the user\n");
   if (op.symmetry != 0)  fprintf(stderr,"mdjeep: only one symmetric half of the tree is explored: ");
   if (op.symmetry == 1)  fprintf(stderr,"left-side subtree\n");
   if (op.symmetry == 2)  fprintf(stderr,"right-side subtree\n");
   fprintf(stderr,"mdjeep: tolerance epsilon = %g, resolution = %5.2lf\n",op.eps,op.r);
   if (op.version != 3)  fprintf(stderr,"mdjeep: the previous file format of versions 0.1 and 0.2 is used for reading input file\n");

   // opening input file
   info.filename = (char*)calloc(strlen(argv[fidx]),sizeof(char));
   strcpy(info.filename,argv[fidx]);
   input = fopen(argv[fidx],"r");
   if (input == NULL)
   {
      fprintf(stderr,"mdjeep: error while reading input file '%s'\n",argv[fidx]);
      return 1;
   };

   // computing the number of vertices and distances in the input file
   n = 0;  m = 0;  n0 = -1;  eof = 1;
   while (eof != EOF)
   {
      if (op.version == 3)
         eof = fscanf(input,"%d %d %d %d %lf %lf %s %s %s %s\n",&i,&j,&gi,&gj,&val1,&val2,ch[0],ch[1],ch[2],ch[3]);
      else
         eof = fscanf(input,"%d %d %lf %lf %s %s %s %s\n",&i,&j,&val1,&val2,ch[0],ch[1],ch[2],ch[3]);
      if (eof != EOF)
      {
         m++;
         if (i > n)  n = i;
         if (j > n)  n = j;
         if (n0 == -1)
            n0 = n;
         else if (n0 > i)
            n0 = i;
         else if (n0 > j)
            n0 = j;
      };
   };
   rewind(input);

   // initial check on the input file
   if (n <= 0 || m <= 0)
   {
      fprintf(stderr,"mdjeep: error: the input file seems to be empty\n");
      return 1;
   };

   // memory allocation for the input instance
   v = (VERTEX*)calloc(n-n0+1,sizeof(VERTEX));
   for (i = 0; i < n; i++)  v[i].Id = -1;

   // loading the instance
   eof = 1;
   while (eof != EOF)
   {
      // reading next line
      if (op.version == 3)
      {
         eof = fscanf(input,"%d %d %d %d %lf %lf %s %s %s %s\n",&i,&j,&gi,&gj,&val1,&val2,ch[0],ch[1],ch[2],ch[3]);
         if (gi < 0)  gi = 0;
         if (gj < 0)  gj = 0;
      }
      else
      {
         eof = fscanf(input,"%d %d %lf %lf %s %s %s %s\n",&i,&j,&val1,&val2,ch[0],ch[1],ch[2],ch[3]);
         gi = 0;  gj = 0;
      };
      if (eof != EOF)
      {
         // checking the vertex integer label
         if (i < 0 || j < 0)
         {
            fprintf(stderr,"mdjeep: error while reading input file, looks like there is a negative integer label\n");
            return 1;
         };

         // loading the first vertex
         // (and verifying that there are no two vertices with the same Id)
         I = i - n0;
         if (v[I].Id == -1)
         {
            initVertex(&v[I],i,gi,ch[0],ch[2]);
         }
         else
         {
            if (v[I].groupId != gi || strcmp(v[I].Name,ch[0]) != 0 || strcmp(v[I].Group,ch[2]) != 0)
            {
               fprintf(stderr,"mdjeep: error while reading input file, two different vertices with the same integer label\n");
               fprintf(stderr,"        vertex1 (%s,%s%d) <> vertex2 (%s,%s)\n",v[I].Name,v[I].Group,v[I].groupId,ch[0],ch[2]);
               return 1;
            };
         };

         // loading the second vertex
         // (and verifying that there are no two vertices with the same Id)
         J = j - n0;
         if (v[J].Id == -1)
         {
            initVertex(&v[J],j,gj,ch[1],ch[3]);
         }
         else
         {
            if (v[J].groupId != gj || strcmp(v[J].Name,ch[1]) != 0 || strcmp(v[J].Group,ch[3]) != 0)
            {
               fprintf(stderr,"mdjeep: error while reading input file, two different vertices with the same integer label\n");
               fprintf(stderr,"        vertex1 (%s,%s%d) <> vertex2 (%s,%s)\n",v[J].Name,v[J].Group,v[J].groupId,ch[1],ch[3]);
               return 1;
            };
         };

         // loading the distance bounds
         if (i < j)
         {
            I = i - n0;  J = j - n0;
         }
         else if (j < i)
         {
            I = j - n0;  J = i - n0;
         }
         else
         {
            fprintf(stderr,"mdjeep: error while reading input file, it looks like the vertex distance to itself is in the list\n");
            return 1;
         };

         if (v[J].ref == NULL)
         {
            v[J].ref = initReference(I,val1,val2);
         }
         else
         {
            addDistance(v[J].ref,I,val1,val2);
         };
      };
   };
   fclose(input);

   // verification of instance consistency
   if (m != totalNumberOfDistances(n-n0+1,v))
   {
      fprintf(stderr,"mdjeep: error while reading the instance; distance set not correctly loaded\n");
      return 1;
   };

   // printing instance details
   fprintf(stderr,"mdjeep: file '%s' read: %d vertices / %d distances\n",argv[fidx],n-n0+1,m);

   // checking whether the first three instance vertices form a clique
   clique = initialClique(n-n0+1,v,op.eps);
   if (!clique)
   {
      fprintf(stderr,"mdjeep: error: the first three vertices of the input instance do not form a clique\n");
      fprintf(stderr,"               the instance cannot be discretized\n");
      return 1;
   };

   // checking whether the input instance is discretizable
   i = isDDGP(n-n0+1,v,op.eps,clique);
   if (i != 0)
   {
      fprintf(stderr,"mdjeep: error: the input instance is not discretizable\n");
      fprintf(stderr,"               not enough references for vertex %d (at least 3, at least 2 exact)\n",n0+i);
      fprintf(stderr,"               stopping here... other necessary distances may be unavailable\n");
      return 1;
   };

   // the input instance is discretizable: message
   fprintf(stderr,"mdjeep: the input instance is discretizable\n");

   // verifying whether all distances are exact
   if (m == totalNumberOfExactDistances(n-n0+1,v,op.eps))
   {
      info.exact = true;
      op.r = 0.0;
      fprintf(stderr,"mdjeep: the instance contains only exact distances (the resolution parameter is disabled)\n");
   };

   // checking the consecutivity assumption (optional)
   if (info.exact || consecutivity)
   {
      consecutivity = isDMDGP(n-n0+1,v,op.eps,true);
      fprintf(stderr,"mdjeep: the instance ");
      if (consecutivity)
         fprintf(stderr,"satisfies ");
      else
         fprintf(stderr,"does not satisfy ");
      fprintf(stderr,"the consecutivity assumption\n");
   };

   // message about memory allocation
   fprintf(stderr,"mdjeep: allocating memory ...");

   // memory allocation for the solution X
   X = allocateMatrix(3,n-n0+1);

   // setting up some extra parameters
   S.be = 0.05;  // bound expansion factor
   S.pi = 3.14159265358979323846;

   // memory allocation for the arrays in SEARCH
   S.sym = (bool*)calloc(n-n0+1,sizeof(bool));
   S.refs = (triplet*)calloc(n-n0+1,sizeof(triplet));
   // S.costheta = allocateVector(n-n0+1);  S.cosomega = allocateVector(n-n0+1);  // not used in versions 0.3.0 and 0.3.1
   S.e = (int*)calloc(n-n0+1,sizeof(int));
   S.pX = allocateMatrix(3,n-n0+1);      S.lX = allocateMatrix(3,n-n0+1);      S.uX = allocateMatrix(3,n-n0+1);
   S.y = allocateVector(m);              S.gy = allocateVector(m);             S.sy = allocateVector(m);
   S.yp = allocateVector(m);             S.gyp = allocateVector(m);
   S.gX = allocateMatrix(3,n-n0+1);      S.sX = allocateMatrix(3,n-n0+1);
   S.Xp = allocateMatrix(3,n-n0+1);      S.gXp = allocateMatrix(3,n-n0+1);
   S.DX = allocateMatrix(3,n-n0+1);      S.YX = allocateMatrix(3,n-n0+1);      S.ZX = allocateMatrix(3,n-n0+1);
   S.Dy = allocateVector(m);             S.Yy = allocateVector(m);             S.Zy = allocateVector(m);
   S.memory = allocateVector(n-n0+1);
   fprintf(stderr," done\n");

   // precomputing all reference vertices
   for (i = 3; i < n-n0+1; i++)  S.refs[i] = findReferences(n,v,i,op.eps,info.exact);

   // looking for symmetries in the search tree
   fprintf(stderr,"mdjeep: checking symmetries ... ");
   findSymmetries(n-n0+1,v,S.sym);
   fprintf(stderr,"layers:");
   for (i = 0; i < n-n0+1; i++)  if (S.sym[i])  fprintf(stderr," %d",n0+i);
   fprintf(stderr,"\n");

   // counting the maximum number of digits necessary to represent the vertex ranks
   // (necessary for monitor)
   info.ndigits = numberOfDigits(n);

   // removing extension from input file
   info.output = removExtension(info.filename);

   // calling BP
   fprintf(stderr,"mdjeep: exploring the search tree ... ");
   if (op.monitor)
   {
      fprintf(stderr,"layer ");
      for (i = 0; i < info.ndigits; i++)  fprintf(stderr," ");
   };
   gettimeofday(&t1,0);
   if (info.exact)
   {
      bp_exact(0,n,m,v,X,S,op,&info);
      if (info.nsols == 0)
      {
         fprintf(stderr,"\nmdjeep: first attempt failed ... trying again with optimized reference vertices ... ");
         for (i = 3; i < n-n0+1; i++)  S.refs[i] = findReferences(n,v,i,op.eps,false);
         if (op.monitor)
         {
            fprintf(stderr,"layer ");
            for (i = 0; i < info.ndigits; i++)  fprintf(stderr," ");
         };
         bp_exact(0,n,m,v,X,S,op,&info);
      };
   }
   else
   {
      bp(0,n,m,v,X,S,op,&info);
   };
   gettimeofday(&t2,0);
   fprintf(stderr,"\n");
 
   // printing the results
   fprintf(stderr,"mdjeep: %d solutions found ",info.nsols);
   if (info.nsols == info.maxsols)  fprintf(stderr,"(max %d)",info.maxsols);
   fprintf(stderr,"\n");
   fprintf(stderr,"mdjeep: %d branches were pruned\n",info.pruning);
   if (!info.exact)  fprintf(stderr,"mdjeep: %d calls to spectral projected gradient (%d successful)\n",info.nspg,info.nspgok);
   if (info.nsols > 0)  fprintf(stderr,"mdjeep: best solution #%d: LDE = %10.8lf, MDE = %10.8lf\n",info.best_sol,info.best_lde,info.best_mde);
   timestring = splitime(t1,t2);
   fprintf(stderr,"mdjeep: CPU time = %s\n",timestring);

   // freeing memory
   free(timestring);
   freeVector(S.memory);
   freeVector(S.Dy);  freeVector(S.Yy);  freeVector(S.Zy);
   freeMatrix(3,S.DX);  freeMatrix(3,S.YX);  freeMatrix(3,S.ZX);
   freeMatrix(3,S.Xp);  freeMatrix(3,S.gXp);
   freeMatrix(3,S.gX);  freeMatrix(3,S.sX);
   freeVector(S.yp);  freeVector(S.gyp);
   freeVector(S.y);  freeVector(S.gy);  freeVector(S.sy);
   freeMatrix(3,S.pX);  freeMatrix(3,S.lX);  freeMatrix(3,S.uX);
   free(S.e);
   // freeVector(S.costheta);  freeVector(S.cosomega);
   free(S.refs);
   free(S.sym);
   free(X);
   freeVertex(n,v);
   free(info.filename);

   return 0;
};

