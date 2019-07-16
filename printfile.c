/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - printing PDB files
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
*****************************************************************************/

#include "bp.h"

int cut = 50;

int printpdb(int n,SOLUTION *sol,double lde,char *filename,int s)
{
   int i,k,nn;
   char file[1000],pdb[1000];
   FILE *output;

   // removing the path from the file name
   k = strlen(filename) - 1;
   while (k > 0 && filename[k] != '/')  k--;
   if (k == 0)  sprintf(file,"%s",filename);
   if (k != 0)  strncpy(file,&filename[k+1],strlen(filename)-k);
   file[strlen(filename)-k] = '\0';

   // removing extensions from the file name
   k = 0;
   while (k < strlen(file) && file[k] != '.')  k++;
   if (k != strlen(file)-1)  strncpy(pdb,&file[0],k);
   pdb[k] = '\0';

   // creating new file name
   if (s == 0)   // only the best solution is printed
   {
      sprintf(file,"%s_best.pdb",pdb);
   }
   else
   {
      sprintf(file,"%s_%d.pdb",pdb,s);
   };

   // opening the pdb file
   output = fopen(file,"w");
   if (output == NULL)  return 1;

   // writing the file
   fprintf(output,"HEADER      MD-jeep version 0.1\n");
   fprintf(output,"REMARK   1  Branch and Prune for the DMDGP\n");
   fprintf(output,"REMARK   1  by:\n");
   fprintf(output,"REMARK   1    A. Mucherino, L. Liberti, C. Lavor, N. Maculan\n");
   fprintf(output,"REMARK   1 \n");
   fprintf(output,"REMARK   1  filename: '%s'\n",file);
   fprintf(output,"REMARK   1 \n");
   if (s == 0)  fprintf(output,"REMARK   1  This is the best found solution\n");
   if (s != 0)  fprintf(output,"REMARK   1  Solution number = %d\n",s);
   fprintf(output,"REMARK   1  LDE function = %g\n",lde);
   fprintf(output,"REMARK   1 \n");

   // printing the choices (left/right,+/-) on the binary tree for this solution 
   nn = n/cut;
   fprintf(output,"REMARK   2  Sequence of choices on the binary tree:\n");
   fprintf(output,"REMARK   2 \n");
   for (k = 0; k < nn+1; k++)
   {
      fprintf(output,"REMARK   2   ");
      for (i = 0; i < cut && (k+1)*i < n; i++) 
      {
         if (sol[(k+1)*i].branch == -1) 
         {
            fprintf(output,"+");
         }
         else
         {
            fprintf(output,"-");
         };
      };
      fprintf(output,"\n");
   };
   fprintf(output,"REMARK   2 \n");
   
   // writing the PDB file
   k = 0;
   for (i = 0; i < n; i++)
   {
      if (sol[i].atom[0] == 'N')  k++;   // if N is not present in the instances, this amino acid counter cannot work ...
      fprintf(output,"ATOM   %4d  %-3s %-3s  %4d    %8.3lf%8.3lf%8.3lf \n",i+1,sol[i].atom,sol[i].amino,k,sol[i].x,sol[i].y,sol[i].z);
   };

   // closing the file
   fclose(output);

   return 0;
};
 

