/********************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - functions to print
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    binary representation of solutions corrected
              Jul 28 2019  v.0.3.0  generic function "printfile" added, some modifications on "printpdb"
              Mar 21 2020  v.0.3.1  no changes
              May 19 2020  v.0.3.2  undefined vertex attributes are not printed
*********************************************************************************************************/

#include "bp.h"

// print solutions with the available vertex attributes in a text file
// (s=0 : prints only one solution; solution number s>0 : prints multiple solutions in the same file)
void printfile(int n,VERTEX *v,double **X,char *filename,int s)
{
   int i,k;
   int groupId;
   char *outfile;
   bool group;
   FILE *output;

   // if n is not positive, there is nothing to print
   if (n <= 0)  return;

   // file name (txt extension)
   outfile = (char*)calloc(strlen(filename)+5,sizeof(char));
   for (k = 0; k < strlen(filename); k++)  outfile[k] = filename[k];
   outfile[k] = '.';  k++;
   outfile[k] = 't';  k++;
   outfile[k] = 'x';  k++;
   outfile[k] = 't';

   // opening output file (rewriting or appending)
   if (s < 2)
      output = fopen(outfile,"w");
   else
      output = fopen(outfile,"a");
   if (output == NULL)
   {
      fprintf(stderr,"printfile: error while opening file '%s' to write\n",filename);
      abort();
   };

   // groupId will not be printed if it always contain the same value
   group = false;  groupId = v[0].groupId;
   for (i = 1; i < n && !group; i++)  group = group || (groupId != v[i].groupId);
 
   // writing file
   if (s != 0)  fprintf(output,"MODEL %d\n",s);
   for (i = 0; i < n; i++)
   {
      // Id
      fprintf(output," %d",v[i].Id);

      // Name
      if (strcmp(v[i].Name,"(no name)"))  fprintf(output," %s",v[i].Name);

      // groupId
      if (group)  fprintf(output," %d",v[i].groupId);

      // Group
      if (strcmp(v[i].Group,"(no group name)"))  fprintf(output," %s",v[i].Group);

      // coordinates
      fprintf(output," %13.9lf %13.9lf %13.9lf\n",X[0][i],X[1][i],X[2][i]);
   };
   fclose(output);
   free(outfile);
};

// print solutions in PDB format
// (s=0 : prints only one solution; solution number s>0 : prints multiple solutions in the same file)
void printpdb(int n,VERTEX *v,double **X,char *filename,int s)
{
   int i,k;
   char *outfile;
   FILE *output;

   // if n is not positive, there is nothing to print
   if (n <= 0)  return;

   // file name (pdb extension)
   outfile = (char*)calloc(strlen(filename)+5,sizeof(char));
   for (k = 0; k < strlen(filename); k++)  outfile[k] = filename[k];
   outfile[k] = '.';  k++;
   outfile[k] = 'p';  k++;
   outfile[k] = 'd';  k++;
   outfile[k] = 'b';

   // opening output file (rewriting or appending)
   if (s < 2)
      output = fopen(outfile,"w");
   else
      output = fopen(outfile,"a");
   if (output == NULL)
   {
      fprintf(stderr,"printfile: error while opening file '%s' to write\n",filename);
      abort();
   };

   // writing header file (only if new file or overwriting)
   if (s < 2)
   {
      fprintf(output,"HEADER      MD-jeep version 0.3.2\n");
      fprintf(output,"REMARK   1 \n");
      fprintf(output,"REMARK   1  Branch and Prune for Discretizable Distance Geometry\n");
      fprintf(output,"REMARK   1 \n");
      fprintf(output,"REMARK   1  by: Mucherino, Goncalves, Lavor, Liberti, Lin, Maculan\n");
      fprintf(output,"REMARK   1 \n");
      fprintf(output,"REMARK   1  filename: '%s'\n",outfile);
      fprintf(output,"REMARK   1 \n");
   };

   // writing the PDB model
   if (s != 0)  fprintf(output,"MODEL%9d\n",s);
   for (i = 0; i < n; i++)
   {
      // Id
      fprintf(output,"%-6s%5d  ","ATOM",v[i].Id);

      // atom code
      if (strcmp(v[i].Name,"(no name)"))
         fprintf(output,"%-4s",v[i].Name);
      else
         fprintf(output,"%-4s","XX");

      // amino acid code
      if (strcmp(v[i].Group,"(no group name)"))
         fprintf(output,"%-3s ",v[i].Group);
      else
         fprintf(output,"%-3s ","UNK");

      // other info, including coordinates
      fprintf(output,"%s%4d    %8.3f%8.3f%8.3f \n","A",v[i].groupId,X[0][i],X[1][i],X[2][i]);
   };
   if (s != 0)  fprintf(output,"ENDMDL%8d\n",s);

   // ending
   fclose(output);
   free(outfile);
};

