/****************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for the DMDGP - utilities
  Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
*****************************************************************************/

#include "bp.h"

// computing the cosine of theta (angle among three atoms)
// -------------------------------------------------------

double costheta(int i,int j,int k,PROBL *p,int **ind)
{
   double d12,d23,d13;
   double val;

   d12 = p[ind[i][j]].l;  d23 = p[ind[j][k]].l;  d13 = p[ind[i][k]].l;
   val = d12*d12 + d23*d23 - d13*d13;
   val = val / (2.0*d12*d23);

   if (val < -1.0)  val = -1.0;
   if (val >  1.0)  val =  1.0;

   return val;
};


// computing the cosine of the torsion angles
// ------------------------------------------

double cosomega(int i,int j,int k,int h,PROBL *p,int **ind)
{
   double d12,d13,d14,d23,d24,d34;
   double a,b,c,e,f;
   double d24q,d23q;
   double sumq,den,val;
   double sqrtarg;

   d12 = p[ind[i][j]].l;  d13 = p[ind[i][k]].l;  d14 = p[ind[i][h]].l;
   d23 = p[ind[j][k]].l;  d24 = p[ind[j][h]].l;  d34 = p[ind[k][h]].l;

   a = d12*d12 + d24*d24 - d14*d14;        a = a / (2.0*d12*d24);
   b = d24*d24 + d23*d23 - d34*d34;        b = b / (2.0*d24*d23);
   c = d12*d12 + d23*d23 - d13*d13;        c = c / (2.0*d12*d23);
   e = 1.0 - b*b;
   f = 1.0 - c*c;
   if (e < 0.0 || f < 0.0)  return -2;  // this is the signal that something went wrong during the computation of the cosine
   e = sqrt(e);  f = sqrt(f);
   val = (a - b*c) / (e*f);

   if (val < -1.0)  val = -1.0;
   if (val >  1.0)  val =  1.0;

   return val;
};


// computing the distance between two points
// -----------------------------------------

double distance(double x1,double y1,double z1,double x2,double y2,double z2)
{
   double px,py,pz,dist;

   px = x2 - x1;
   py = y2 - y1;
   pz = z2 - z1;
   dist = sqrt( px*px + py*py + pz*pz );

   return dist;
};


// BP usage - function invoked if there are too few arguments
// ----------------------------------------------------------

void bp_usage(void)
{
   fprintf(stderr,"bp: too few arguments\n");
   fprintf(stderr,"    syntax: ./bp [options] instance.nmr\n");
   fprintf(stderr,"    NMR file format: list of distances [i j iatom jatom iamino jamino l u]\n");
   fprintf(stderr,"       where:\n");
   fprintf(stderr,"    * i is the label of the first atom\n");
   fprintf(stderr,"    * j is the label of the second atom\n");
   fprintf(stderr,"    * iatom is the name of the atom i (3 letter format)\n");
   fprintf(stderr,"    * jatom is the name of the atom j (3 letter format)\n");
   fprintf(stderr,"    * iamino is the amino acid where atom i is contained\n");
   fprintf(stderr,"    * jamino is the amino acid where atom j is contained\n");
   fprintf(stderr,"    * l is the lower bound on the distance (i,j)\n");
   fprintf(stderr,"    * u is the upper bound on the distance (i,j)\n");
   fprintf(stderr,"      (in this version, l and u must be the same)\n");
   fprintf(stderr,"    Options:\n");
   fprintf(stderr,"    -e sets the tolerance epsilon\n");
   fprintf(stderr,"    -p prints the best found solution in a PDB file\n");
   fprintf(stderr,"    -P prints all the solutions in different PDB files\n");
   fprintf(stderr,"    -s symmetric solutions are not computed\n");
   fprintf(stderr,"    -t triangular inequalities are checked\n");
   fprintf(stderr,"    -1 the algorithm stops at the first solution\n");
   fprintf(stderr,"       (when using -1, options -p and -P have the same effect)\n");
};



