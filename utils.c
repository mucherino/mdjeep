/*************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - utilities
  Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.2
  History:    May 01 2010  v.0.1  first release
              May 08 2014  v.0.2  functions costheta and cosomega updated
                                  functions for computing cosines (omega,theta) from coords added
                                  functions norm and cross_prod added
**************************************************************************************************/

#include "bp.h"

// computing the cosine of theta (angle among three atoms)
// -------------------------------------------------------

double costheta(int i,int j,int k,PROBL *p,int **ind)
{
   int k12,k23,k13;
   double d12,d23,d13;
   double val;

   k12 = ind[i][j];  k23 = ind[j][k];  k13 = ind[i][k];
   if (k12 != -1 && k23 != -1 && k13 != -1)
   {
      d12 = p[k12].l;  d23 = p[k23].l;  d13 = p[k13].l;
      val = d12*d12 + d23*d23 - d13*d13;
      val = val / (2.0*d12*d23);
      if (val < -1.0)  val = -1.0;
      if (val >  1.0)  val =  1.0;
   }
   else
   {
      val = -2.0;
   };

   return val;
};

double costheta_coords(int i,int j,int k,SOLUTION *sol,PROBL *p,int **ind)
{
   int k12,k23,k13;
   double d12,d23,d13;
   double val;

   k12 = ind[i][j];  k23 = ind[j][k];  k13 = ind[i][k];

   if (k12 != -1)
   {
      d12 = p[k12].l;
   }
   else
   {
      d12 = distance(sol[j-1].x,sol[j-1].y,sol[j-1].z,sol[i-1].x,sol[i-1].y,sol[i-1].z);
   };

   d23 = p[k23].l;  d13 = p[k13].l;  // known because of discr. assumptions

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
   int k12,k13,k14,k23,k24,k34;
   double d12,d13,d14,d23,d24,d34;
   double a,b,c,e,f;
   double d24q,d23q;
   double sumq,den,val;
   double sqrtarg;

   k12 = ind[i][j];  k13 = ind[i][k];  k14 = ind[i][h];
   k23 = ind[j][k];  k24 = ind[j][h];  k34 = ind[k][h];

   if (k12 != -1 && k13 != -1 && k14 != -1 && k23 != -1 && k24 != -1 && k34 != -1)
   {
      d12 = p[k12].l;  d13 = p[k13].l;  d14 = p[k14].l;
      d23 = p[k23].l;  d24 = p[k24].l;  d34 = p[k34].l;

      a = d12*d12 + d24*d24 - d14*d14;        a = a / (2.0*d12*d24);
      b = d24*d24 + d23*d23 - d34*d34;        b = b / (2.0*d24*d23);
      c = d12*d12 + d23*d23 - d13*d13;        c = c / (2.0*d12*d23);
      e = 1.0 - b*b;
      f = 1.0 - c*c;
      if (e < 0.0 || f < 0.0)  return -3.0;  // something wrong!
      e = sqrt(e);  f = sqrt(f);
      val = (a - b*c) / (e*f);
      if (val < -1.0)  val = -1.0;
      if (val >  1.0)  val =  1.0;
   }
   else
   {
      val = -2.0;
   };

   return val;
};

double cosomega_coords(int i,int j,int k,int h,SOLUTION *sol,PROBL *p,int **ind)
{
   int k12,k13,k14,k23,k24,k34;
   double d12,d13,d14,d23,d24,d34;
   double a,b,c,e,f;
   double d24q,d23q;
   double sumq,den,val;
   double sqrtarg;

   k12 = ind[i][j];  k13 = ind[i][k];  k23 = ind[j][k];  // may not be known a priori
   k14 = ind[i][h];  k24 = ind[j][h];  k34 = ind[k][h];  // known because of discr. ass.

   if (k12 != -1)
   {
      d12 = p[k12].l;
   }
   else
   {
      d12 = distance(sol[j-1].x,sol[j-1].y,sol[j-1].z,sol[i-1].x,sol[i-1].y,sol[i-1].z);
   };

   if (k13 != -1)
   {
      d13 = p[k13].l;
   }
   else
   {
      d13 = distance(sol[k-1].x,sol[k-1].y,sol[k-1].z,sol[i-1].x,sol[i-1].y,sol[i-1].z);
   };

   if (k23 != -1)
   {
      d23 = p[k23].l;
   }
   else
   {
      d23 = distance(sol[k-1].x,sol[k-1].y,sol[k-1].z,sol[j-1].x,sol[j-1].y,sol[j-1].z);
   };
 
   d14 = p[k14].l;  d24 = p[k24].l;  d34 = p[k34].l;

   a = d12*d12 + d24*d24 - d14*d14;        a = a / (2.0*d12*d24);
   b = d24*d24 + d23*d23 - d34*d34;        b = b / (2.0*d24*d23);
   c = d12*d12 + d23*d23 - d13*d13;        c = c / (2.0*d12*d23);
   e = 1.0 - b*b;
   f = 1.0 - c*c;
   if (e < 0.0 || f < 0.0)  return -3.0;  // something wrong!
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

// computing the norm of a 3-dimensional vector
// --------------------------------------------

double norm(double *v)
{
   return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
};

// computing the cross product between two 3-dimensional vectors
// -------------------------------------------------------------

void cross_prod(double *v1,double *v2,double *res)
{
   res[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]); 
   res[1] = (v1[2]*v2[0]) - (v1[0]*v2[2]);
   res[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
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
   fprintf(stderr,"    -m1 method1 used for coordinate computation (based on cumulative matrices)\n");
   fprintf(stderr,"    -m2 method2 used for coordinate computation (based on change of basis)\n");
   fprintf(stderr,"     -e sets the tolerance epsilon\n");
   fprintf(stderr,"     -p prints the best found solution in a PDB file\n");
   fprintf(stderr,"     -P prints all the solutions in different PDB files\n");
   fprintf(stderr,"     -s only one symmetric half of the tree is explored\n");
   fprintf(stderr,"     -t triangular inequalities are checked\n");
   fprintf(stderr,"     -1 the algorithm stops at the first solution\n");
   fprintf(stderr,"        (when using -1, options -p and -P have the same effect)\n");
};

