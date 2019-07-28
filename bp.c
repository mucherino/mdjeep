/********************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - BP
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    implementation of 2nd method for coordinate computation
              Jul 28 2019  v.0.3.0  extension to interval distances, use of coarse-grained representation, 
                                    and resolution parameter
*********************************************************************************************************/

#include "bp.h"

extern double INFTY;
bool keep_going = true;
bool PRINTED = false;

// signal catcher
void intHandler(int a)
{
   fprintf(stderr," signal caught: stopping (partial solution printed if -P or -p options used)");
   keep_going = false;
};

// branch-and-prune
void bp(int i,int n,int m,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info)
{
   int j;
   int nb,it;
   int pru;
   double A,B,U[9];
   double cdist;
   double ctheta,stheta;
   double cosomega00,cosomega01;
   double sinomega00,sinomega01;
   double lomega0,uomega0;
   double lomega1,uomega1;
   double omega;
   double dist,alpha,opt;
   double obj,lde,mde;
   Omega *current;
   omegaList omegaL;
   REFERENCE *rtmp,*r1,*r2,*r3;

   // signal handler
   signal(SIGINT,intHandler);

   // first call to BP?
   if (i == 0)
   {
      // initializing the BP call counter
      info->ncalls = 0;

      // The first three vertices can be positioned by using the initial clique

      // vertex 0
         X[0][0] = 0.0;         X[1][0] = 0.0;         X[2][0] = 0.0;
      S.lX[0][0] = -op.eps;  S.lX[1][0] = -op.eps;  S.lX[2][0] = -op.eps;
      S.uX[0][0] =  op.eps;  S.uX[1][0] =  op.eps;  S.uX[2][0] =  op.eps;

      // vertex 1
      r1 = getReference(v,0,1);
         X[0][1] = -lowerBound(r1);       X[1][1] = 0.0;         X[2][1] = 0.0;
      S.lX[0][1] =  X[0][1] - op.eps;  S.lX[1][1] = -op.eps;  S.lX[2][1] = -op.eps; 
      S.uX[0][1] =  X[0][1] + op.eps;  S.uX[1][1] =  op.eps;  S.uX[2][1] =  op.eps;

      // vertex 2
      r2 = getReference(v,1,2);
      ctheta = costheta(0,1,2,v,X);  stheta = sqrt(1.0 - ctheta*ctheta);
         X[0][2] = -lowerBound(r1) + lowerBound(r2)*ctheta;  X[1][2] = lowerBound(r2)*stheta;  X[2][2] = 0.0; 
      S.lX[0][2] = X[0][2] - op.eps;  S.lX[1][2] = X[1][2] - op.eps;  S.lX[2][2] = -op.eps;
      S.uX[0][2] = X[0][2] + op.eps;  S.uX[1][2] = X[1][2] + op.eps;  S.uX[2][2] =  op.eps;

      i = i + 3;  // branching starts at vertex i+3
   };

   // updating BP call counter
   info->ncalls++;

   // searching for the references
   r1 = v[i].ref;
   if (isIntervalDistance(r1,op.eps))  r1 = nextExactDistance(v[i].ref,op.eps);
   rtmp = nextExactDistance(r1,op.eps);
   if (rtmp == NULL || r1 == NULL)
   {
      fprintf(stderr,"BP internal error: it looks like the discretization assumptions are actually not satisfied\n");
      abort();
   };
   if (r1->otherId > rtmp->otherId)
   {
      r2 = rtmp;
   }
   else
   {
      r2 = r1;  r1 = rtmp;
   };
   rtmp = nextExactDistance(rtmp,op.eps);
   if (rtmp != NULL)
   {
      if (r2->otherId > rtmp->otherId)
      {
         r3 = rtmp;
      }
      else if (r1->otherId > rtmp->otherId)
      {
         r3 = r2;  r2 = rtmp;
      }
      else
      {
         r3 = r2;  r2 = r1;  r1 = rtmp;
      };

      rtmp = nextExactDistance(rtmp,op.eps);
      while (rtmp != NULL)
      {
         if (rtmp->otherId > r1->otherId)
         {
            r3 = r2;  r2 = r1;  r1 = rtmp;
         }
         else if (rtmp->otherId > r2->otherId)
         {
            r3 = r2;  r2 = rtmp;
         }
         else if (rtmp->otherId > r1->otherId)
         {
            r3 = rtmp;
         };
         rtmp = nextExactDistance(rtmp,op.eps);
      };
   }
   else
   {
      if (isIntervalDistance(v[i].ref,op.eps))
         r3 = v[i].ref;
      else
         r3 = nextIntervalDistance(v[i].ref,op.eps);
      if (r3 == NULL)
      {
         fprintf(stderr,"BP internal error: it looks like the discretization assumptions are actually not satisfied\n");
         abort();
      };
      rtmp = nextIntervalDistance(r3,op.eps);
      while (rtmp != NULL)
      {
         if (rtmp->otherId > r3->otherId)  r3 = rtmp;
         rtmp = nextIntervalDistance(rtmp,op.eps);
      };
   };
   cdist = lowerBound(r1);
   
   // theta angles ("bond" angles)
   ctheta = costheta(r2->otherId,r1->otherId,i,v,X);
   stheta = sqrt(1.0 - ctheta*ctheta);

   // generating U matrix (only once)
   UMatrix(r3->otherId,r2->otherId,r1->otherId,i,X,U);

   // omega angles (torsion angles)
   nb = 2;
   cosomega00 = cosomega(r3->otherId,r2->otherId,r1->otherId,i,v,X,0.0);
   cosomega01 = cosomega(r3->otherId,r2->otherId,r1->otherId,i,v,X,1.0);
   sinomega00 = sqrt(1.0 - cosomega00*cosomega00);
   sinomega01 = sqrt(1.0 - cosomega01*cosomega01);
   lomega0 = atan2(+sinomega00,cosomega00);  uomega0 = atan2(+sinomega01,cosomega01);
   lomega1 = atan2(-sinomega00,cosomega00);  uomega1 = atan2(-sinomega01,cosomega01);

   // if the two omega intervals are adjacent, we can consider the union
   if (fabs(uomega0 - lomega1) < op.eps)
   {
      nb = 1;
      uomega0 = uomega1;
   }
   else if (fabs(uomega1 - lomega0) < op.eps)
   {
      nb = 1;
      lomega0 = lomega1;
   };

   // if the layer is symmetric, it is not necessary to consider the entire omega intervals
   if (S.sym[i])
   {
      lomega0 = 0.5*(lomega0 + uomega0);
      uomega0 = lomega0;
      if (nb == 2)
      {
         lomega1 = 0.5*(lomega1 + uomega1);
         uomega1 = lomega1;
      };
   };

   // initializing omega list
   omegaL = initOmegaList(lomega0,uomega0);
   if (nb == 2)  attachNewOmegaInterval(firstOmegaInterval(omegaL),lomega1,uomega1);

   // verifying the "arclength" of every arc wrt the given resolution parameter
   if (!info->exact)  splitOmegaIntervals(firstOmegaInterval(omegaL),cdist,op.r);

   // starting point for iterating over omega angles (it depends on op.symmetry)
   if (op.symmetry < 2)
      current = firstOmegaInterval(omegaL);
   else
      current = lastOmegaInterval(omegaL);

   // branching over the obtained omega sub-intervals (using omegaList iterators)
   while (current != NULL && keep_going)
   {
      // monitor
      if (op.monitor)  fprintf(stderr,"\b\b\b%3d",i);

      // the vertex position is initially placed at the center of the arc
      if (info->exact)
      {
         omega = omegaIntervalLowerBound(current);
      }
      else
      {
         lomega0 = omegaIntervalLowerBound(current);
         uomega0 = omegaIntervalUpperBound(current);
         omega = 0.5*(lomega0 + uomega0);
      };
      genCoordinates(r1->otherId,i,X,U,cdist,ctheta,stheta,cos(omega),sin(omega));

      // generation of the box inscribing the arc
      if (isExactDistance(r3,op.eps))
      {
         // the box has the size of the tolerance on the three dimensions
         S.lX[0][i] = X[0][i] - op.eps;  S.uX[0][i] = X[0][i] + op.eps;
         S.lX[1][i] = X[1][i] - op.eps;  S.uX[1][i] = X[1][i] + op.eps;
         S.lX[2][i] = X[2][i] - op.eps;  S.uX[2][i] = X[2][i] + op.eps;
      }
      else
      {
         // computing min and max x coordinates over the arc
         A = U[3]*cdist*stheta;  B = U[6]*cdist*stheta;
         if (A != 0.0)
         {
            alpha = atan2(B,A);  opt = S.pi + alpha;
            opt = projection(opt,lomega0,uomega0,op.eps);
            S.lX[0][i] = S.lX[0][r1->otherId] - U[0]*cdist*ctheta + A*cos(opt) + B*sin(opt) - op.eps;
            opt = projection(alpha,lomega0,uomega0,op.eps);
            S.uX[0][i] = S.uX[0][r1->otherId] - U[0]*cdist*ctheta + A*cos(opt) + B*sin(opt) + op.eps;
         }
         else
         {
            S.lX[0][i] = X[0][i] - op.eps;  S.uX[0][i] = X[0][i] + op.eps;
         };

         // computing min and max y coordinates over the arc
         A = U[4]*cdist*stheta;  B = U[7]*cdist*stheta;
         if (A != 0.0)
         {
            alpha = atan2(B,A);  opt = S.pi + alpha;
            opt = projection(opt,lomega0,uomega0,op.eps);
            S.lX[1][i] = S.lX[1][r1->otherId] - U[1]*cdist*ctheta + A*cos(opt) + B*sin(opt) - op.eps;
            opt = projection(alpha,lomega0,uomega0,op.eps);
            S.uX[1][i] = S.uX[1][r1->otherId] - U[1]*cdist*ctheta + A*cos(opt) + B*sin(opt) + op.eps;
         }
         else
         {
            S.lX[1][i] = X[1][i] - op.eps;  S.uX[1][i] = X[1][i] + op.eps;
         };

         // computing min and max z coordinates over the arc
         A = U[5]*cdist*stheta;  B = U[8]*cdist*stheta;
         if (A != 0.0)
         {
            alpha = atan2(B,A);  opt = S.pi + alpha;
            opt = projection(opt,lomega0,uomega0,op.eps);
            S.lX[2][i] = S.lX[2][r1->otherId] - U[2]*cdist*ctheta + A*cos(opt) + B*sin(opt) - op.eps;
            opt = projection(alpha,lomega0,uomega0,op.eps);
            S.uX[2][i] = S.uX[2][r1->otherId] - U[2]*cdist*ctheta + A*cos(opt) + B*sin(opt) + op.eps;
         }
         else
         {
            S.lX[2][i] = X[2][i] - op.eps;  S.uX[2][i] = X[2][i] + op.eps;
         };
      };

      // performing the DDF pruning device
      pru = DDF(i,v,X,op.eps);

      // improving the current solution by local optimization (if necessary)
      if (!info->exact && pru > 0)
      {
         // call the spectral projected gradient
         S.be = 0.5;
         while (pru > 0 && S.be < 5.0)
         {
            spg(i+1,v,X,S,&it,&obj);
            info->nspg++;
            pru = DDF(i,v,X,op.eps);
            S.be = 2.0*S.be;
         };
         if (pru == 0)  info->nspgok++;
      };

      if (pru > 0)  info->pruning++;

      // if the current partial solution is OK (either from the beginning, or after the SPG call)
      if (pru == 0)
      {
         // verifying whether the partial solution is too close to the previous computed one
         if (info->nsols > 0)
         {
            dist = 0.0;
            for (j = 0; j <= i; j++)
            {
               dist = dist + pairwise_distance(X[0][j],X[1][j],X[2][j],S.pX[0][j],S.pX[1][j],S.pX[2][j]);
            };
            dist = dist/(i+1);
            if (dist < op.r)
            {
               goto NEXT;  // skip this branch
            };
         };

         if (i < n - 1)  
         {
            // next vertex
            bp(i+1,n,m,v,X,S,op,info);
         }
         else
         {
            // solution found
            info->nsols = info->nsols + 1;

            // printing the solution (if requested)
            if (op.print > 1)
            {
               if (op.format == 0)
                  printfile(n,v,X,info->output,info->nsols);
               else
                  printpdb(n,v,X,info->output,info->nsols);
            };

            // evaluating the quality of the solution
            lde = compute_lde(n,v,X,op.eps);
            mde = compute_mde(n,v,X,op.eps);

            // best solution found so far
            if (mde < info->best_mde)
            {
               info->best_sol = info->nsols;
               info->best_lde = lde;
               info->best_mde = mde;
               if (op.print == 1)
               {
                  if (op.format == 0)
                     printfile(n,v,X,info->output,0);
                  else
                     printpdb(n,v,X,info->output,0);
               };
            };
            copyMatrix(3,n,X,S.pX);
         };
      };

      // if only one solution is requested, bp stops as soon as the first solution is found
      if (op.allone == 1)  if (info->nsols > 0)  break;

      // the search stops after maxsols solutions
      if (info->nsols >= info->maxsols)  break;

      // skipping one half of the tree (optional)
      if (i == 4)  if (op.symmetry > 0)  break;

      // preparing for next iteration
NEXT: if (omegaIntervalHasNextAlongDirection(current,op.symmetry<2))
         current = omegaIntervalNextAlongDirection(current,op.symmetry<2);
      else
         current = NULL;
   };

   // handling ^C signal catcher
   if (!keep_going)
   {
      if (!PRINTED)
      {
         if(op.print > 0 && info->nsols == 0)
         {
            if (op.format == 0)
               printfile(n,v,X,info->output,0);
            else
               printpdb(n,v,X,info->output,0);
            PRINTED = true;
         };
      };
   };

   // freeing memory space for omega list
   freeOmegaList(omegaL);

   return;
};

