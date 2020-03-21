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
              Mar 21 2020  v.0.3.1  improving method for the computation of the box bounds
                                    adding bp_exact function, for instances consisting of exact distances
*********************************************************************************************************/

#include "bp.h"

extern double INFTY;
bool keep_going = true;
bool check = false;
bool PRINTED = false;

// signal catcher
void intHandler(int a)
{
   fprintf(stderr," signal caught: stopping (partial solution printed if -P or -p options used)");
   keep_going = false;
};

// branch-and-prune (general version)
void bp(int i,int n,int m,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info)
{
   int j,k;
   int nb,it;
   int pru;
   int ldigits;
   double A,B,U[9];
   double cdist;
   double cTheta,sTheta;
   double cosOmega00,cosOmega01;
   double sinOmega00,sinOmega01;
   double lomega0,uomega0;
   double lomega1,uomega1;
   double omega;
   double dist,alpha,opt;
   double obj,lde,mde;
   Omega *current;
   omegaList omegaL;
   REFERENCE *r1,*r2,*r3;

   // signal handler
   signal(SIGINT,intHandler);

   // first call to BP?
   if (i == 0)
   {
      // initializing the BP call counter
      info->ncalls = 0;

      // The first three vertices can be positioned by using the initial clique

      // vertex 0
         X[0][0] = 0.0;     X[1][0] = 0.0;     X[2][0] = 0.0;
      S.lX[0][0] = 0.0;  S.lX[1][0] = 0.0;  S.lX[2][0] = 0.0;
      S.uX[0][0] = 0.0;  S.uX[1][0] = 0.0;  S.uX[2][0] = 0.0;
      S.e[0] = 0;

      // vertex 1
      r1 = getReference(v,0,1);
         X[0][1] = -lowerBound(r1);  X[1][1] = 0.0;     X[2][1] = 0.0;
      S.lX[0][1] =  X[0][1];      S.lX[1][1] = 0.0;  S.lX[2][1] = 0.0;
      S.uX[0][1] =  X[0][1];      S.uX[1][1] = 0.0;  S.uX[2][1] = 0.0;
      S.e[1] = 0;

      // vertex 2
      r2 = getReference(v,1,2);
      cTheta = costheta(0,1,2,v,X);  sTheta = sqrt(1.0 - cTheta*cTheta);
         X[0][2] = -lowerBound(r1) + lowerBound(r2)*cTheta;  X[1][2] = lowerBound(r2)*sTheta;  X[2][2] = 0.0; 
      S.lX[0][2] = X[0][2];  S.lX[1][2] = X[1][2];  S.lX[2][2] = 0.0;
      S.uX[0][2] = X[0][2];  S.uX[1][2] = X[1][2];  S.uX[2][2] = 0.0;
      S.e[2] = 0;

      i = i + 3;  // branching starts at vertex i+3
   };

   // updating BP call counter
   info->ncalls++;
   it = 0;

   // reference vertices
   r3 = S.refs[i].r3;  r2 = S.refs[i].r2;  r1 = S.refs[i].r1;
   cdist = lowerBound(r1);
   
   // theta angles ("bond" angles)
   cTheta = costheta(r2->otherId,r1->otherId,i,v,X);
   sTheta = sqrt(1.0 - cTheta*cTheta);

   // generating U matrix (only once)
   UMatrix(r3->otherId,r2->otherId,r1->otherId,i,X,U);

   // omega angles (torsion angles)
   nb = 2;
   cosOmega00 = cosomega(r3->otherId,r2->otherId,r1->otherId,i,v,X,0.0);
   cosOmega01 = cosomega(r3->otherId,r2->otherId,r1->otherId,i,v,X,1.0);
   sinOmega00 = sqrt(1.0 - cosOmega00*cosOmega00);
   sinOmega01 = sqrt(1.0 - cosOmega01*cosOmega01);
   lomega0 = atan2(+sinOmega00,cosOmega00);  uomega0 = atan2(+sinOmega01,cosOmega01);
   lomega1 = atan2(-sinOmega00,cosOmega00);  uomega1 = atan2(-sinOmega01,cosOmega01);

   // if the two omega intervals are adjacent, we can consider the union
   if (i > 3)
   {
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
   splitOmegaIntervals(firstOmegaInterval(omegaL),cdist,op.r);

   // counting total number of omega intervals (necessary only at layer 3)
   if (i == 3)  nb = numberOfOmegaIntervals(firstOmegaInterval(omegaL));

   // starting point for iterating over omega angles (it depends on op.symmetry)
   if (op.symmetry < 2)
      current = firstOmegaInterval(omegaL);
   else
      current = lastOmegaInterval(omegaL);

   // branching over the obtained omega sub-intervals (using omegaList iterators)
   while (current != NULL && keep_going)
   {
      // monitor
      if (op.monitor)  
      {
         ldigits = numberOfDigits(i);
         for (k = 0; k < info->ndigits; k++)  fprintf(stderr,"\b");
         for (k = 0; k < info->ndigits - ldigits; k++)  fprintf(stderr," ");
         fprintf(stderr,"%d",i);
      };
      it++;

      // disabling the comparison with previous solutions when we move 
      // from the left to the right side of the tree
      if (op.symmetry == 0)  if (i == 3)  if (check)  if (it == nb/2 + 1)  check = false;

      // the vertex position is initially placed at the center of the arc
      lomega0 = omegaIntervalLowerBound(current);
      uomega0 = omegaIntervalUpperBound(current);
      omega = 0.5*(lomega0 + uomega0);
      genCoordinates(r1->otherId,i,X,U,cdist,cTheta,sTheta,cos(omega),sin(omega));

      // generation of the box inscribing the arc
      S.e[i] = 0;
      if (isExactDistance(r3,op.eps))
      {
         // the box has the size of the tolerance on the three dimensions
         S.lX[0][i] = X[0][i];  S.uX[0][i] = X[0][i];
         S.lX[1][i] = X[1][i];  S.uX[1][i] = X[1][i];
         S.lX[2][i] = X[2][i];  S.uX[2][i] = X[2][i];
      }
      else
      {
         // computing min and max x coordinates over the arc
         A = U[3]*cdist*sTheta;  B = U[6]*cdist*sTheta;
         if (A != 0.0)
         {
            lomega1 = A*cos(lomega0) + B*sin(lomega0);
            uomega1 = A*cos(uomega0) + B*sin(uomega0);
            alpha = atan2(B,A);
            if (alpha < 0.0)
               opt = alpha + S.pi;
            else
               opt = alpha - S.pi;
            alpha = A*cos(alpha) + B*sin(alpha);
            alpha = maximum(alpha,lomega1,uomega1);
            opt = A*cos(opt) + B*sin(opt);
            opt = minimum(opt,lomega1,uomega1);
            S.lX[0][i] = S.lX[0][r1->otherId] - U[0]*cdist*cTheta + opt - op.eps;
            S.uX[0][i] = S.uX[0][r1->otherId] - U[0]*cdist*cTheta + alpha + op.eps;
         }
         else
         {
            S.lX[0][i] = X[0][i] - op.eps;  S.uX[0][i] = X[0][i] + op.eps;
         };

         // computing min and max y coordinates over the arc
         A = U[4]*cdist*sTheta;  B = U[7]*cdist*sTheta;
         if (A != 0.0)
         {
            lomega1 = A*cos(lomega0) + B*sin(lomega0);
            uomega1 = A*cos(uomega0) + B*sin(uomega0);
            alpha = atan2(B,A);
            if (alpha < 0.0)
               opt = alpha + S.pi;
            else
               opt = alpha - S.pi;
            alpha = A*cos(alpha) + B*sin(alpha);
            alpha = maximum(alpha,lomega1,uomega1);
            opt = A*cos(opt) + B*sin(opt);
            opt = minimum(opt,lomega1,uomega1);
            S.lX[1][i] = S.lX[1][r1->otherId] - U[1]*cdist*cTheta + opt - op.eps;
            S.uX[1][i] = S.uX[1][r1->otherId] - U[1]*cdist*cTheta + alpha + op.eps;
         }
         else
         {
            S.lX[1][i] = X[1][i] - op.eps;  S.uX[1][i] = X[1][i] + op.eps;
         };

         // computing min and max z coordinates over the arc
         A = U[5]*cdist*sTheta;  B = U[8]*cdist*sTheta;
         if (A != 0.0)
         {
            lomega1 = A*cos(lomega0) + B*sin(lomega0);
            uomega1 = A*cos(uomega0) + B*sin(uomega0);
            alpha = atan2(B,A);
            if (alpha < 0.0)
               opt = alpha + S.pi;
            else
               opt = alpha - S.pi;
            alpha = A*cos(alpha) + B*sin(alpha);
            alpha = maximum(alpha,lomega1,uomega1);
            opt = A*cos(opt) + B*sin(opt);
            opt = minimum(opt,lomega1,uomega1);
            S.lX[2][i] = S.lX[2][r1->otherId] - U[2]*cdist*cTheta + opt - op.eps;
            S.uX[2][i] = S.uX[2][r1->otherId] - U[2]*cdist*cTheta + alpha + op.eps;
         }
         else
         {
            S.lX[2][i] = X[2][i] - op.eps;  S.uX[2][i] = X[2][i] + op.eps;
         };
      };

      // performing the DDF pruning device
      pru = DDF(i,v,X,op.eps);

      // if it is necessary to refine the current solution
      if (pru > 0)
      {
         // verification of distance between the boxes (if DDF gave a negative result)
         if (BoxDDF(i,v,S.lX,S.uX,op.eps) == 0)
         {
            // if the distance between the boxes is feasible, then we can try to improve the current solution by local optimization
            k = 0;
            while (pru > 0 && k < 5)  // max 5 bound expansions allowed per BP call (introduced in version 0.3.1)
            {
               expandBounds(i+1,S.e,X,S.lX,S.uX,S.be);
               spg(i+1,v,X,S,&it,&obj);
               info->nspg++;
               pru = DDF(i,v,X,op.eps);
               k++;
            };
            if (pru == 0)  info->nspgok++;
         };
      };
      if (pru > 0)  info->pruning++;

      // if the current partial solution is OK (either since the beginning, or after the SPG call)
      if (pru == 0)
      {
         // verifying whether the partial solution is too close to the previous computed one
         if (check)
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

            // we will start to compare new solutions with previous ones
            check = true;

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
      if (i == 3)  if (op.symmetry > 0)  break;

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

// branch-and-prune (specific version for instances consisting of exact distances only)
void bp_exact(int i,int n,int m,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info)
{
   int h,k;
   int ldigits;
   double U[9];
   double tmp;
   double cdist,lde,mde;
   double cTheta,sTheta;
   double cosOmega,sinOmega[2];
   REFERENCE *r1,*r2,*r3;

   // signal handler
   signal(SIGINT,intHandler);

   // first call to BP (exact) ?
   if (i == 0)
   {
      // initializing the BP call counter
      info->ncalls = 0;

      // The first three vertices can be positioned by using the initial clique

      // vertex 0
      X[0][0] = 0.0;  X[1][0] = 0.0;  X[2][0] = 0.0;

      // vertex 1
      r1 = getReference(v,0,1);
      X[0][1] = -lowerBound(r1);  X[1][1] = 0.0;  X[2][1] = 0.0;

      // vertex 2
      r2 = getReference(v,1,2);
      cTheta = costheta(0,1,2,v,X);  sTheta = sqrt(1.0 - cTheta*cTheta);
      X[0][2] = -lowerBound(r1) + lowerBound(r2)*cTheta;  X[1][2] = lowerBound(r2)*sTheta;  X[2][2] = 0.0;

      i = i + 3;  // branching starts at vertex i+3
   };

   // updating BP call counter
   info->ncalls++;

   // reference vertices
   r3 = S.refs[i].r3;  r2 = S.refs[i].r2;  r1 = S.refs[i].r1;
   cdist = lowerBound(r1);

   // theta angles ("bond" angles)
   cTheta = costheta(r2->otherId,r1->otherId,i,v,X);
   sTheta = sqrt(1.0 - cTheta*cTheta);

   // generating U matrix (only once)
   UMatrix(r3->otherId,r2->otherId,r1->otherId,i,X,U);

   // omega angles (torsion angles)
   cosOmega = cosomega(r3->otherId,r2->otherId,r1->otherId,i,v,X,0.0);
   sinOmega[0] = sqrt(1.0 - cosOmega*cosOmega);
   sinOmega[1] = -sinOmega[0];
   if (op.symmetry == 2)
   {
      tmp = sinOmega[0];
      sinOmega[0] = sinOmega[1];
      sinOmega[1] = tmp;
   };

   // branching
   for (h = 0; h < 2 && keep_going; h++)
   {
      // monitor
      if (op.monitor)
      {
         ldigits = numberOfDigits(i);
         for (k = 0; k < info->ndigits; k++)  fprintf(stderr,"\b");
         for (k = 0; k < info->ndigits - ldigits; k++)  fprintf(stderr," ");
         fprintf(stderr,"%d",i);
      };

      // generating the coordinates for the current vertex (left-handed branch)
      genCoordinates(r1->otherId,i,X,U,cdist,cTheta,sTheta,cosOmega,sinOmega[h]);

      // performing the DDF pruning device
      if (DDF(i,v,X,op.eps) == 0)
      {
         // all distances are satisfied at the current layer
         if (i < n - 1)
         {
            // next vertex
            bp_exact(i+1,n,m,v,X,S,op,info);
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
         };
      }
      else
      {
         info->pruning++;
      };

      // if the sine of omega is zero, there is actually no branching
      if (fabs(sinOmega[h]) < 1.e-6)  break;

      // if only one solution is requested, bp stops as soon as the first solution is found
      if (op.allone == 1)  if (info->nsols > 0)  break;

      // the search stops after maxsols solutions
      if (info->nsols >= info->maxsols)  break;

      // skipping one half of the tree (optional)
      if (i == 3)  if (op.symmetry > 0)  break;
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

   return;
};

