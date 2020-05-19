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
              May 19 2020  v.0.3.2  new bound expansion technique in bp version for interval distances
                                    bp_exact may choose the "best" triplet of discretization vertices
                                    a time limit for both bp implementations can now be set up
*********************************************************************************************************/

#include "bp.h"

extern double INFTY;
bool keep_going = true;
bool newsol = false;
bool backtracking = false;
bool check = false;
bool PRINTED = false;
struct timeval startime;

// signal catcher
void intHandler(int a)
{
   fprintf(stderr," signal caught: stopping (partial solution printed if -P or -p options used)");
   keep_going = false;
};

// branch-and-prune (general version)
// -> i, current vertex of be realized
// -> n, total number of vertices forming the instance
// -> v, array of vertices
// -> X, current matrix of coordinates
// -> S, SEARCH structure
// -> op, OPTION structure
// -> info, INFORMATION structure
void bp(int i,int n,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info)
{
   int j,k;
   int it,nb;
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
   double pperr,perr;
   Omega *current;
   omegaList omegaL;
   REFERENCE *r1,*r2,*r3;
   struct timeval currentime;

   // signal handler
   signal(SIGINT,intHandler);

   // first call to BP?
   if (i == 0)
   {
      // initializing the BP call counter
      info->ncalls = 0;

      // vertex 0
      X[0][0] =  0.0;  X[1][0] = 0.0;  X[2][0] = 0.0;
      createBox(0,X,op.eps,S.lX,S.uX);

      // vertex 1
      r1 = getReference(v,0,1);
      X[0][1] = -lowerBound(r1);  X[1][1] = 0.0;  X[2][1] = 0.0;
      createBox(1,X,op.eps,S.lX,S.uX);

      // vertex 2
      r2 = getReference(v,1,2);
      cTheta = costheta(0,1,2,v,X);  sTheta = sqrt(1.0 - cTheta*cTheta);
      X[0][2] = -lowerBound(r1) + lowerBound(r2)*cTheta;  X[1][2] = lowerBound(r2)*sTheta;  X[2][2] = 0.0; 
      createBox(2,X,op.eps,S.lX,S.uX);

      // we start to count the time for BP from this point
      gettimeofday(&startime,0);

      // branching starts at vertex i+3
      i = i + 3;
   };

   // updating BP call counter
   info->ncalls++;
   it = 0;

   // reference vertices
   r3 = S.refs[i].r3;  r2 = S.refs[i].r2;  r1 = S.refs[i].r1;
   cdist = lowerBound(r1);
   
   // theta angle ("bond" angles)
   cTheta = costheta(otherVertexId(r2),otherVertexId(r1),i,v,X);
   sTheta = sqrt(1.0 - cTheta*cTheta);

   // generating U matrix (only once)
   UMatrix(otherVertexId(r3),otherVertexId(r2),otherVertexId(r1),i,X,U);

   // omega angle (torsion angles)
   nb = 2;
   cosOmega00 = cosomega(otherVertexId(r3),otherVertexId(r2),otherVertexId(r1),i,v,X,0.0,op.eps);
   cosOmega01 = cosomega(otherVertexId(r3),otherVertexId(r2),otherVertexId(r1),i,v,X,1.0,op.eps);
   if (cosOmega00 == -2.0 || cosOmega01 == -2.0)  return;  // infeasibility already detected
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
      genCoordinates(otherVertexId(r1),i,X,U,cdist,cTheta,sTheta,cos(omega),sin(omega));

      // generation of the box inscribing the arc
      if (isExactDistance(r3,op.eps))
      {
         // the box has the size equal to the tolerance over the three dimensions
         createBox(i,X,op.eps,S.lX,S.uX);
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
            S.lX[0][i] = S.lX[0][otherVertexId(r1)] - U[0]*cdist*cTheta + opt - op.eps;
            S.uX[0][i] = S.uX[0][otherVertexId(r1)] - U[0]*cdist*cTheta + alpha + op.eps;
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
            S.lX[1][i] = S.lX[1][otherVertexId(r1)] - U[1]*cdist*cTheta + opt - op.eps;
            S.uX[1][i] = S.uX[1][otherVertexId(r1)] - U[1]*cdist*cTheta + alpha + op.eps;
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
            S.lX[2][i] = S.lX[2][otherVertexId(r1)] - U[2]*cdist*cTheta + opt - op.eps;
            S.uX[2][i] = S.uX[2][otherVertexId(r1)] - U[2]*cdist*cTheta + alpha + op.eps;
         }
         else
         {
            S.lX[2][i] = X[2][i] - op.eps;  S.uX[2][i] = X[2][i] + op.eps;
         };
      };

      // expanding the box until some reference distances are not satisfied
      expandBounds(i,v,S.lX,S.uX,op.be,op.eps);

      // performing the DDF pruning device
      perr = DDF(i,v,X);

      // if it is necessary to refine the current solution
      if (perr > op.eps)
      {
         // verification of distance between the boxes (if DDF gave a negative result)
         if (BoxDDF(i,v,S.lX,S.uX) < op.eps)
         {
            k = 0;
            do // if the distance between the boxes is feasible, 
            {     // then we can try to improve the current solution by local optimization
               pperr = perr;
               spg(i+1,v,X,S,op,info,&it,&obj);
               info->nspg++;
               perr = DDF(i,v,X);
               reCenterBounds(i+1,v,X,S.lX,S.uX,op.be,op.eps);
               if (perr < op.eps)  info->nspgok++;
               k++;
            }
            while (perr > op.eps && pperr - perr > op.eps && k < 20 && keep_going);
         };
      };
      if (perr > op.eps)  info->pruning++;

      // if the current partial solution is OK (either since the beginning, or after local optimization)
      if (perr < op.eps)
      {
         // verifying whether the partial solution is too close to the previous one
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
            bp(i+1,n,v,X,S,op,info);
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

      // maxtime limit reached?
      gettimeofday(&currentime,0);
      if (currentime.tv_sec - startime.tv_sec > op.maxtime)  keep_going = false;

      // skipping one half of the tree (optional)
      if (i == 3)  if (op.symmetry > 0)
      {
         info->pruning++;
         break;
      };

      // if only one solution is requested, bp stops as soon as the first solution is found
      if (op.allone == 1)  if (info->nsols > 0)  break;

      // the search stops after maxsols solutions
      if (info->nsols >= info->maxsols)  break;

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
               printfile(i,v,X,info->output,0);
            else
               printpdb(i,v,X,info->output,0);
            PRINTED = true;
         };
      };
   };

   // freeing memory space for omega list
   freeOmegaList(omegaL);

   return;
};

// branch-and-prune (specific version for instances consisting of exact, and precise, distances)
// -> i, current vertex of be realized
// -> n, total number of vertices forming the instance
// -> v, array of vertices
// -> X, current matrix of coordinates
// -> S, SEARCH structure
// -> op, OPTION structure
// -> info, INFORMATION structure
void bp_exact(int i,int n,VERTEX *v,double **X,SEARCH S,OPTION op,INFORMATION *info)
{
   int h,k;
   int ldigits;
   double U[9];
   double tmp;
   double cdist,lde,mde;
   double cTheta,sTheta;
   double cosOmega,sinOmega[2];
   double perr,berr;
   triplet t,best;
   struct timeval currentime;

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
      t.r1 = getReference(v,0,1);
      X[0][1] = -lowerBound(t.r1);  X[1][1] = 0.0;  X[2][1] = 0.0;

      // vertex 2
      t.r2 = getReference(v,1,2);
      cTheta = costheta(0,1,2,v,X);  sTheta = sqrt(1.0 - cTheta*cTheta);
      X[0][2] = -lowerBound(t.r1) + lowerBound(t.r2)*cTheta;  X[1][2] = lowerBound(t.r2)*sTheta;  X[2][2] = 0.0;

      // we start to count the time for BP from this point
      gettimeofday(&startime,0);

      // branching starts at vertex i+3
      i = i + 3;
   };

   // updating BP call counter
   info->ncalls++;

   // if we're not backtracking, we are exploring the tree for a new solution
   if (!backtracking)  newsol = false;

   // selection of the discretization vertices
   cTheta = 0.0;
   if (info->consec)
   {
      // the consecutivity assumption is satisfied
      best.r1 = getReference(v,i,i-1);
      best.r2 = getReference(v,i,i-2);
      best.r3 = getReference(v,i,i-3);

      // is this triplet too flat?
      cTheta = costheta(otherVertexId(best.r2),otherVertexId(best.r1),i,v,X);
   };

   // we select the discretization vertices leading to the smallest error
   // (if the one above cannot be defined or it is too flat)
   if (fabs(cTheta) < op.eps)
   {
      berr = INFTY;
      t = nullTriplet();  best = nullTriplet();

      do // trying out all possible triplets
      {
         t = nextTripletRef(v[i].ref,t,op.eps);
         if (isValidTriplet(t,op.eps))
         {
            // theta angle ("bond" angles)
            cTheta = costheta(otherVertexId(t.r2),otherVertexId(t.r1),i,v,X);
            if (fabs(cTheta) < op.eps)  continue;  // before invoking bp, it was verified
            sTheta = sqrt(1.0 - cTheta*cTheta);    // that "flattest" triplet is not too flat!
            if (sTheta < op.eps)  continue;
            cdist = lowerBound(t.r1);

            // generating U matrix
            UMatrix(otherVertexId(t.r3),otherVertexId(t.r2),otherVertexId(t.r1),i,X,U);

            // omega angle (torsion angles)
            cosOmega = cosomega(otherVertexId(t.r3),otherVertexId(t.r2),otherVertexId(t.r1),i,v,X,0.0,op.eps);
            if (cosOmega == -2.0)  continue;
            sinOmega[0] = sqrt(1.0 - cosOmega*cosOmega);

            // generating the coordinates for the vertex by using the current triplet t
            genCoordinates(otherVertexId(t.r1),i,X,U,cdist,cTheta,sTheta,cosOmega,sinOmega[0]);

            // verifying the error over the entire set of reference distances
            perr = DDF(i,v,X);
            if (perr < berr)
            {
               best.r1 = t.r1;
               best.r2 = t.r2;
               best.r3 = t.r3;
               berr = perr;
            };
         };
      }
      while (!isNullTriplet(t) && keep_going);
   };

   // using the best found triplet to compute the coordinates
   if (isValidTriplet(best,op.eps))
   {
      // theta angle of best
      cTheta = costheta(otherVertexId(best.r2),otherVertexId(best.r1),i,v,X);
      sTheta = sqrt(1.0 - cTheta*cTheta);
      cdist = lowerBound(best.r1);

      // generating U matrix
      UMatrix(otherVertexId(best.r3),otherVertexId(best.r2),otherVertexId(best.r1),i,X,U);

      // omega angle of best
      cosOmega = cosomega(otherVertexId(best.r3),otherVertexId(best.r2),otherVertexId(best.r1),i,v,X,0.0,op.eps);
      if (cosOmega == -2.0)  return;  // infeasibility already detected
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
         if (op.monitor && (i == 4 || i%10 == 0 || i == n - 1))
         {
            ldigits = numberOfDigits(i);
            for (k = 0; k < info->ndigits; k++)  fprintf(stderr,"\b");
            for (k = 0; k < info->ndigits - ldigits; k++)  fprintf(stderr," ");
            fprintf(stderr,"%d",i);
         };

         // generating the coordinates for the vertex by using the best triplet
         genCoordinates(otherVertexId(best.r1),i,X,U,cdist,cTheta,sTheta,cosOmega,sinOmega[h]);

         // performing the DDF pruning device
         if (DDF(i,v,X) < op.eps)
         {
            // all distances are satisfied at the current layer
            if (i < n - 1)
            {
               // next vertex
               backtracking = false;
               bp_exact(i+1,n,v,X,S,op,info);
               backtracking = true;
            }
            else
            {
               // solution found
               newsol = true;
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

         // maxtime limit reached?
         gettimeofday(&currentime,0);
         if (currentime.tv_sec - startime.tv_sec > op.maxtime)  keep_going = false;

         // skipping one half of the tree (optional)
         if (i == 3)  if (op.symmetry > 0)
         {
            info->pruning++;
            break;
         };

         // if we are backtracking after a solution was found, we don't explore the second half of the branch if:
         // - the consecutivity assumption is satisfied and the current vertex is not symmetric
         // - the sine of the omega angle is too small (fixed tolerance)
         if (i > 3)  if (newsol)  if (sinOmega[0] < 0.05 || (info->consec && !S.sym[i]))
         {
            info->pruning++;
            break;
         };

         // if only one solution is requested, bp stops as soon as the first solution is found
         if (op.allone == 1)  if (info->nsols > 0)  break;

         // the search stops after maxsols solutions
         if (info->nsols >= info->maxsols)  break;
      };
   };

   // handling ^C signal catcher
   if (!keep_going)
   {
      if (!PRINTED)
      {
         if(op.print > 0 && info->nsols == 0)
         {
            if (op.format == 0)
               printfile(i,v,X,info->output,0);
            else
               printpdb(i,v,X,info->output,0);
            PRINTED = true;
         };
      };
   };

   return;
};

