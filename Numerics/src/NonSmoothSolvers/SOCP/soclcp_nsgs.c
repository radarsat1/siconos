/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "soclcp_projection.h"
/* #include "soclcp_Newton.h" */
/* #include "soclcp_Path.h" */
/* #include "soclcp_NCPGlockerFixedPoint.h" */
/* #include "soclcp_unitary_enumerative.h" */
#include "soclcp_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES

#include "debug.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void soclcp_nsgs_update(int cone, SecondOrderConeLinearComplementarityProblem* problem, SecondOrderConeLinearComplementarityProblem* localproblem, double * r, SolverOptions* options)
{
  /* Build a local problem for a specific cone
     r corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific cone
   r corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  soclcp_nsgs_fillMLocal(problem, localproblem, cone);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.rBlock,
     excluding the block corresponding to the current cone. ****/
  soclcp_nsgs_computeqLocal(problem, localproblem, r, cone);

  /* coefficient for current block*/
  localproblem->mu[0] = problem->mu[cone];

  /* index for current block*/
  localproblem->coneIndex[0] = 0;


  /* coefficient for current block*/
  localproblem->n = problem->coneIndex[cone+1] - problem->coneIndex[cone];


}
void soclcp_initializeLocalSolver_nsgs(Solver_soclcp_Ptr* solve, Update_soclcp_Ptr* update, FreeSolverNSGS_soclcp_Ptr* freeSolver, ComputeError_soclcp_Ptr* computeError, SecondOrderConeLinearComplementarityProblem* problem, SecondOrderConeLinearComplementarityProblem* localproblem, SolverOptions * localsolver_options)
{


  /** Connect to local solver */
  switch(localsolver_options->solverId)
  {
  /* Projection */
  case SICONOS_SOCLCP_ProjectionOnCone:
  {
    *solve = &soclcp_projectionOnCone_solve;
    *update = &soclcp_nsgs_update;
    *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_free;
    *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error;
    soclcp_projection_initialize(problem, localproblem);
    break;
  }
  /* case SICONOS_SOCLCP_ProjectionOnConeWithDiagonalization: */
  /* { */
  /*   *solve = &soclcp_projectionWithDiagonalization_solve; */
  /*   *update = &soclcp_projectionWithDiagonalization_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   soclcp_projection_initialize(problem, localproblem); */
  /*   break; */
  /* } */
  case SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration:
  {
    *solve = &soclcp_projectionOnConeWithLocalIteration_solve;
    *update = &soclcp_nsgs_update;
    *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projectionOnConeWithLocalIteration_free;
    *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error;
    soclcp_projectionOnConeWithLocalIteration_initialize(problem, localproblem,localsolver_options );
    break;
  }
  /* case SICONOS_SOCLCP_projectionOnConeWithRegularization: */
  /* { */
  /*   *solve = &soclcp_projectionOnCone_solve; */
  /*   *update = &soclcp_projection_update_with_regularization; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_with_regularization_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   soclcp_projection_initialize_with_regularization(problem, localproblem); */
  /*   break; */
  /* } */
  /* /\* Newton solver (Alart-Curnier) *\/ */
  /* case SICONOS_SOCLCP_AlartCurnierNewton: */
  /* { */
  /*   *solve = &soclcp_Newton_solve; */
  /*   *update = &soclcp_AC_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Newton_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   soclcp_Newton_initialize(problem, localproblem, localsolver_options); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_DampedAlartCurnierNewton: */
  /* { */
  /*   *solve = &soclcp_Newton_solve; */
  /*   *update = &soclcp_AC_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Newton_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   soclcp_Newton_initialize(problem, localproblem, localsolver_options); */
  /*   break; */
  /* } */
  /* /\* Newton solver (Glocker-Fischer-Burmeister)*\/ */
  /* case SICONOS_SOCLCP_NCPGlockerFBNewton: */
  /* { */
  /*   *solve = &soclcp_Newton_solve; */
  /*   *update = &NCPGlocker_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Newton_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   // *computeError = &fake_compute_error; */
  /*   soclcp_Newton_initialize(problem, localproblem, localsolver_options); */
  /*   break; */
  /* } */
  /* /\* Path solver (Glocker Formulation) *\/ */
  /* case SICONOS_SOCLCP_NCPGlockerFBPATH: */
  /* { */
  /*   *solve = &soclcp_Path_solve; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Path_free; */
  /*   *update = &NCPGlocker_update; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   // *computeError = &fake_compute_error; */
  /*   soclcp_Path_initialize(problem, localproblem, localsolver_options); */
  /*   break; */
  /* } */

  /* /\* Fixed Point solver (Glocker Formulation) *\/ */
  /* case SICONOS_SOCLCP_NCPGlockerFBFixedPoint: */
  /* { */
  /*   *solve = &soclcp_FixedP_solve; */
  /*   *update = &NCPGlocker_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_FixedP_free; */
  /*   *computeError = &fake_compute_error_nsgs; */
  /*   soclcp_FixedP_initialize(problem, localproblem, localsolver_options); */
  /*   break; */
  /* } */
  /* /\*iparam[4] > 10 are reserved for Tresca resolution *\/ */
  /* case SICONOS_SOCLCP_projectionOnCylinder: */
  /* { */
  /*   *solve = &soclcp_projectionOnCylinder_solve; */
  /*   *update = &soclcp_projectionOnCylinder_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_Tresca_compute_error; */
  /*   soclcp_projection_initialize(problem, localproblem); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_QUARTIC: */
  /* { */
  /*   *solve = &soclcp_unitary_enumerative_solve; */
  /*   *update = &soclcp_nsgs_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_unitary_enumerative_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   soclcp_unitary_enumerative_initialize(localproblem); */
  /*   break; */
  /* } */
  /* case SICONOS_SOCLCP_QUARTIC_NU: */
  /* { */
  /*   *solve = &soclcp_unitary_enumerative_solve; */
  /*   *update = &soclcp_nsgs_update; */
  /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_unitary_enumerative_free; */
  /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
  /*   soclcp_unitary_enumerative_initialize(localproblem); */
  /*   break; */
  /* } */
  default:
  {
    fprintf(stderr, "Numerics, soclcp_nsgs failed. Unknown internal solver : %s.\n", idToName(localsolver_options->solverId));
    exit(EXIT_FAILURE);
  }
  }
}
void soclcp_nsgs_computeqLocal(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, double *r, int cone)
{

  double *qLocal = localproblem->q;

  int n = problem->n;

  int normal = problem->coneIndex[cone];

  int dim = problem->coneIndex[cone+1]-problem->coneIndex[cone];


  /* r current block set to zero, to exclude current cone block */
  int i;
  double * rsave = (double *) malloc(dim*sizeof(double));
  for(i = 0; i < dim; i++)
  {
    rsave[i] = r[normal + i];
    r[normal + i] = 0.0;
  }

  /* qLocal computation*/

  for(i = 0; i < dim; i++)
  {
    qLocal[i] = problem->q[normal +i];
  }
  if(problem->M->storageType == 0)
  {
    double * MM = problem->M->matrix0;
    int incx = n, incy = 1;

    for(i = 0; i < dim; i++)
    {
      qLocal[i] += cblas_ddot(n , &MM[normal+i] , incx , r , incy);
    }
  }
  else if(problem->M->storageType == 1)
  {
    /* qLocal += rowMB * r
    with rowMB the row of blocks of MGlobal which corresponds to the current cone
    */
    DEBUG_PRINTF("dim= %i\n", dim);
    rowProdNoDiagSBM(n, dim, cone, problem->M->matrix1, r, qLocal, 0);
  }
  for(int i = 0; i < dim; i++)
  {
    r[normal + i]=  rsave[i];
  }
  free(rsave);
}

void soclcp_nsgs_fillMLocal(SecondOrderConeLinearComplementarityProblem * problem, SecondOrderConeLinearComplementarityProblem * localproblem, int cone)
{

  NumericsMatrix * MGlobal = problem->M;

  int n = problem->n;

  int storageType = MGlobal->storageType;
  if(storageType == 0)
    // Dense storage
  {
    int normal = problem->coneIndex[cone];
    int dim = problem->coneIndex[cone+1]-problem->coneIndex[cone+1];
    int inc = n * normal;
    double * MM = MGlobal->matrix0;
    double * MLocal =  localproblem->M->matrix0;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    for(int j =0; j< dim; j++)
    {
      for(int i = 0; i < dim; i++)  MLocal[i+j*dim] = MM[inc + normal+i];
      inc += n;
    }
  }
  else if(storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, cone);
    localproblem->M->matrix0 = MGlobal->matrix1->block[diagPos];
  }
  else
    numericsError("soclcp_projection -", "unknown storage type for matrix M");
}


/* swap two indices */
void uint_swap(unsigned int *a, unsigned int *b);


/* shuffle an unsigned array */
void uint_shuffle(unsigned int *a, unsigned int n);


void soclcp_nsgs(SecondOrderConeLinearComplementarityProblem* problem, double *r, double *v, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of cones */
  unsigned int nc = problem->nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  if(*info == 0)
    return;

  if(options->numberOfInternalSolvers < 1)
  {
    numericsError("soclcp_nsgs", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }
  assert(&options[1]);

  SolverOptions * localsolver_options = options->internalSolvers;


  Solver_soclcp_Ptr local_solver = NULL;
  Update_soclcp_Ptr update_localproblem = NULL;
  FreeSolverNSGS_soclcp_Ptr freeSolver = NULL;
  ComputeError_soclcp_Ptr computeError = NULL;

  unsigned int cone ;
  int isConeDimensionsEqual=1;
  /* Connect local solver and local problem*/
  unsigned int dim = problem->coneIndex[1]-problem->coneIndex[0];
  for(cone = 1; cone < nc; cone++)
  {
    if(problem->coneIndex[cone+1]-problem->coneIndex[cone] != dim)
    {
      isConeDimensionsEqual=0;
      break;
    }
  }
  SecondOrderConeLinearComplementarityProblem* localproblem;


  if(isConeDimensionsEqual)
  {
    localproblem = (SecondOrderConeLinearComplementarityProblem*)malloc(sizeof(SecondOrderConeLinearComplementarityProblem));
    localproblem->nc = 1;
    localproblem->n = dim;
    localproblem->q = (double*)malloc(dim * sizeof(double));
    localproblem->mu = (double*)malloc(sizeof(double));
    localproblem->coneIndex = (unsigned int*)malloc(2*sizeof(unsigned int));
    localproblem->coneIndex[0]=0;
    localproblem->coneIndex[1]=dim;

    if(problem->M->storageType == 0)
    {
      localproblem->M = createNumericsMatrixFromData(NM_DENSE, dim, dim,
                        malloc(dim*dim* sizeof(double)));
    }
    else
    {
      localproblem->M = createNumericsMatrix(NM_DENSE, dim, dim);
    }
  }
  else
  {
    fprintf(stderr, "soclcp_nsgs error: not yet implemented.\n");
    exit(EXIT_FAILURE);
  }

  soclcp_initializeLocalSolver_nsgs(&local_solver, &update_localproblem,
                                    (FreeSolverNSGS_soclcp_Ptr *)&freeSolver, &computeError,
                                    problem , localproblem, localsolver_options);

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  unsigned int *scones = NULL;

  if(iparam[9])  /* shuffle */
  {
    scones = (unsigned int *) malloc(nc * sizeof(unsigned int));
    for(unsigned int i = 0; i<nc ; ++i)
    {
      scones[i] = i;
    }
    uint_shuffle(scones, nc);
  }




  /*  dparam[0]= dparam[2]; // set the tolerance for the local solver */

  if(iparam[1] == 1 || iparam[1] == 2)
  {
    int n =problem->n;
    double * rold = (double*)malloc(n*sizeof(double)); // save memory if isConeDimensionsEqual
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /* Loop through the cone  */
      //cblas_dcopy( n , q , incx , v , incy );
      error = 0.0;
      for (int i =0; i < n; i++) rold[i] = r[i];
      for(unsigned int i= 0 ; i < nc ; ++i)
      {
        if(iparam[9])
        {
          cone = scones[i];
        }
        else
        {
          cone = i;
        }

        if(verbose > 1) printf("----------------------------------- Cone Number %i\n", cone);
        (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
        
        localsolver_options->iparam[4] = cone;

        (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]) , localsolver_options);
      }
      for (int i=0;i < n; i++)
      {
        error += pow(r[i] - rold[i], 2);
      }

      /* **** Criterium convergence **** */
      error = sqrt(error);
      if(verbose > 0)
        printf("----------------------------------- SOCLCP - NSGS - Iteration %i Error = %14.7e\n", iter, error);
      if(error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }

    if(iparam[1] == 1)  /* Full criterium */
    {
      double absolute_error;
      (*computeError)(problem, r , v, tolerance, options, &absolute_error);
      if(verbose > 0)
      {
        if(absolute_error > error)
        {
          printf("----------------------------------- SOCLCP - NSGS - Warning absolute Error = %14.7e is larger than incremental error = %14.7e\n", absolute_error, error);
        }
      }
    }
    free(rold);
  }
  else
  {

    if(iparam[9])
    {
      int withRelaxation=iparam[8];
      if(withRelaxation)
      {
        int n = problem->n;
        double * rold = (double*)malloc(n*sizeof(double)); // save memory if isConeDimensionsEqual
        double omega = dparam[8];
        unsigned int dim;
        while((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          for (int i =0; i < n; i++) rold[i] = r[i];
          /* Loop through the cone points */
          //cblas_dcopy( n , q , incx , v , incy );
          for(unsigned int i= 0 ; i < nc ; ++i)
          {
            cone = scones[i];

            if(verbose > 1) printf("----------------------------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[4] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);

            dim = problem->coneIndex[cone+1]-problem->coneIndex[cone];
            for (unsigned int i =0; i <dim; ++i)
            {
              r[problem->coneIndex[cone]+i] = omega*r[problem->coneIndex[cone]+i]
                +(1.0-omega)*rold[problem->coneIndex[cone]+i];
            }
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r , v, tolerance, options, &error);

          if(verbose > 0)
            printf("----------------------------------- SOCLCP - NSGS - Iteration %i Error = %14.7e\n", iter, error);

          if(error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if(options->callback)
          {
            options->callback->collectStatsIteration(options->callback->env, problem->n,
                                                     r, v,
                                                     error, NULL);
          }
        }
        free(rold);
      }
      else
      {
        while((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the cone  */
          //cblas_dcopy( n , q , incx , v , incy );
          for(unsigned int i= 0 ; i < nc ; ++i)
          {
            
            cone = scones[i];


            if(verbose > 1) printf("----------------------------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[4] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);

          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r , v, tolerance, options, &error);

          if(verbose > 0)
            printf("----------------------------------- SOCLP - NSGS - Iteration %i Error = %14.7e\n", iter, error);

          if(error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if(options->callback)
          {
            options->callback->collectStatsIteration(options->callback->env, problem->n,
                                                     r, v,
                                                     error, NULL);
          }
        }

      }


    }
    else
    {
      int withRelaxation=iparam[8];
      if(withRelaxation)
      {
        int n = problem->n;
        double * rold = (double*)malloc(n*sizeof(double)); // save memory if isConeDimensionsEqual
        double omega = dparam[8];
        while((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          for (int i =0; i < n; i++) rold[i] = r[i];

          /* Loop through the cones  */
          //cblas_dcopy( n , q , incx , v , incy );
          for(cone= 0 ; cone< nc ; ++cone)
          {

            if(verbose > 1) printf("----------------------------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[4] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);

            dim = problem->coneIndex[cone+1]-problem->coneIndex[cone];
            for (unsigned int i =0; i <dim; ++i)
            {
              r[problem->coneIndex[cone]+i] = omega*r[problem->coneIndex[cone]+i]
                +(1.0-omega)*rold[problem->coneIndex[cone]+i];
            } 
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r , v, tolerance, options, &error);

          if(verbose > 0)
            printf("----------------------------------- SOCLCP - NSGS - Iteration %i Error = %14.7e\n", iter, error);

          if(error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if(options->callback)
          {
            options->callback->collectStatsIteration(options->callback->env, problem->n,
                                                     r, v,
                                                     error, NULL);
          }
        }

      }
      else
      {
        while((iter < itermax) && (hasNotConverged > 0))
        {
          ++iter;
          /* Loop through the cones  */
          //cblas_dcopy( n , q , incx , v , incy );
          for(cone= 0 ; cone < nc ; ++cone)
          {

            if(verbose > 1) printf("----------------------------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[4] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);

          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r , v, tolerance, options, &error);

          if(verbose > 0)
            printf("----------------------------------- SOCLCP - NSGS - Iteration %i Error = %14.7e\n", iter, error);

          if(error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if(options->callback)
          {
            options->callback->collectStatsIteration(options->callback->env, problem->n,
                r, v,
                error, NULL);
          }
        }

      }
    }





  }
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /***** Free memory *****/
  (*freeSolver)(problem,localproblem,localsolver_options);
  if(problem->M->storageType == 0 && localproblem->M->matrix0 != NULL)
  {
    free(localproblem->M->matrix0);
  }
  localproblem->M->matrix0 = NULL;
  freeSecondOrderConeLinearComplementarityProblem(localproblem);

  if(scones)  /* shuffle */
  {
    free(scones);
  }

}

int soclcp_nsgs_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if(verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  /*  strcpy(options->solverName,"NSGS");*/
  options->solverId = SICONOS_SOCLCP_NSGS;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  options->callback = NULL;
  options->numericsOptions = NULL;
  for(i = 0; i < 10; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  soclcp_projection_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}