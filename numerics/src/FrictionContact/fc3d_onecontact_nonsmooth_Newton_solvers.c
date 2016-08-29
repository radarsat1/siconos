/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "NonSmoothNewton.h"
#include "NCP_Solvers.h"
#include "fc3d_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include "FischerBurmeister.h"
#include "SiconosBlas.h"
#include "AlartCurnierGenerated.h"
#include "fc3d_GlockerFischerBurmeister_functions.h"
#include "op3x3.h"


#define DEBUG_CHECK


#define OPTI_RHO
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

static computeNonsmoothFunction  Function = NULL;
static NewtonFunctionPtr F = NULL;
static NewtonFunctionPtr jacobianF = NULL;
static UpdateSolverPtr updateSolver = NULL;
static PostSolverPtr postSolver = NULL;
static FreeSolverNSGSPtr freeSolver = NULL;

/* size of a block */
static int Fsize;
void fc3d_AC_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * options);
void fc3d_AC_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * options)
{
  /*
    In initialize, these operators are "connected" to their corresponding static variables, that will be used to build local problem
    for each considered contact.
    Local problem is built during call to update (which depends on the storage type for M).
  */

  /* localFC3D = localproblem; */
  /* globalFC3D = problem; */
  DEBUG_PRINTF("fc3d_AC_initialize starts with options->iparam[10] = %i\n", options->iparam[10]);
  if (options->iparam[10] == 0 )
  {
    Function = &(computeAlartCurnierSTD);
  }
  else if (options->iparam[10] == 1 )
  {
    Function = &(computeAlartCurnierJeanMoreau);
  }
  else if (options->iparam[10] == 2 )
  {
    Function = &(fc3d_AlartCurnierFunctionGenerated);
  }
  else if (options->iparam[10] == 3 )
  {
    Function = &fc3d_AlartCurnierJeanMoreauFunctionGenerated;;
  }
  else if (options->iparam[10] == 4 )
  {
    Function = NULL;
  }


/* #ifdef AC_STD */
/* computeNonsmoothFunction  Function = &(computeAlartCurnierSTD); */
/* #endif */
/* #ifdef AC_JeanMoreau */
/* computeNonsmoothFunction  Function = &(computeAlartCurnierJeanMoreau); */
/* #endif */

/* // computeAlartCurnier[JeanMoreau] == AC_Generated */

/* #ifdef AC_Generated */
/* computeNonsmoothFunction  Function = &(fc3d_AlartCurnierFunctionGenerated); */
/* #endif */

/* // HandMade not done */
/* #ifdef AC_HandMade */
/* computeNonsmoothFunction  Function = &(fc3d_AlartCurnierFunctionHandMade); */
/* #endif */

}

void fc3d_AC_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);

void fc3d_AC_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options)
{
}
void fc3d_AC_free_P(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);
void fc3d_AC_free_P(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options)
{
  free(localsolver_options->dWork);
  localsolver_options->dWork=NULL;
}
void fc3d_AC_post(int contact, double* reaction);
void fc3d_AC_post(int contact, double* reaction)
{
  /* This function is required in the interface but useless in Alart-Curnier case */
}

void fc3d_onecontact_nonsmooth_Newton_solvers_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem,    SolverOptions * localsolver_options)
{

  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */

  /* Alart-Curnier formulation */
  if (localsolver_options->solverId == SICONOS_FRICTION_3D_ONECONTACT_NSN_AC)
  {
    fc3d_AC_initialize(problem, localproblem,localsolver_options);
    /* Fsize = 3; */
    /* F = &F_AC; */
    /* jacobianF = &jacobianF_AC; */
    /*     updateSolver = &fc3d_AC_update; */
    postSolver = &fc3d_AC_post;
    freeSolver = &fc3d_AC_free;
  }
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP)
  {
    fc3d_AC_initialize(problem, localproblem,localsolver_options);
    /* Fsize = 3; */
    /* F = &F_AC; */
    /* jacobianF = &jacobianF_AC; */
    /*     updateSolver = &fc3d_AC_update; */
    postSolver = &fc3d_AC_post;
    freeSolver = &fc3d_AC_free;

  }
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_P)
  {
    fc3d_AC_initialize(problem, localproblem,localsolver_options);
    fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem, localsolver_options );
    /* Fsize = 3; */
    /* F = &F_AC; */
    /* jacobianF = &jacobianF_AC; */
    /*     updateSolver = &fc3d_AC_update; */
    postSolver = &fc3d_AC_post;
    freeSolver = &fc3d_AC_free_P;

  }




  /* Glocker formulation - Fischer-Burmeister function used in Newton */
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_NCPGlockerFBNewton)
  {
    Fsize = 5;
    NCPGlocker_initialize(problem, localproblem);
    F = &F_GlockerFischerBurmeister;
    jacobianF = &jacobianF_GlockerFischerBurmeister;
    /*     updateSolver = &NCPGlocker_update; */
    postSolver = &NCPGlocker_post;
    freeSolver = (FreeSolverNSGSPtr)&NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

int fc3d_onecontact_nonsmooth_Newton_solvers_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{

  /*  (*updateSolver)(contact, reaction); */
  double * reactionBlock = reaction;

  int * iparam = options->iparam;
  double * dparam = options->dparam;

  int info;
  if (options->solverId == SICONOS_FRICTION_3D_ONECONTACT_NSN_AC)
  {
    info = fc3d_onecontact_nonsmooth_Newton_solvers_solve_direct(localproblem, reactionBlock, iparam, dparam);
  }
  else if (options->solverId == SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP)
  {
    info = fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped(localproblem, reactionBlock, iparam, dparam);
  }
  else if (options->solverId == SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_P)
  {
    info = fc3d_projectionOnConeWithLocalIteration_solve (localproblem, reactionBlock , options);
    info = fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped(localproblem, reactionBlock, iparam, dparam);
  }
  else
  {
    info = nonSmoothDirectNewton(Fsize, reactionBlock, &F, &jacobianF, iparam, dparam);
  }
  if (info > 0)
  {
    if (verbose > 0)
    {
      printf("Numerics, fc3d_onecontact_nonsmooth_Newton_solvers_solve, warning. reached max. number of iterations (%i) without convergence for contact %i. Residual = %12.8e\n", options->iparam[0], options->iparam[4], dparam[1]);
      /* note : exit on failure should be done in DefaultCheckSolverOutput */
    }
  }
  return info;
  /*  (*postSolver)(contact,reaction); */
}

void fc3d_onecontact_nonsmooth_Newton_solvers_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options)
{
  F = NULL;
  jacobianF = NULL;
  updateSolver = NULL;
  postSolver = NULL;
  (*freeSolver)(problem, localproblem, localsolver_options);
}

void fc3d_onecontact_nonsmooth_Newton_solvers_computeError(int n, double* velocity, double*reaction, double * error)
{
  /*   int numberOfContacts = n/3; */
  /*   int sizeGlobal = numberOfContacts*FSize; */
  /*   //  double * FGlobal = (double*)malloc(sizeGlobal*sizeof(*FGlobal));  */
  /*   (*computeFGlobal)(reaction,velocity); */
  /*   int i; */
  /*   double Fz; */
  /*   *error = 0; */
  /*   for(i=0;i<sizeGlobal;++i) */
  /*     { */
  /*       Fz = velocity[i]*reaction[i]; */
  /*       if(Fz>0) */
  /*  *error+=Fz; */
  /*       if(reaction[i]<0) */
  /*  *error+=reaction[i]; */
  /*       if(velocity[i]<0) */
  /*  *error+=velocity[i]; */
  /*     } */

  /*   // (*computeVelocity)(FGlobal); */

  /*   free(FGlobal); */

}


static void AC_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
{

  NumericsMatrix * MGlobal = problem->M;
  int n = 3 * problem->numberOfContacts;

  // Dense storage
  int storageType = MGlobal->storageType;
  if (storageType == 0)
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    double * MM = MGlobal->matrix0;
    double * MLocal =  localproblem->M->matrix0;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    MLocal[0] = MM[inc + in];
    MLocal[1] = MM[inc + it];
    MLocal[2] = MM[inc + is];
    inc += n;
    MLocal[3] = MM[inc + in];
    MLocal[4] = MM[inc + it];
    MLocal[5] = MM[inc + is];
    inc += n;
    MLocal[6] = MM[inc + in];
    MLocal[7] = MM[inc + it];
    MLocal[8] = MM[inc + is];
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, contact);
    localproblem->M->matrix0 = MGlobal->matrix1->block[diagPos];
    /*     cblas_dcopy(9, MGlobal->matrix1->block[diagPos], 1, localproblem->M->matrix0 , 1); */

  }
  else if (storageType == 2)
  {
    /* ok, we maintain the sparseblock storage from the sparse one */
    if (!problem->M->matrix1)
    {
      problem->M->matrix1 = (SparseBlockStructuredMatrix*) malloc(sizeof(SparseBlockStructuredMatrix));
      problem->M->matrix1->block = NULL;
      problem->M->matrix1->index1_data = NULL;
      problem->M->matrix1->index2_data = NULL;
    }
    sparseToSBM(problem->dimension, NM_triplet(problem->M), problem->M->matrix1);
    int diagPos = getDiagonalBlockPos(problem->M->matrix1, contact);
    localproblem->M->matrix0 = problem->M->matrix1->block[diagPos];

    /* Direct access is below, but there is a bug, find it! */

    /* CSparseMatrix* MM = NM_triplet(MGlobal); */
    /* int in = 3 * contact, it = in + 1, is = it + 1; */
    /* int inc = n * in; */
    /* double * MLocal =  localproblem->M->matrix0; */
    /* /\* The part of MM which corresponds to the current block is copied into MLocal *\/ */
    /* MLocal[0] = MM->x[inc + in]; */
    /* MLocal[1] = MM->x[inc + it]; */
    /* MLocal[2] = MM->x[inc + is]; */
    /* inc += n; */
    /* MLocal[3] = MM->x[inc + in]; */
    /* MLocal[4] = MM->x[inc + it]; */
    /* MLocal[5] = MM->x[inc + is]; */
    /* inc += n; */
    /* MLocal[6] = MM->x[inc + in]; */
    /* MLocal[7] = MM->x[inc + it]; */
    /* MLocal[8] = MM->x[inc + is]; */
  }
  else
  {
    numericsError("fc3d_AlartCurnier:AC_fillMLocal() -", "unknown storage type for matrix M");
  }

}

void fc3d_onecontact_nonsmooth_Newton_AC_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double * reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  AC_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  fc3d_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}


int fc3d_onecontact_nonsmooth_Newtow_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 14;
  options->dSize = 14;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);

  options->iparam[0] = 10;

  options->dparam[0] = 1e-16;
  options->iparam[10] = 0;     /* 0 STD AlartCurnier,
                                  1 JeanMoreau,
                                  2 STD generated,
                                  3 JeanMoreau generated */
  options->iparam[11] = 0;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[12] = 10;   /* max iter line search */

  return 0;
}


int fc3d_onecontact_nonsmooth_Newton_solvers_solve_direct(FrictionContactProblem* localproblem, double * R, int *iparam, double *dparam)
{
  double mu = localproblem->mu[0];
  double * qLocal = localproblem->q;

  double * MLocal = localproblem->M->matrix0;

  double Tol = dparam[0];
  double itermax = iparam[0];


  int i, j, k, inew;

  // store the increment
  double dR[3] = {0., 0., 0.};

  // store the value fo the function
  double F[3] = {0., 0., 0.};

  // Store the (sub)-gradient of the function
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Compute values of Rho (should be here ?)
  double rho[3] = {1., 1., 1.};
#ifdef OPTI_RHO
  computerho(localproblem, rho);
  DEBUG_PRINTF("rho[0] = %4.2e, rho[1] = %4.2e, rho[2] = %4.2e \n", rho[0], rho[1], rho[2]);
#endif

  // compute the velocity
  double velocity[3] = {0., 0., 0.};

  for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                          + MLocal[i + 1 * 3] * R[1] +
                                          + MLocal[i + 2 * 3] * R[2] ;

  for (inew = 0 ; inew < itermax ; ++inew) // Newton iteration
  {
    //Update function and gradient

    Function(R, velocity, mu, rho, F, A, B);

/* #ifndef AC_Generated */
#ifndef DEBUG_CHECK
    double AWpB[9];
    if (iparam[10] != 3 && iparam[10] != 4)
    {
      double Fg[3] = {0., 0., 0.};
      double Ag[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
      double Bg[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};



      assert(*rho > 0. && *(rho + 1) > 0. && *(rho + 2) > 0.);

/* #ifdef AC_STD */
      if  (iparam[10] == 0 )
      {
        fc3d_AlartCurnierFunctionGenerated(R, velocity, mu, rho, Fg, Ag, Bg);
      }

/* #endif */

/* #ifdef AC_JeanMoreau */
      if  (iparam[10] == 1 )
      {
        fc3d_AlartCurnierJeanMoreauFunctionGenerated(R, velocity, mu, rho, Fg, Ag, Bg);
      }
/* #endif */

      sub3(F, Fg);
      sub3x3(A, Ag);
      sub3x3(B, Bg);

      assert(hypot3(Fg) <= 1e-7);
      assert(hypot9(Ag) <= 1e-7);
      assert(hypot9(Bg) <= 1e-7);
      cpy3x3(A, Ag);
      cpy3x3(B, Bg);
      mm3x3(A, MLocal, AWpB);
      add3x3(B, AWpB);
    }
#endif
/* #endif */


    // compute -(A MLocal +B)
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        AWplusB[i + 3 * j] = 0.0;
        for (k = 0; k < 3; k++)
        {
          AWplusB[i + 3 * j] -= A[i + 3 * k] * MLocal[k + j * 3];
        }
        AWplusB[i + 3 * j] -= B[i + 3 * j];
      }
    }

/* #ifdef AC_STD */

#ifndef DEBUG_CHECK
    if (iparam[10]==0)
    {
      scal3x3(-1., AWpB);
      sub3x3(AWplusB, AWpB);
      assert(hypot9(AWpB) <= 1e-7);
    }
#endif
/* #endif */

    // Solve the linear system
    if ( solv3x3(AWplusB, dR, F) )
    {
      // if determinant is zero, replace dR=NaN with zero (i.e. don't
      // modify R) and return early
      dR[0] = 0; dR[1] = 0; dR[2] = 0;
      inew = itermax;
      DEBUG_EXPR(
        assert(0 && "solv3x3 returned error, bad determinant found."));
    }

    // upate iterates
    R[0] += dR[0];
    R[1] += dR[1];
    R[2] += dR[2];
    // compute new residue
    for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                            + MLocal[i + 1 * 3] * R[1] +
                                            + MLocal[i + 2 * 3] * R[2] ;
    Function(R, velocity, mu, rho, F, NULL, NULL);
    dparam[1] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) / (1.0 + sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])) ; // improve with relative tolerance

    /*      dparam[2] =0.0;
            fc3d_unitary_compute_and_add_error( R , velocity,mu, &(dparam[2]));*/




    if (verbose > 1) printf("-----------------------------------    fc3d_onecontact_nonsmooth_Newton_solvers_solve_direct number of iteration = %i  error = %.10e \n", inew, dparam[1]);

    if (dparam[1] < Tol)
    {
      /*    printf("-----------------------------------    fc3d_onecontact_nonsmooth_Newton_solvers_solve_direct number of iteration = %i  error = %.10e \t error2 = %.10e \n",inew,dparam[1], dparam[2]); */

      return 0;
    }

  }// End of the Newton iteration

  /*  printf("-----------------------------------    LocalNewtonSolver number of iteration = %i  error = %.10e \t error2 = %.10e \n",inew,dparam[1], dparam[2]); */
  return 1;

}



static int LineSearchGP(FrictionContactProblem* localproblem,
                  computeNonsmoothFunction  Function,
                  double * t_opt,
                  double R[3],
                  double dR[3],
                  double *rho,
                  int LSitermax,
                  double * F,
                  double * A,
                  double * B,
                  double * velocity)
{
  DEBUG_PRINT("LineSearchGP -- Start Line search\n");

  double alpha = *t_opt;

  double inf = 1e20;

  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.1, m2 = 0.9;


  /*     // store the value fo the function */
  /*     double F[3]={0.,0.,0.}; */

  /*     // Store the (sub)-gradient of the function */
  /*     double A[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.}; */
  /*     double B[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.}; */

  /*     double velocity[3]={0.,0.,0.}; */

  double mu = localproblem->mu[0];
  double * qLocal = localproblem->q;
  double * MLocal = localproblem->M->matrix0;

  /*     for (int i=0; i<3; i++) velocity[i] = MLocal[i+0*3]*R[0] + qLocal[i] */
  /*          + MLocal[i+1*3]*R[1] + */
  /*          + MLocal[i+2*3]*R[2] ; */

  /*     Function(R,velocity,mu,rho,F,A,B); */


  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(3 , F , 1 , F , 1);

  double tmp[3] = {0., 0., 0.};

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // compute A MLocal +B
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      AWplusB[i + 3 * j] = 0.0;
      for (int k = 0; k < 3; k++)
      {
        AWplusB[i + 3 * j] += A[i + 3 * k] * MLocal[k + j * 3];
      }
      AWplusB[i + 3 * j] += B[i + 3 * j];
    }
  }

#ifdef DEBUG_MESSAGES
  for (int l = 0; l < 3; l++)
  {
    for (int k = 0; k < 3; k++)
    {
      printf("AWplusB[%i+3*%i] = %le\t", l, k, AWplusB[l + 3 * k]);
    }
    printf("\n");
  }
#endif

  for (int i = 0; i < 3; i++)
  {
    tmp[i] = 0.0;
    for (int j = 0; j < 3; j++)
    {
      tmp[i] += AWplusB[i + 3 * j] * dR[j]  ;
    }
  }




  double dqdt0 = 0.0;
  for (int i = 0; i < 3; i++)
  {
    dqdt0 += F[i] * tmp[i];
  }
#ifdef DEBUG_MESSAGES
  printf("q0 = %12.8e \n", q0);
  printf("dqdt0 = %12.8e \n", dqdt0);
  for (int i = 0; i < 3; i++)
  {
    printf("tmp[%i] = %12.8e \t", i, tmp[i]);
  }
  printf("\n");
  for (int i = 0; i < 3; i++)
  {
    printf("dR[%i] = %12.8e \t", i, dR[i]);
  }
  printf("\n");
#endif

  for (int iter = 0; iter < LSitermax; iter++)
  {

    for (int i = 0; i < 3; i++)  tmp[i] = R[i] + alpha * dR[i];

    for (int i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * tmp[0] + qLocal[i]
          + MLocal[i + 1 * 3] * tmp[1] +
          + MLocal[i + 2 * 3] * tmp[2] ;

    Function(tmp, velocity, mu, rho, F, NULL, NULL);

    double q  = 0.5 * cblas_ddot(3 , F , 1 , F , 1);

    double slope = (q - q0) / alpha;

#ifdef DEBUG_MESSAGES
    printf("q = %12.8e \n", q);
    printf("slope = %12.8e \n", slope);
#endif


    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2)
    {
      DEBUG_PRINTF("Success in LS: alpha = %12.8e\n", alpha);
      *t_opt = alpha;
      if (verbose > 1)
      {
        printf("-----------------------------------------    LineSearchGP success number of iteration = %i  alpha = %.10e \n", iter, alpha);
      }
      return 0;

    }
    else if (!C1)
    {
#ifdef DEBUG_MESSAGES
      printf("LS: alpha too small = %12.8e\t, slope =%12.8e\n", alpha, slope);
      printf(" m1*dqdt0 =%12.8e\t, m2*dqdt0 =%12.8e\n ", m1 * dqdt0 , m2 * dqdt0);
#endif
      //std::cout << "t = " << t << " is too small : slope = " << slope << ", m2*qp0 = " << m2*qp0 << std::endl;
      alphamin = alpha;
    }
    else   // not(C2)
    {
#ifdef DEBUG_MESSAGES
      printf("LS: alpha too big = %12.8e\t, slope =%12.8e\n", alpha, slope);
      printf(" m1*dqdt0 =%12.8e\t, m2*dqdt0 =%12.8e\n ", m1 * dqdt0 , m2 * dqdt0);
#endif
      //std::cout << "t = " << t << " is too big : slope = " << slope << ", m1*qp0 = " << m1*qp0 << std::endl;
      alphamax = alpha;
    }
    if (alpha < inf)
    {
      alpha = 0.5 * (alphamin + alphamax);
    }
    else
    {
      alpha = 10 * alpha;
    }


  }
  if (verbose > 1)
  {
    printf("-----------------------------------------    LineSearchGP failed max number of iteration reached  = %i  with alpha = %.10e \n", LSitermax, alpha);
  }
  *t_opt = alpha;
  return -1;
}

double fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped_oneitereval(
  FrictionContactProblem* localproblem, double * R, int *iparam,
  double *dparam)
{
  double mu = localproblem->mu[0];
  double rho[3] = {1., 1., 1.};
  double velocity[3] = {0., 0., 0.};
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double F[3] = {0., 0., 0.};
  double * qLocal = localproblem->q;
  double * MLocal = localproblem->M->matrix0;
  for (int i = 0; i < 3; i++)
    velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                + MLocal[i + 1 * 3] * R[1]
                + MLocal[i + 2 * 3] * R[2] ;
  Function(R, velocity, mu, rho, F, A, B);
  return 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2])
    / (1.0 + sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]));
}

int fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped(FrictionContactProblem* localproblem, double * R, int *iparam, double *dparam)
{
  DEBUG_PRINT("fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped() starts \n");
  DEBUG_EXPR(verbose=3;);

  assert(localproblem);
  assert(localproblem->q);
  assert(localproblem->mu);
  assert(localproblem->M);
  assert(localproblem->M->matrix0);

  double mu = localproblem->mu[0];
  double * qLocal = localproblem->q;
  double * MLocal = localproblem->M->matrix0;


  double Tol = dparam[0];
  int itermax = iparam[0];
  int LSitermax = iparam[12];


  int i, j, k, inew;

  // store the value fo the function
  double F[3] = {0., 0., 0.};

  // Store the (sub)-gradient of the function
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // store the search direction
  double dR[3] = {0., 0., 0.};

  // path length
  double t = 1.;
  double t_opt = 1.;
  double t_init = 1.;
  int NumberofLSfailed = 0;

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Compute values of Rho (should be here ?)
  double rho[3] = {1., 1., 1.};
#ifdef OPTI_RHO
//  computerho(localproblem, rho);
  DEBUG_PRINTF("rho[0] = %4.2e, rho[1] = %4.2e, rho[2] = %4.2e \n", rho[0], rho[1], rho[2]);
#endif

  // compute the velocity
  double velocity[3] = {0., 0., 0.};
  //cpy3(qLocal,velocity);
  //mvp3x3(MLocal,velocity)

  for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                          + MLocal[i + 1 * 3] * R[1] +
                                          + MLocal[i + 2 * 3] * R[2] ;

  DEBUG_EXPR_WE(for (int i =0 ; i < 3; i++) printf("R[%i]= %12.8e,\t velocity[%i]= %12.8e,\n",i,R[i],i,velocity[i]););
  DEBUG_PRINT("fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped -- Start Newton iteration\n");
  assert(Function);

  for (inew = 0 ; inew < itermax ; ++inew) // Newton iteration
  {
    //Update function and gradient
    Function(R, velocity, mu, rho, F, A, B);

    DEBUG_EXPR_WE(for ( int i =0 ; i < 3; i++) printf("F[%i]=%12.8e\t",i,F[i]); printf("\n"););

    // compute -(A MLocal +B)
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        AWplusB[i + 3 * j] = 0.0;
        for (k = 0; k < 3; k++)
        {
          AWplusB[i + 3 * j] -= A[i + 3 * k] * MLocal[k + j * 3];
        }
        AWplusB[i + 3 * j] -= B[i + 3 * j];
      }
    }

    if ( solv3x3(AWplusB, dR, F) )
    {
      // if determinant is zero, replace dR=NaN with zero (i.e. don't
      // modify R) and return early
      dR[0] = 0; dR[1] = 0; dR[2] = 0;
      inew = itermax;
      DEBUG_EXPR(
        assert(0 && "solv3x3 returned error, bad determinant found."));
    }
    else
    {
      // Perform Line Search

      t_opt = t_init;
      int infoLS = LineSearchGP(localproblem, Function, &t_opt, R, dR,
                                rho, LSitermax, F, A, B, velocity);

      if (infoLS == 0)
        t = t_opt;
      else
      {
        NumberofLSfailed++;
        if (NumberofLSfailed > 5)
        {
          t = 100.0;
          if (verbose > 1)
            printf("-----------------------------------------  "
                   "Max Number of LineSearchGP failed =%i Tilt point\n ",
                   NumberofLSfailed);
          NumberofLSfailed = 0;
        }
      }
    }

    // upate iterates
    R[0] = R[0] + t * dR[0];
    R[1] = R[1] + t * dR[1];
    R[2] = R[2] + t * dR[2];

    // compute new residue
    for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                            + MLocal[i + 1 * 3] * R[1] +
                                            + MLocal[i + 2 * 3] * R[2] ;

    Function(R, velocity, mu, rho, F, NULL, NULL);
    dparam[1] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) / (1.0 + sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])) ; // improve with relative tolerance

    if (verbose > 1) printf("-----------------------------------  fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped.  number of iteration = %i  error = %.10e \n", inew, dparam[1]);
    if (dparam[1] < Tol) return 0;


  }// End of the Newton iteration


  return 1;

}
