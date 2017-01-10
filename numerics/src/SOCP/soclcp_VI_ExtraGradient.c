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

#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SecondOrderConeLinearComplementarityProblem_as_VI.h"
#include "VariationalInequality_Solvers.h"
#include "SOCLCP_Solvers.h"
#include "soclcp_compute_error.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"




void soclcp_VI_ExtraGradient(SecondOrderConeLinearComplementarityProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{

  /* Dimension of the problem */
  int n = problem->n;


  VariationalInequality *vi = (VariationalInequality *)malloc(sizeof(VariationalInequality));

  //vi.self = &vi;
  vi->F = &Function_VI_SOCLCP;
  vi->ProjectionOnX = &Projection_VI_SOCLCP;

  int iter=0;
  double error=1e24;

  SecondOrderConeLinearComplementarityProblem_as_VI *soclcp_as_vi= (SecondOrderConeLinearComplementarityProblem_as_VI*)malloc(sizeof(SecondOrderConeLinearComplementarityProblem_as_VI));
  vi->env =soclcp_as_vi ;
  vi->size =  n;


  /*Set the norm of the VI to the norm of problem->q  */
  vi->normVI= cblas_dnrm2(n , problem->q , 1);
  vi->istheNormVIset=1;

  soclcp_as_vi->vi = vi;
  soclcp_as_vi->soclcp = problem;
  /* soclcp_display(fc3d_as_vi->fc3d); */

  SolverOptions * visolver_options = (SolverOptions *) malloc(sizeof(SolverOptions));
  variationalInequality_setDefaultSolverOptions(visolver_options,
                                                SICONOS_VI_EG);

  visolver_options->params = options->params;

  variationalInequality_ExtraGradient(vi, reaction, velocity , info , visolver_options);



  /* **** Criterium convergence **** */
  soclcp_compute_error(problem, reaction, velocity, options, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  error = visolver_options->params.common.residu;
  iter = visolver_options->params.common.extra_iter_done;

  options->params.common.residu = error;
  options->params.common.extra_iter_done = iter;


  if (verbose > 0)
  {
    printf("----------------------------------- SOCLCP - VI Extra Gradient (VI_EG) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  free(vi);

  solver_options_delete(visolver_options);
  free(visolver_options);
  visolver_options=NULL;
  free(soclcp_as_vi);



}


int soclcp_VI_ExtraGradient_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ExtraGradient Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_SOCLCP_VI_EG;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->dWork = NULL;
  solver_options_nullify(options);

  memset(&options->params, 0, sizeof(options->params));

  options->params.common.max_iter = 20000;
  options->params.common.tolerance = 1e-3;
  options->params.common.rho = 1e-3;
  options->params.common.rho = -1.0; // rho is variable by default
  options->internalSolvers = NULL;

  return 0;
}
