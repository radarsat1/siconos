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



#ifndef SolverOptions_H
#define SolverOptions_H

/*!\file SolverOptions.h
  Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
  \author Franck Perignon
*/
#include "SiconosConfig.h"
#include "NumericsFwd.h"

/** \struct Callback SolverOptions.h
Structure used to store user callbacks inside solvers
*/
typedef struct
{
  void *env; /**< general user environment */
  void (*collectStatsIteration)(void *env, int size, double*reaction,
                       double*velocity, double error, void* extra_data);/**< pointer on a function
* Its signature is: user env, problem size, reaction,
* velocity, error at end of solver iteration (when this makes sense) and an
* extra data structure */
} Callback;

/* Solver Options Parameters */

struct SiconosSolverCommonParams
{
  int max_iter;
  int iter_done;
  int prealloc;
  double tolerance;
  double residu;
};

struct SiconosSolverPivotBasedParams
{
  struct SiconosSolverCommonParams common;
  int pivot_rule;
};

struct SiconosSolverPathSearchParams
{
  struct SiconosSolverCommonParams common;
  int stack_size;
};

struct SiconosGoldsteinParams
{
  int itermax;
  double c;
  double alphamax;
};

struct SiconosNonMonotoneParams
{
  int watchdog_type;
  int projected_gradient_type;
  int n_max;
  double delta;
  double delta_var;
  double sigma;
  double alpha_min_watchdog;
  double alpha_min_pgrad;
  double merit_incr;
};

struct SiconosLineSearchParams
{
  struct SiconosSolverCommonParams common;
  int nonmonotone_ls;
  int nonmonotone_ls_m;
  int force_arcsearch;
  int criterion;
  double alpha_min;
  struct SiconosNonMonotoneParams nm;
  struct SiconosGoldsteinParams goldstein;
};

struct SolverOptionsParams
{
  union {
    struct SiconosSolverCommonParams common;
    struct SiconosSolverPivotBasedParams pivot_based;
    struct SiconosSolverPathSearchParams path_search;
    struct SiconosLineSearchParams line_search;
  };
};

/** \struct SolverOptions_ SolverOptions.h
    Structure used to send options (name, parameters and so on) to a specific solver (mainly from Kernel to Numerics).
*/
struct SolverOptions
{
  int solverId;                            /**< solverId Id of the solver (see ) */
  int isSet;                               /**< isSet int equal to false(0) if the parameters below have not been set (ie need to read default values) else true(1)*/
  struct SolverOptionsParams params;       /**< solver parameters */
  int filterOn;                            /**< filterOn 1 to check solution validity after the driver call, else 0. Default = 1. (For example if
                                            * filterOn = 1 for a LCP, lcp_compute_error() will be called at the end of the process) */
  int dWorkSize;                           /**< dWorkSize size of vector iWork */
  double * dWork;                          /**< dWork is a pointer on a working memory zone (for doubles) reserved for the solver .*/
  int iWorkSize;                           /**< iWorkSize size of vector iWork */
  int * iWork;                             /**< iWork is a pointer on a working memory zone (for integers) reserved for the solver .*/
  int numberOfInternalSolvers;             /**< numberOfInternalSolvers the number of internal or local 'sub-solvers' used by the solver*/
  struct SolverOptions * internalSolvers; /**< internalSolvers pointer to sub-solvers*/
  Callback * callback;                     /**< callback a pointer to user Callback*/

  void * solverParameters;                 /**< additional parameters specific to the solver */

  void * solverData;                       /**< additional data specific to the solver */

};

enum SICONOS_NUMERICS_PROBLEM_TYPE
{
  SICONOS_NUMERICS_PROBLEM_LCP = 0,
  SICONOS_NUMERICS_PROBLEM_MLCP = 1,
  SICONOS_NUMERICS_PROBLEM_EQUALITY = 2,
  SICONOS_NUMERICS_PROBLEM_FC2D = 3,
  SICONOS_NUMERICS_PROBLEM_FC3D = 4,
  SICONOS_NUMERICS_PROBLEM_NCP = 5,
  SICONOS_NUMERICS_PROBLEM_MCP = 6,
  SICONOS_NUMERICS_PROBLEM_VI = 7,
  SICONOS_NUMERICS_PROBLEM_AVI = 8
};


extern const char* const SICONOS_NUMERICS_PROBLEM_LCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_MLCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_NCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_MCP_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_EQUALITY_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_FC2D_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_FC3D_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_VI_STR;
extern const char* const SICONOS_NUMERICS_PROBLEM_AVI_STR;


#include "SolverOptions_helpers.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** screen display of solver parameters
      \param options the structure to be displayed
  */
  void solver_options_print(SolverOptions* options);

  /** free some SolverOptions fields;
   *   \param options the structure to clean
   */
  void solver_options_delete(SolverOptions * options);

  /* Set all pointer fields to NULL, except iparam and dparam
   * \param options the struct to initialize
   */
  void solver_options_nullify(SolverOptions* options);

  /** fill a SolverOptions struct: set fields, allocate memory and set common
   * values
   * \param options struct to fill
   * \param solverId identity of the solver
   * \param iSize size of the iparam field (integer parameters)
   * \param dSize size of the dparam field (double parameters)
   * \param iter_max maximum number of iterations before the solver stops
   * if this does not make sense or is unwanted, give inf as value
   * \param tol tolerance for the solution.
   * if this does not make sense or is unwanted, give inf as value
   */
  void solver_options_fill(SolverOptions* options, int solverId, int iSize, int dSize, int iter_max, double tol);

  /** set parameters in SolverOption. This function should be used instead of
   * rewrittent each time a new function for setting the parameters
   * \param options the struct to set
   * \param solverId the id of the solver
   */
  void solver_options_set(SolverOptions* options, int solverId);

  /** return the id of a solver based on its name
   * \param pName the name of the solver
   * \return the id of the solver or 0 if it failed
   */
  int solver_options_name_to_id(char * pName);

  /** return the name of a solver given its id
   * \param Id the id of the solver
   * \return the name of the solver
   */
  const char * solver_options_id_to_name(int Id);

  /** return the name of a problem type (LCP, NCP, VI, ...) based on its id
   * \param id the id of the problem
   * \return the name of the problem
   */
  const char * ns_problem_id_to_name(int id);

  /** free the solverData structure
   * \param options the structure to free
   */
  void solver_options_free_solver_specific_data(SolverOptions* options);

  /** copy SolverOptions
   * \param options_ori the structure to copy
   * \param options the output structure 
   */
  void solver_options_copy(SolverOptions* options_ori, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif



#endif
