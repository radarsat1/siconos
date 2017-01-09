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

#include "PathSearch.h"

#include <assert.h>

#include "SolverOptions.h"

void free_solverData_PathSearch(void* solverData)
{
  assert(solverData);
  pathsearch_data* solverData_PathSearch = (pathsearch_data*) solverData;
  free_NMS_data(solverData_PathSearch->data_NMS);
  free(solverData_PathSearch->lsa_functions);
}

void pathsearch_default_SolverOption(SolverOptions* options)
{
  options->params.line_search.nonmonotone_ls = NM_LS_MEAN;
  options->params.line_search.nonmonotone_ls_m = 10;
  options->params.line_search.path.stack_size = 5;
  options->params.line_search.nm.watchdog_type = LINESEARCH;
  options->params.line_search.nm.projected_gradient_type = ARCSEARCH;
  options->params.line_search.nm.n_max = 10;
  options->params.line_search.nm.delta = 20;
  options->params.line_search.nm.delta_var = 0.8;
  options->params.line_search.nm.sigma = 0.01;
  options->params.line_search.nm.alpha_min_watchdog = 1e-12;
  options->params.line_search.nm.alpha_min_pgrad = 1e-12;
  options->params.line_search.nm.merit_incr = 1.1;
}
