/* Siconos-Numerics version 1.1.3, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

/*!\file dr_nlgs.c


This subroutine allows the resolution of Dual Relay (DR) problem.\n
Try \f$(z,w)\f$ such that:\n

\f$
\left\lbrace
\begin{array}{l}
w - M z = q\\
-z \in \partial\psi_{[-b,a]}(w)\\
\end{array}
\right.
\f$

here M is an (nn \f$\times\f$nn)- matrix, q an nn-dimensional vector, w, z, a and b are nn-dimensional vectors.

*/

/*!\fn void dr_nlgs( double *vec, double *q, int *nn, double *a, double *b, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)

   dr_nlgs is a specific nlgs (Non Linear Gauss Seidel) solver for dual relay problems.\n

   \param vec      On enter, a (nn \f$\times\f$nn)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
   \param q        On enter, a nn-vector of doubles containing the components of the double vector.
   \param nn       On enter, an integer, the dimension of the second member.
   \param a        On enter, a nn-vector of doubles, the upper bound.
   \param b        On enter, a nn-vector of doubles, the lower bound.
   \param itermax  On enter, an integer, the maximum iterations required.
   \param tol      On enter, a double, the tolerance required.
   \param chat     On enter, an integer, the output log identifiant:\n
                    0  =  no output, \n
                    >0 =  active screen output\n

   \param it_end   On return, an integer, the number of iterations carried out.
   \param res      On return, a double, the error value.
   \param z        On return, a nn-vector of doubles, the solution of the problem.
   \param w        On return, a nn-vector of doubles, the solution of the problem.
   \param info     On return, an integer, the termination reason:\n
                    0 = convergence,\n
        1 = no convergence,\n
        2 = nul diagonal term.\n


   \author Nineb Sheherazade.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blaslapack.h"




void dr_nlgs(double *vec, double *q, int *nn, double *a, double *b, int * itermax, double * tol, int *chat, double *z, double *w, int *it_end, double * res, int *info)
{



  int i, j, iter1, k, ispeak = *chat;
  int n = *nn,  itt = *itermax;
  integer incx = 1, incy = 1;
  double errmax = *tol, alpha, beta, mina;
  double err1, num, den, avn, xn, apn;
  double *zt, *wnum1;

  char trans = 'N';





  wnum1    = (double*) malloc(n * sizeof(double));
  zt       = (double*) malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    w[i]     = 0.;
    z[i]     = 0.;
    zt[i]    = 0.;
    wnum1[i] = 0.;
  }


  iter1 = 1;
  err1 = 1.;


  while ((iter1 < itt) && (err1 > errmax))
  {

    iter1 = iter1 + 1;

    for (i = 0; i < n; i++)
    {
      avn = 0.;
      apn = 0.;

      for (j = 0; j <= i - 1; j++)
        avn = avn + vec[j * n + i] * z[j];

      for (k = i + 1; k < n; k++)
        apn = apn + vec[k * n + i] * z[k];

      xn = -q[i] - avn - apn;

      zt[i] = -xn;

      if (a[i] < zt[i])
      {
        mina = a[i];
      }
      else
      {
        mina = zt[i];
      }

      if (-b[i] < mina)
      {
        w[i] = mina;
      }
      else
      {
        w[i] = -b[i];
      }

      if (fabs(vec[i * n + i]) < 1e-12)
      {
        printf("\n Warning nul diagonal term of M \n");

        free(zt);
        free(wnum1);

        *info = 2;

        return;

      }
      else
        z[i] = 1 / vec[i * n + i] * (w[i] + xn);


    }

    /*              Convergence criterium              */

    dcopy_((integer *)&n, w, &incx, wnum1, &incy);

    alpha = -1.;
    daxpy_((integer *)&n, &alpha, q, &incx, wnum1, &incy);

    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, (integer *)&n, (integer *)&n, &alpha, vec, (integer *)&n, z, &incx, &beta, wnum1, &incy);

    num = ddot_((integer *)&n, wnum1, &incx, wnum1, &incy);

    den = ddot_((integer *)&n, q, &incx, q, &incy);

    err1 = sqrt(num) / sqrt(den);
    *it_end = iter1;
    *res = err1;




  }


  if (err1 > errmax)
  {
    if (ispeak > 0)
      printf("No convergence after %d iterations, the residue is %g\n", iter1, err1);

    *info = 1;
  }
  else
  {
    if (ispeak > 0)
      printf("Convergence after %d iterations, the residue is %g \n", iter1, err1);

    *info = 0;
  }



  free(wnum1);
  free(zt);



}
