#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*!\file rd_nlgs.c



*/



double ddot_(int *, double [], int *, double [], int*);


/*!\fn  rd_nlgs(double vec[],double *qq,int *nn,double a[],int * itermax, double * tol,double z[],double w[],int *it_end,double * res,int *info)

   rd_nlgs is a specific nlgs (Gauss Seidel Non Linear)solver for relay problems.

   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param a On enter a pointer over doubles, the bound.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param tol On enter a pointer over doubles, the tolerance required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
 */


rd_nlgs(double vec[], double *q, int *nn, double a[], double b[], int * itermax, double * tol, double z[], double w[], int *it_end, double * res, int *info)
{

  FILE *f101;
  int i, j, iter1, k;
  int n = *nn, incx = 1, incy = 1, itt = *itermax;
  double errmax = *tol, alpha, beta, mina;
  double err1, num, den, avn, avt, apn, apt, xn;
  double *zt, *wnum1;
  char trans = 'T';
  /*  double M[n][n];*/
  double(*M)[n];


  f101 = fopen("resultat_nlgs.dat", "w+");

  M =  malloc(n * n * sizeof(double));

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      M[i][j] = vec[i * n + j];

  wnum1 = (double*) malloc(n * sizeof(double));
  zt = (double*) malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    w[i] = 0.;
    z[i] = 0.;
    zt[i] = 0.;
  }

  iter1 = 1;
  err1 = 1.;

  while ((iter1 < itt) && (err1 > errmax))
  {
    iter1 = iter1 + 1;
    /*  dcopy_( &n, q, &incx, wnum1, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_( &trans, &n, &n, &alpha, M, &n, z, &incx, &beta, wnum1, &incy);*/

    for (i = 0; i < n; i++)
    {
      avn = 0.;
      apn = 0.;
      for (j = 0; j <= i - 1; j++)
        avn = avn + M[i][j] * z[j];
      for (k = i + 1; k < n; k++)
        apn = apn + M[i][k] * z[k];
      xn = q[i] - avn - apn;
      /*      zt[i] = 1/M[i][i] *xn;*/
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
      z[i] = 1 / M[i][i] * (w[i] + xn);


    }

    /* ///////// convergence criterium ///////////// */

    dcopy_(&n, w, &incx, wnum1, &incy);
    alpha = 1.;
    daxpy_(&n, &alpha, q, &incx, wnum1, &incy);
    alpha = 1.;
    beta = -1.;
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, wnum1, &incy);
    num = ddot_(&n, wnum1, &incx, wnum1, &incy);
    den = ddot_(&n, q, &incx, q, &incy);

    err1 = sqrt(num) / sqrt(den);
    *it_end = iter1;
    *res = err1;

    for (i = 0; i < n; i++)
    {
      /*result_gs[i][iter1-1] = z[i]; */
      fprintf(f101, "%d  %d  %14.7e %14.7e\n", iter1 - 1, i, z[i], w[i]);
    }



  }


  if (err1 > errmax)
  {
    printf("no convergence after %d iterations, the residue is %g\n", iter1, err1);
    *info = 1;
  }
  else
  {
    printf("there is convergence after %d iterations, the residue is %g \n", iter1, err1);
    *info = 0;
  }

  free(wnum1);
  free(zt);
  free(M);
  fclose(f101);

}
