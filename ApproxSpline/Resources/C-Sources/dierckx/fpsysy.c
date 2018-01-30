/* fpsysy.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/*<       subroutine fpsysy(a,n,g) >*/
/* Subroutine */ int fpsysy_(doublereal *a, integer *n, doublereal *g)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, i1;
    static doublereal fac;

/* subroutine fpsysy solves a linear n x n symmetric system */
/*    (a) * (b) = (g) */
/* on input, vector g contains the right hand side ; on output it will */
/* contain the solution (b). */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer n >*/
/*  ..array arguments.. */
/*<       real a(6,6),g(6) >*/
/*  ..local scalars.. */
/*<       real fac >*/
/*<       integer i,i1,j,k >*/
/*  .. */
/*<       g(1) = g(1)/a(1,1) >*/
    /* Parameter adjustments */
    --g;
    a -= 7;

    /* Function Body */
    g[1] /= a[7];
/*<       if(n.eq.1) return >*/
    if (*n == 1) {
	return 0;
    }
/*  decomposition of the symmetric matrix (a) = (l) * (d) *(l)' */
/*  with (l) a unit lower triangular matrix and (d) a diagonal */
/*  matrix */
/*<       do 10 k=2,n >*/
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
/*<          a(k,1) = a(k,1)/a(1,1) >*/
	a[k + 6] /= a[7];
/*<   10  continue >*/
/* L10: */
    }
/*<       do 40 i=2,n >*/
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          i1 = i-1 >*/
	i1 = i__ - 1;
/*<          do 30 k=i,n >*/
	i__2 = *n;
	for (k = i__; k <= i__2; ++k) {
/*<             fac = a(k,i) >*/
	    fac = a[k + i__ * 6];
/*<             do 20 j=1,i1 >*/
	    i__3 = i1;
	    for (j = 1; j <= i__3; ++j) {
/*<                fac = fac-a(j,j)*a(k,j)*a(i,j) >*/
		fac -= a[j + j * 6] * a[k + j * 6] * a[i__ + j * 6];
/*<   20        continue >*/
/* L20: */
	    }
/*<             a(k,i) = fac >*/
	    a[k + i__ * 6] = fac;
/*<             if(k.gt.i) a(k,i) = fac/a(i,i) >*/
	    if (k > i__) {
		a[k + i__ * 6] = fac / a[i__ + i__ * 6];
	    }
/*<   30     continue >*/
/* L30: */
	}
/*<   40  continue >*/
/* L40: */
    }
/*  solve the system (l)*(d)*(l)'*(b) = (g). */
/*  first step : solve (l)*(d)*(c) = (g). */
/*<       do 60 i=2,n >*/
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          i1 = i-1 >*/
	i1 = i__ - 1;
/*<          fac = g(i) >*/
	fac = g[i__];
/*<          do 50 j=1,i1 >*/
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
/*<             fac = fac-g(j)*a(j,j)*a(i,j) >*/
	    fac -= g[j] * a[j + j * 6] * a[i__ + j * 6];
/*<   50     continue >*/
/* L50: */
	}
/*<          g(i) = fac/a(i,i) >*/
	g[i__] = fac / a[i__ + i__ * 6];
/*<   60  continue >*/
/* L60: */
    }
/*  second step : solve (l)'*(b) = (c) */
/*<       i = n >*/
    i__ = *n;
/*<       do 80 j=2,n >*/
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/*<          i1 = i >*/
	i1 = i__;
/*<          i = i-1 >*/
	--i__;
/*<          fac = g(i) >*/
	fac = g[i__];
/*<          do 70 k=i1,n >*/
	i__2 = *n;
	for (k = i1; k <= i__2; ++k) {
/*<             fac = fac-g(k)*a(k,i) >*/
	    fac -= g[k] * a[k + i__ * 6];
/*<   70     continue >*/
/* L70: */
	}
/*<          g(i) = fac >*/
	g[i__] = fac;
/*<   80  continue >*/
/* L80: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpsysy_ */

