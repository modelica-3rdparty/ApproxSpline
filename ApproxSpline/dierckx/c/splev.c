/* splev.f -- translated by f2c (version 20061008).
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

/*<       subroutine splev(t,n,c,k,x,y,m,ier) >*/
/* Subroutine */ int splev_(doublereal *t, integer *n, doublereal *c__, 
	integer *k, doublereal *x, doublereal *y, integer *m, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, k1, l1;
    static doublereal tb;
    static integer ll;
    static doublereal te, sp;
    static integer nk1;
    static doublereal arg;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  subroutine splev evaluates in a number of points x(i),i=1,2,...,m */
/*  a spline s(x) of degree k, given in its b-spline representation. */

/*  calling sequence: */
/*     call splev(t,n,c,k,x,y,m,ier) */

/*  input parameters: */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length n, which contains the b-spline coefficients. */
/*    k    : integer, giving the degree of s(x). */
/*    x    : array,length m, which contains the points where s(x) must */
/*           be evaluated. */
/*    m    : integer, giving the number of points where s(x) must be */
/*           evaluated. */

/*  output parameter: */
/*    y    : array,length m, giving the value of s(x) at the different */
/*           points. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    m >= 1 */
/*    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1. */

/*  other subroutines required: fpbspl. */

/*  references : */
/*    de boor c  : on calculating with b-splines, j. approximation theory */
/*                 6 (1972) 50-62. */
/*    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths */
/*                 applics 10 (1972) 134-149. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/*  ..scalar arguments.. */
/*<       integer n,k,m,ier >*/
/*  ..array arguments.. */
/*<       real t(n),c(n),x(m),y(m) >*/
/*  ..local scalars.. */
/*<       integer i,j,k1,l,ll,l1,nk1 >*/
/*<       real arg,sp,tb,te >*/
/*  ..local array.. */
/*<       real h(6) >*/
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
/*<       ier = 10 >*/
    /* Parameter adjustments */
    --c__;
    --t;
    --y;
    --x;

    /* Function Body */
    *ier = 10;
/*<       if(m-1) 100,30,10 >*/
    if ((i__1 = *m - 1) < 0) {
	goto L100;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L10;
    }
/*<   10  do 20 i=2,m >*/
L10:
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         if(x(i).lt.x(i-1)) go to 100 >*/
	if (x[i__] < x[i__ - 1]) {
	    goto L100;
	}
/*<   20  continue >*/
/* L20: */
    }
/*<   30  ier = 0 >*/
L30:
    *ier = 0;
/*  fetch tb and te, the boundaries of the approximation interval. */
/*<       k1 = k+1 >*/
    k1 = *k + 1;
/*<       nk1 = n-k1 >*/
    nk1 = *n - k1;
/*<       tb = t(k1) >*/
    tb = t[k1];
/*<       te = t(nk1+1) >*/
    te = t[nk1 + 1];
/*<       l = k1 >*/
    l = k1;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*  main loop for the different points. */
/*<       do 80 i=1,m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch a new x-value arg. */
/*<         arg = x(i) >*/
	arg = x[i__];
/*<         if(arg.lt.tb) arg = tb >*/
	if (arg < tb) {
	    arg = tb;
	}
/*<         if(arg.gt.te) arg = te >*/
	if (arg > te) {
	    arg = te;
	}
/*  search for knot interval t(l) <= arg < t(l+1) */
/*<   40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50 >*/
L40:
	if (arg < t[l1] || l == nk1) {
	    goto L50;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 40 >*/
	goto L40;
/*  evaluate the non-zero b-splines at arg. */
/*<   50    call fpbspl(t,n,k,arg,l,h) >*/
L50:
	fpbspl_(&t[1], n, k, &arg, &l, h__);
/*  find the value of s(x) at x=arg. */
/*<         sp = 0. >*/
	sp = 0.;
/*<         ll = l-k1 >*/
	ll = l - k1;
/*<         do 60 j=1,k1 >*/
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
/*<           ll = ll+1 >*/
	    ++ll;
/*<           sp = sp+c(ll)*h(j) >*/
	    sp += c__[ll] * h__[j - 1];
/*<   60    continue >*/
/* L60: */
	}
/*<         y(i) = sp >*/
	y[i__] = sp;
/*<   80  continue >*/
/* L80: */
    }
/*<  100  return >*/
L100:
    return 0;
/*<       end >*/
} /* splev_ */

