/* splder.f -- translated by f2c (version 20061008).
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

/*<       subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier) >*/
/* Subroutine */ int splder_(doublereal *t, integer *n, doublereal *c__, 
	integer *k, integer *nu, doublereal *x, doublereal *y, integer *m, 
	doublereal *wrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, k1, k2, l1, l2;
    static doublereal ak;
    static integer kk;
    static doublereal tb;
    static integer ll;
    static doublereal te;
    static integer nn;
    static doublereal sp;
    static integer nk1, nk2;
    static doublereal fac, arg;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  subroutine splder evaluates in a number of points x(i),i=1,2,...,m */
/*  the derivative of order nu of a spline s(x) of degree k,given in */
/*  its b-spline representation. */

/*  calling sequence: */
/*     call splder(t,n,c,k,nu,x,y,m,wrk,ier) */

/*  input parameters: */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length n, which contains the b-spline coefficients. */
/*    k    : integer, giving the degree of s(x). */
/*    nu   : integer, specifying the order of the derivative. 0<=nu<=k */
/*    x    : array,length m, which contains the points where the deriv- */
/*           ative of s(x) must be evaluated. */
/*    m    : integer, giving the number of points where the derivative */
/*           of s(x) must be evaluated */
/*    wrk  : real array of dimension n. used as working space. */

/*  output parameters: */
/*    y    : array,length m, giving the value of the derivative of s(x) */
/*           at the different points. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    0 <= nu <= k */
/*    m >= 1 */
/*    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1. */

/*  other subroutines required: fpbspl */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths */
/*                applics 10 (1972) 134-149. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/*  ..scalar arguments.. */
/*<       integer n,k,nu,m,ier >*/
/*  ..array arguments.. */
/*<       real t(n),c(n),x(m),y(m),wrk(n) >*/
/*  ..local scalars.. */
/*<       integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn >*/
/*<       real ak,arg,fac,sp,tb,te >*/
/*  ..local arrays .. */
/*<       real h(6) >*/
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
/*<       ier = 10 >*/
    /* Parameter adjustments */
    --wrk;
    --c__;
    --t;
    --y;
    --x;

    /* Function Body */
    *ier = 10;
/*<       if(nu.lt.0 .or. nu.gt.k) go to 200 >*/
    if (*nu < 0 || *nu > *k) {
	goto L200;
    }
/*<       if(m-1) 200,30,10 >*/
    if ((i__1 = *m - 1) < 0) {
	goto L200;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L10;
    }
/*<   10  do 20 i=2,m >*/
L10:
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         if(x(i).lt.x(i-1)) go to 200 >*/
	if (x[i__] < x[i__ - 1]) {
	    goto L200;
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
/*  the derivative of order nu of a spline of degree k is a spline of */
/*  degree k-nu,the b-spline coefficients wrk(i) of which can be found */
/*  using the recurrence scheme of de boor. */
/*<       l = 1 >*/
    l = 1;
/*<       kk = k >*/
    kk = *k;
/*<       nn = n >*/
    nn = *n;
/*<       do 40 i=1,nk1 >*/
    i__1 = nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          wrk(i) = c(i) >*/
	wrk[i__] = c__[i__];
/*<   40  continue >*/
/* L40: */
    }
/*<       if(nu.eq.0) go to 100 >*/
    if (*nu == 0) {
	goto L100;
    }
/*<       nk2 = nk1 >*/
    nk2 = nk1;
/*<       do 60 j=1,nu >*/
    i__1 = *nu;
    for (j = 1; j <= i__1; ++j) {
/*<          ak = kk >*/
	ak = (doublereal) kk;
/*<          nk2 = nk2-1 >*/
	--nk2;
/*<          l1 = l >*/
	l1 = l;
/*<          do 50 i=1,nk2 >*/
	i__2 = nk2;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             l1 = l1+1 >*/
	    ++l1;
/*<             l2 = l1+kk >*/
	    l2 = l1 + kk;
/*<             fac = t(l2)-t(l1) >*/
	    fac = t[l2] - t[l1];
/*<             if(fac.le.0.) go to 50 >*/
	    if (fac <= 0.) {
		goto L50;
	    }
/*<             wrk(i) = ak*(wrk(i+1)-wrk(i))/fac >*/
	    wrk[i__] = ak * (wrk[i__ + 1] - wrk[i__]) / fac;
/*<   50     continue >*/
L50:
	    ;
	}
/*<          l = l+1 >*/
	++l;
/*<          kk = kk-1 >*/
	--kk;
/*<   60  continue >*/
/* L60: */
    }
/*<       if(kk.ne.0) go to 100 >*/
    if (kk != 0) {
	goto L100;
    }
/*  if nu=k the derivative is a piecewise constant function */
/*<       j = 1 >*/
    j = 1;
/*<       do 90 i=1,m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          arg = x(i) >*/
	arg = x[i__];
/*<   70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80 >*/
L70:
	if (arg < t[l + 1] || l == nk1) {
	    goto L80;
	}
/*<          l = l+1 >*/
	++l;
/*<          j = j+1 >*/
	++j;
/*<          go to 70 >*/
	goto L70;
/*<   80     y(i) = wrk(j) >*/
L80:
	y[i__] = wrk[j];
/*<   90  continue >*/
/* L90: */
    }
/*<       go to 200 >*/
    goto L200;
/*<  100  l = k1 >*/
L100:
    l = k1;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*<       k2 = k1-nu >*/
    k2 = k1 - *nu;
/*  main loop for the different points. */
/*<       do 180 i=1,m >*/
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
/*<  140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150 >*/
L140:
	if (arg < t[l1] || l == nk1) {
	    goto L150;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 140 >*/
	goto L140;
/*  evaluate the non-zero b-splines of degree k-nu at arg. */
/*<  150    call fpbspl(t,n,kk,arg,l,h) >*/
L150:
	fpbspl_(&t[1], n, &kk, &arg, &l, h__);
/*  find the value of the derivative at x=arg. */
/*<         sp = 0. >*/
	sp = 0.;
/*<         ll = l-k1 >*/
	ll = l - k1;
/*<         do 160 j=1,k2 >*/
	i__2 = k2;
	for (j = 1; j <= i__2; ++j) {
/*<           ll = ll+1 >*/
	    ++ll;
/*<           sp = sp+wrk(ll)*h(j) >*/
	    sp += wrk[ll] * h__[j - 1];
/*<  160    continue >*/
/* L160: */
	}
/*<         y(i) = sp >*/
	y[i__] = sp;
/*<  180  continue >*/
/* L180: */
    }
/*<  200  return >*/
L200:
    return 0;
/*<       end >*/
} /* splder_ */

