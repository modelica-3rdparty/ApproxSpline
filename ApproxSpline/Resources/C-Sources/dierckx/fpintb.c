/* fpintb.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpintb(t,n,bint,nk1,x,y) >*/
/* Subroutine */ int fpintb_(doublereal *t, integer *n, doublereal *bint, 
	integer *nk1, doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal a, b, f, h__[6];
    static integer i__, j, k, l;
    static doublereal h1[6];
    static integer j1, k1, l0, ia, ib;
    static doublereal ak;
    static integer li, lj, lk, it;
    static doublereal arg, one;
    static integer min__;
    static doublereal aint[6];

/*  subroutine fpintb calculates integrals of the normalized b-splines */
/*  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n. */
/*  it makes use of the formulae of gaffney for the calculation of */
/*  indefinite integrals of b-splines. */

/*  calling sequence: */
/*     call fpintb(t,n,bint,nk1,x,y) */

/*  input parameters: */
/*    t    : real array,length n, containing the position of the knots. */
/*    n    : integer value, giving the number of knots. */
/*    nk1  : integer value, giving the number of b-splines of degree k, */
/*           defined on the set of knots ,i.e. nk1 = n-k-1. */
/*    x,y  : real values, containing the end points of the integration */
/*           interval. */
/*  output parameter: */
/*    bint : array,length nk1, containing the integrals of the b-splines. */
/*  .. */
/*  ..scalars arguments.. */
/*<       integer n,nk1 >*/
/*<       real x,y >*/
/*  ..array arguments.. */
/*<       real t(n),bint(nk1) >*/
/*  ..local scalars.. */
/*<       integer i,ia,ib,it,j,j1,k,k1,l,li,lj,lk,l0,min >*/
/*<       real a,ak,arg,b,f,one >*/
/*  ..local arrays.. */
/*<       real aint(6),h(6),h1(6) >*/
/*  initialization. */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --t;
    --bint;

    /* Function Body */
    one = 1.;
/*<       k1 = n-nk1 >*/
    k1 = *n - *nk1;
/*<       ak = k1 >*/
    ak = (doublereal) k1;
/*<       k = k1-1 >*/
    k = k1 - 1;
/*<       do 10 i=1,nk1 >*/
    i__1 = *nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         bint(i) = 0. >*/
	bint[i__] = 0.;
/*<   10  continue >*/
/* L10: */
    }
/*  the integration limits are arranged in increasing order. */
/*<       a = x >*/
    a = *x;
/*<       b = y >*/
    b = *y;
/*<       min = 0 >*/
    min__ = 0;
/*<       if(a-b) 30,160,20 >*/
    if ((d__1 = a - b) < 0.) {
	goto L30;
    } else if (d__1 == 0) {
	goto L160;
    } else {
	goto L20;
    }
/*<   20  a = y >*/
L20:
    a = *y;
/*<       b = x >*/
    b = *x;
/*<       min = 1 >*/
    min__ = 1;
/*<   30  if(a.lt.t(k1)) a = t(k1) >*/
L30:
    if (a < t[k1]) {
	a = t[k1];
    }
/*<       if(b.gt.t(nk1+1)) b = t(nk1+1) >*/
    if (b > t[*nk1 + 1]) {
	b = t[*nk1 + 1];
    }
/*  using the expression of gaffney for the indefinite integral of a */
/*  b-spline we find that */
/*  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1) */
/*    where for t(l) <= x < t(l+1) */
/*    res(j,x) = 0, j=1,2,...,l-k-1 */
/*             = 1, j=l+1,l+2,...,nk1 */
/*             = aint(j+k-l+1), j=l-k,l-k+1,...,l */
/*               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i))) */
/*                 i=0,1,...,k */
/*<       l = k1 >*/
    l = k1;
/*<       l0 = l+1 >*/
    l0 = l + 1;
/*  set arg = a. */
/*<       arg = a >*/
    arg = a;
/*<       do 90 it=1,2 >*/
    for (it = 1; it <= 2; ++it) {
/*  search for the knot interval t(l) <= arg < t(l+1). */
/*<   40    if(arg.lt.t(l0) .or. l.eq.nk1) go to 50 >*/
L40:
	if (arg < t[l0] || l == *nk1) {
	    goto L50;
	}
/*<         l = l0 >*/
	l = l0;
/*<         l0 = l+1 >*/
	l0 = l + 1;
/*<         go to 40 >*/
	goto L40;
/*  calculation of aint(j), j=1,2,...,k+1. */
/*  initialization. */
/*<   50    do 55 j=1,k1 >*/
L50:
	i__1 = k1;
	for (j = 1; j <= i__1; ++j) {
/*<           aint(j) = 0. >*/
	    aint[j - 1] = 0.;
/*<   55    continue >*/
/* L55: */
	}
/*<         aint(1) = (arg-t(l))/(t(l+1)-t(l)) >*/
	aint[0] = (arg - t[l]) / (t[l + 1] - t[l]);
/*<         h1(1) = one >*/
	h1[0] = one;
/*<         do 70 j=1,k >*/
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
/*  evaluation of the non-zero b-splines of degree j at arg,i.e. */
/*    h(i+1) = nl-j+i,j(arg), i=0,1,...,j. */
/*<           h(1) = 0. >*/
	    h__[0] = 0.;
/*<           do 60 i=1,j >*/
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             li = l+i >*/
		li = l + i__;
/*<             lj = li-j >*/
		lj = li - j;
/*<             f = h1(i)/(t(li)-t(lj)) >*/
		f = h1[i__ - 1] / (t[li] - t[lj]);
/*<             h(i) = h(i)+f*(t(li)-arg) >*/
		h__[i__ - 1] += f * (t[li] - arg);
/*<             h(i+1) = f*(arg-t(lj)) >*/
		h__[i__] = f * (arg - t[lj]);
/*<   60      continue >*/
/* L60: */
	    }
/*  updating of the integrals aint. */
/*<           j1 = j+1 >*/
	    j1 = j + 1;
/*<           do 70 i=1,j1 >*/
	    i__2 = j1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             li = l+i >*/
		li = l + i__;
/*<             lj = li-j1 >*/
		lj = li - j1;
/*<             aint(i) = aint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj)) >*/
		aint[i__ - 1] += h__[i__ - 1] * (arg - t[lj]) / (t[li] - t[lj]
			);
/*<             h1(i) = h(i) >*/
		h1[i__ - 1] = h__[i__ - 1];
/*<   70    continue >*/
/* L70: */
	    }
	}
/*<         if(it.eq.2) go to 100 >*/
	if (it == 2) {
	    goto L100;
	}
/*  updating of the integrals bint */
/*<         lk = l-k >*/
	lk = l - k;
/*<         ia = lk >*/
	ia = lk;
/*<         do 80 i=1,k1 >*/
	i__2 = k1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           bint(lk) = -aint(i) >*/
	    bint[lk] = -aint[i__ - 1];
/*<           lk = lk+1 >*/
	    ++lk;
/*<   80    continue >*/
/* L80: */
	}
/*  set arg = b. */
/*<         arg = b >*/
	arg = b;
/*<   90  continue >*/
/* L90: */
    }
/*  updating of the integrals bint. */
/*<  100  lk = l-k >*/
L100:
    lk = l - k;
/*<       ib = lk-1 >*/
    ib = lk - 1;
/*<       do 110 i=1,k1 >*/
    i__2 = k1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         bint(lk) = bint(lk)+aint(i) >*/
	bint[lk] += aint[i__ - 1];
/*<         lk = lk+1 >*/
	++lk;
/*<  110  continue >*/
/* L110: */
    }
/*<       if(ib.lt.ia) go to 130 >*/
    if (ib < ia) {
	goto L130;
    }
/*<       do 120 i=ia,ib >*/
    i__2 = ib;
    for (i__ = ia; i__ <= i__2; ++i__) {
/*<         bint(i) = bint(i)+one >*/
	bint[i__] += one;
/*<  120  continue >*/
/* L120: */
    }
/*  the scaling factors are taken into account. */
/*<  130  f = one/ak >*/
L130:
    f = one / ak;
/*<       do 140 i=1,nk1 >*/
    i__2 = *nk1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         j = i+k1 >*/
	j = i__ + k1;
/*<         bint(i) = bint(i)*(t(j)-t(i))*f >*/
	bint[i__] = bint[i__] * (t[j] - t[i__]) * f;
/*<  140  continue >*/
/* L140: */
    }
/*  the order of the integration limits is taken into account. */
/*<       if(min.eq.0) go to 160 >*/
    if (min__ == 0) {
	goto L160;
    }
/*<       do 150 i=1,nk1 >*/
    i__2 = *nk1;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         bint(i) = -bint(i) >*/
	bint[i__] = -bint[i__];
/*<  150  continue >*/
/* L150: */
    }
/*<  160  return >*/
L160:
    return 0;
/*<       end >*/
} /* fpintb_ */

