/* fpinst.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest) >*/
/* Subroutine */ int fpinst_(integer *iopt, doublereal *t, integer *n, 
	doublereal *c__, integer *k, doublereal *x, integer *l, doublereal *
	tt, integer *nn, doublereal *cc, integer *nest)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, m, i1, k1, mk, nk, nl, ll, nk1;
    static doublereal fac, one, per;

/*  given the b-spline representation (knots t(j),j=1,2,...,n, b-spline */
/*  coefficients c(j),j=1,2,...,n-k-1) of a spline of degree k, fpinst */
/*  calculates the b-spline representation (knots tt(j),j=1,2,...,nn, */
/*  b-spline coefficients cc(j),j=1,2,...,nn-k-1) of the same spline if */
/*  an additional knot is inserted at the point x situated in the inter- */
/*  val t(l)<=x<t(l+1). iopt denotes whether (iopt.ne.0) or not (iopt=0) */
/*  the given spline is periodic. in case of a periodic spline at least */
/*  one of the following conditions must be fulfilled: l>2*k or l<n-2*k. */

/*  ..scalar arguments.. */
/*<       integer k,n,l,nn,iopt,nest >*/
/*<       real x >*/
/*  ..array arguments.. */
/*<       real t(nest),c(nest),tt(nest),cc(nest) >*/
/*  ..local scalars.. */
/*<       real fac,per,one >*/
/*<       integer i,i1,j,k1,m,mk,nk,nk1,nl,ll >*/
/*  .. */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --cc;
    --tt;
    --c__;
    --t;

    /* Function Body */
    one = 1.;
/*<       k1 = k+1 >*/
    k1 = *k + 1;
/*<       nk1 = n-k1 >*/
    nk1 = *n - k1;
/*  the new knots */
/*<       ll = l+1 >*/
    ll = *l + 1;
/*<       i = n >*/
    i__ = *n;
/*<       do 10 j=ll,n >*/
    i__1 = *n;
    for (j = ll; j <= i__1; ++j) {
/*<          tt(i+1) = t(i) >*/
	tt[i__ + 1] = t[i__];
/*<          i = i-1 >*/
	--i__;
/*<   10  continue >*/
/* L10: */
    }
/*<       tt(ll) = x >*/
    tt[ll] = *x;
/*<       do 20 j=1,l >*/
    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
/*<          tt(j) = t(j) >*/
	tt[j] = t[j];
/*<   20  continue >*/
/* L20: */
    }
/*  the new b-spline coefficients */
/*<       i = nk1 >*/
    i__ = nk1;
/*<       do 30 j=l,nk1 >*/
    i__1 = nk1;
    for (j = *l; j <= i__1; ++j) {
/*<          cc(i+1) = c(i) >*/
	cc[i__ + 1] = c__[i__];
/*<          i = i-1 >*/
	--i__;
/*<   30  continue >*/
/* L30: */
    }
/*<       i = l >*/
    i__ = *l;
/*<       do 40 j=1,k >*/
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
/*<          m = i+k1 >*/
	m = i__ + k1;
/*<          fac = (x-tt(i))/(tt(m)-tt(i)) >*/
	fac = (*x - tt[i__]) / (tt[m] - tt[i__]);
/*<          i1 = i-1 >*/
	i1 = i__ - 1;
/*<          cc(i) = fac*c(i)+(one-fac)*c(i1) >*/
	cc[i__] = fac * c__[i__] + (one - fac) * c__[i1];
/*<          i = i1 >*/
	i__ = i1;
/*<   40  continue >*/
/* L40: */
    }
/*<       do 50 j=1,i >*/
    i__1 = i__;
    for (j = 1; j <= i__1; ++j) {
/*<          cc(j) = c(j) >*/
	cc[j] = c__[j];
/*<   50  continue >*/
/* L50: */
    }
/*<       nn = n+1 >*/
    *nn = *n + 1;
/*<       if(iopt.eq.0) return >*/
    if (*iopt == 0) {
	return 0;
    }
/*   incorporate the boundary conditions for a periodic spline. */
/*<       nk = nn-k >*/
    nk = *nn - *k;
/*<       nl = nk-k1 >*/
    nl = nk - k1;
/*<       per = tt(nk)-tt(k1) >*/
    per = tt[nk] - tt[k1];
/*<       i = k1 >*/
    i__ = k1;
/*<       j = nk >*/
    j = nk;
/*<       if(ll.le.nl) go to 70 >*/
    if (ll <= nl) {
	goto L70;
    }
/*<       do 60 m=1,k >*/
    i__1 = *k;
    for (m = 1; m <= i__1; ++m) {
/*<          mk = m+nl >*/
	mk = m + nl;
/*<          cc(m) = cc(mk) >*/
	cc[m] = cc[mk];
/*<          i = i-1 >*/
	--i__;
/*<          j = j-1 >*/
	--j;
/*<          tt(i) = tt(j)-per >*/
	tt[i__] = tt[j] - per;
/*<   60  continue >*/
/* L60: */
    }
/*<       return >*/
    return 0;
/*<   70  if(ll.gt.(k1+k)) return >*/
L70:
    if (ll > k1 + *k) {
	return 0;
    }
/*<       do 80 m=1,k >*/
    i__1 = *k;
    for (m = 1; m <= i__1; ++m) {
/*<          mk = m+nl >*/
	mk = m + nl;
/*<          cc(mk) = cc(m) >*/
	cc[mk] = cc[m];
/*<          i = i+1 >*/
	++i__;
/*<          j = j+1 >*/
	++j;
/*<          tt(j) = tt(i)+per >*/
	tt[j] = tt[i__] + per;
/*<   80  continue >*/
/* L80: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpinst_ */

