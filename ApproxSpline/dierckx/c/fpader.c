/* fpader.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpader(t,n,c,k1,x,l,d) >*/
/* Subroutine */ int fpader_(doublereal *t, integer *n, doublereal *c__, 
	integer *k1, doublereal *x, integer *l, doublereal *d__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, j1, j2;
    static doublereal ak;
    static integer ik, jj, ki, kj, li, lj, lk;
    static doublereal fac, one;

/*  subroutine fpader calculates the derivatives */
/*             (j-1) */
/*     d(j) = s     (x) , j=1,2,...,k1 */
/*  of a spline of order k1 at the point t(l)<=x<t(l+1), using the */
/*  stable recurrence scheme of de boor */
/*  .. */
/*  ..scalar arguments.. */
/*<       real x >*/
/*<       integer n,k1,l >*/
/*  ..array arguments.. */
/*<       real t(n),c(n),d(k1) >*/
/*  ..local scalars.. */
/*<       integer i,ik,j,jj,j1,j2,ki,kj,li,lj,lk >*/
/*<       real ak,fac,one >*/
/*  ..local array.. */
/*<       real h(6) >*/
/*  .. */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --c__;
    --t;
    --d__;

    /* Function Body */
    one = 1.;
/*<       lk = l-k1 >*/
    lk = *l - *k1;
/*<       do 100 i=1,k1 >*/
    i__1 = *k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         ik = i+lk >*/
	ik = i__ + lk;
/*<         h(i) = c(ik) >*/
	h__[i__ - 1] = c__[ik];
/*<  100  continue >*/
/* L100: */
    }
/*<       kj = k1 >*/
    kj = *k1;
/*<       fac = one >*/
    fac = one;
/*<       do 700 j=1,k1 >*/
    i__1 = *k1;
    for (j = 1; j <= i__1; ++j) {
/*<         ki = kj >*/
	ki = kj;
/*<         j1 = j+1 >*/
	j1 = j + 1;
/*<         if(j.eq.1) go to 300 >*/
	if (j == 1) {
	    goto L300;
	}
/*<         i = k1 >*/
	i__ = *k1;
/*<         do 200 jj=j,k1 >*/
	i__2 = *k1;
	for (jj = j; jj <= i__2; ++jj) {
/*<           li = i+lk >*/
	    li = i__ + lk;
/*<           lj = li+kj >*/
	    lj = li + kj;
/*<           h(i) = (h(i)-h(i-1))/(t(lj)-t(li)) >*/
	    h__[i__ - 1] = (h__[i__ - 1] - h__[i__ - 2]) / (t[lj] - t[li]);
/*<           i = i-1 >*/
	    --i__;
/*<  200    continue >*/
/* L200: */
	}
/*<  300    do 400 i=j,k1 >*/
L300:
	i__2 = *k1;
	for (i__ = j; i__ <= i__2; ++i__) {
/*<           d(i) = h(i) >*/
	    d__[i__] = h__[i__ - 1];
/*<  400    continue >*/
/* L400: */
	}
/*<         if(j.eq.k1) go to 600 >*/
	if (j == *k1) {
	    goto L600;
	}
/*<         do 500 jj=j1,k1 >*/
	i__2 = *k1;
	for (jj = j1; jj <= i__2; ++jj) {
/*<           ki = ki-1 >*/
	    --ki;
/*<           i = k1 >*/
	    i__ = *k1;
/*<           do 500 j2=jj,k1 >*/
	    i__3 = *k1;
	    for (j2 = jj; j2 <= i__3; ++j2) {
/*<             li = i+lk >*/
		li = i__ + lk;
/*<             lj = li+ki >*/
		lj = li + ki;
/*<             d(i) = ((x-t(li))*d(i)+(t(lj)-x)*d(i-1))/(t(lj)-t(li)) >*/
		d__[i__] = ((*x - t[li]) * d__[i__] + (t[lj] - *x) * d__[i__ 
			- 1]) / (t[lj] - t[li]);
/*<             i = i-1 >*/
		--i__;
/*<  500    continue >*/
/* L500: */
	    }
	}
/*<  600    d(j) = d(k1)*fac >*/
L600:
	d__[j] = d__[*k1] * fac;
/*<         ak = k1-j >*/
	ak = (doublereal) (*k1 - j);
/*<         fac = fac*ak >*/
	fac *= ak;
/*<         kj = kj-1 >*/
	--kj;
/*<  700  continue >*/
/* L700: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpader_ */

