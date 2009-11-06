/* fpdisc.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpdisc(t,n,k2,b,nest) >*/
/* Subroutine */ int fpdisc_(doublereal *t, integer *n, integer *k2, 
	doublereal *b, integer *nest)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[12];
    static integer i__, j, k, l, k1;
    static doublereal an;
    static integer ik, jk, lj, lk, lp, nk1;
    static doublereal fac;
    static integer lmk;
    static doublereal prod;
    static integer nrint;

/*  subroutine fpdisc calculates the discontinuity jumps of the kth */
/*  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1) */
/*  ..scalar arguments.. */
/*<       integer n,k2,nest >*/
/*  ..array arguments.. */
/*<       real t(n),b(nest,k2) >*/
/*  ..local scalars.. */
/*<       real an,fac,prod >*/
/*<       integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint >*/
/*  ..local array.. */
/*<       real h(12) >*/
/*  .. */
/*<       k1 = k2-1 >*/
    /* Parameter adjustments */
    --t;
    b_dim1 = *nest;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    k1 = *k2 - 1;
/*<       k = k1-1 >*/
    k = k1 - 1;
/*<       nk1 = n-k1 >*/
    nk1 = *n - k1;
/*<       nrint = nk1-k >*/
    nrint = nk1 - k;
/*<       an = nrint >*/
    an = (doublereal) nrint;
/*<       fac = an/(t(nk1+1)-t(k1)) >*/
    fac = an / (t[nk1 + 1] - t[k1]);
/*<       do 40 l=k2,nk1 >*/
    i__1 = nk1;
    for (l = *k2; l <= i__1; ++l) {
/*<         lmk = l-k1 >*/
	lmk = l - k1;
/*<         do 10 j=1,k1 >*/
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
/*<           ik = j+k1 >*/
	    ik = j + k1;
/*<           lj = l+j >*/
	    lj = l + j;
/*<           lk = lj-k2 >*/
	    lk = lj - *k2;
/*<           h(j) = t(l)-t(lk) >*/
	    h__[j - 1] = t[l] - t[lk];
/*<           h(ik) = t(l)-t(lj) >*/
	    h__[ik - 1] = t[l] - t[lj];
/*<   10    continue >*/
/* L10: */
	}
/*<         lp = lmk >*/
	lp = lmk;
/*<         do 30 j=1,k2 >*/
	i__2 = *k2;
	for (j = 1; j <= i__2; ++j) {
/*<           jk = j >*/
	    jk = j;
/*<           prod = h(j) >*/
	    prod = h__[j - 1];
/*<           do 20 i=1,k >*/
	    i__3 = k;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             jk = jk+1 >*/
		++jk;
/*<             prod = prod*h(jk)*fac >*/
		prod = prod * h__[jk - 1] * fac;
/*<   20      continue >*/
/* L20: */
	    }
/*<           lk = lp+k1 >*/
	    lk = lp + k1;
/*<           b(lmk,j) = (t(lk)-t(lp))/prod >*/
	    b[lmk + j * b_dim1] = (t[lk] - t[lp]) / prod;
/*<           lp = lp+1 >*/
	    ++lp;
/*<   30    continue >*/
/* L30: */
	}
/*<   40  continue >*/
/* L40: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpdisc_ */

