/* fpback.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpback(a,z,n,k,c,nest) >*/
/* Subroutine */ int fpback_(doublereal *a, doublereal *z__, integer *n, 
	integer *k, doublereal *c__, integer *nest)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, m, i1, k1;
    static doublereal store;

/*  subroutine fpback calculates the solution of the system of */
/*  equations a*c = z with a a n x n upper triangular matrix */
/*  of bandwidth k. */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer n,k,nest >*/
/*  ..array arguments.. */
/*<       real a(nest,k),z(n),c(n) >*/
/*  ..local scalars.. */
/*<       real store >*/
/*<       integer i,i1,j,k1,l,m >*/
/*  .. */
/*<       k1 = k-1 >*/
    /* Parameter adjustments */
    --c__;
    --z__;
    a_dim1 = *nest;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    k1 = *k - 1;
/*<       c(n) = z(n)/a(n,1) >*/
    c__[*n] = z__[*n] / a[*n + a_dim1];
/*<       i = n-1 >*/
    i__ = *n - 1;
/*<       if(i.eq.0) go to 30 >*/
    if (i__ == 0) {
	goto L30;
    }
/*<       do 20 j=2,n >*/
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/*<         store = z(i) >*/
	store = z__[i__];
/*<         i1 = k1 >*/
	i1 = k1;
/*<         if(j.le.k1) i1 = j-1 >*/
	if (j <= k1) {
	    i1 = j - 1;
	}
/*<         m = i >*/
	m = i__;
/*<         do 10 l=1,i1 >*/
	i__2 = i1;
	for (l = 1; l <= i__2; ++l) {
/*<           m = m+1 >*/
	    ++m;
/*<           store = store-c(m)*a(i,l+1) >*/
	    store -= c__[m] * a[i__ + (l + 1) * a_dim1];
/*<   10    continue >*/
/* L10: */
	}
/*<         c(i) = store/a(i,1) >*/
	c__[i__] = store / a[i__ + a_dim1];
/*<         i = i-1 >*/
	--i__;
/*<   20  continue >*/
/* L20: */
    }
/*<   30  return >*/
L30:
    return 0;
/*<       end >*/
} /* fpback_ */

