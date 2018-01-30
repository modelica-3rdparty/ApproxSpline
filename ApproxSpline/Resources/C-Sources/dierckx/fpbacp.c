/* fpbacp.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpbacp(a,b,z,n,k,c,k1,nest) >*/
/* Subroutine */ int fpbacp_(doublereal *a, doublereal *b, doublereal *z__, 
	integer *n, integer *k, doublereal *c__, integer *k1, integer *nest)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, i1, l0, l1, n2;
    static doublereal store;

/*  subroutine fpbacp calculates the solution of the system of equations */
/*  g * c = z  with g  a n x n upper triangular matrix of the form */
/*            ! a '   ! */
/*        g = !   ' b ! */
/*            ! 0 '   ! */
/*  with b a n x k matrix and a a (n-k) x (n-k) upper triangular */
/*  matrix of bandwidth k1. */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer n,k,k1,nest >*/
/*  ..array arguments.. */
/*<       real a(nest,k1),b(nest,k),z(n),c(n) >*/
/*  ..local scalars.. */
/*<       integer i,i1,j,l,l0,l1,n2 >*/
/*<       real store >*/
/*  .. */
/*<       n2 = n-k >*/
    /* Parameter adjustments */
    --c__;
    --z__;
    b_dim1 = *nest;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *nest;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    n2 = *n - *k;
/*<       l = n >*/
    l = *n;
/*<       do 30 i=1,k >*/
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         store = z(l) >*/
	store = z__[l];
/*<         j = k+2-i >*/
	j = *k + 2 - i__;
/*<         if(i.eq.1) go to 20 >*/
	if (i__ == 1) {
	    goto L20;
	}
/*<         l0 = l >*/
	l0 = l;
/*<         do 10 l1=j,k >*/
	i__2 = *k;
	for (l1 = j; l1 <= i__2; ++l1) {
/*<           l0 = l0+1 >*/
	    ++l0;
/*<           store = store-c(l0)*b(l,l1) >*/
	    store -= c__[l0] * b[l + l1 * b_dim1];
/*<   10    continue >*/
/* L10: */
	}
/*<   20    c(l) = store/b(l,j-1) >*/
L20:
	c__[l] = store / b[l + (j - 1) * b_dim1];
/*<         l = l-1 >*/
	--l;
/*<         if(l.eq.0) go to 80 >*/
	if (l == 0) {
	    goto L80;
	}
/*<   30  continue >*/
/* L30: */
    }
/*<       do 50 i=1,n2 >*/
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         store = z(i) >*/
	store = z__[i__];
/*<         l = n2 >*/
	l = n2;
/*<         do 40 j=1,k >*/
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
/*<           l = l+1 >*/
	    ++l;
/*<           store = store-c(l)*b(i,j) >*/
	    store -= c__[l] * b[i__ + j * b_dim1];
/*<   40    continue >*/
/* L40: */
	}
/*<         c(i) = store >*/
	c__[i__] = store;
/*<   50  continue >*/
/* L50: */
    }
/*<       i = n2 >*/
    i__ = n2;
/*<       c(i) = c(i)/a(i,1) >*/
    c__[i__] /= a[i__ + a_dim1];
/*<       if(i.eq.1) go to 80 >*/
    if (i__ == 1) {
	goto L80;
    }
/*<       do 70 j=2,n2 >*/
    i__1 = n2;
    for (j = 2; j <= i__1; ++j) {
/*<         i = i-1 >*/
	--i__;
/*<         store = c(i) >*/
	store = c__[i__];
/*<         i1 = k >*/
	i1 = *k;
/*<         if(j.le.k) i1=j-1 >*/
	if (j <= *k) {
	    i1 = j - 1;
	}
/*<         l = i >*/
	l = i__;
/*<         do 60 l0=1,i1 >*/
	i__2 = i1;
	for (l0 = 1; l0 <= i__2; ++l0) {
/*<           l = l+1 >*/
	    ++l;
/*<           store = store-c(l)*a(i,l0+1) >*/
	    store -= c__[l] * a[i__ + (l0 + 1) * a_dim1];
/*<   60    continue >*/
/* L60: */
	}
/*<         c(i) = store/a(i,1) >*/
	c__[i__] = store / a[i__ + a_dim1];
/*<   70  continue >*/
/* L70: */
    }
/*<   80  return >*/
L80:
    return 0;
/*<       end >*/
} /* fpbacp_ */

