/* fpbspl.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpbspl(t,n,k,x,l,h) >*/
/* Subroutine */ int fpbspl_(doublereal *t, integer *n, integer *k, 
	doublereal *x, integer *l, doublereal *h__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal f;
    static integer i__, j;
    static doublereal hh[5];
    static integer li, lj;
    static doublereal one;

/*  subroutine fpbspl evaluates the (k+1) non-zero b-splines of */
/*  degree k at t(l) <= x < t(l+1) using the stable recurrence */
/*  relation of de boor and cox. */
/*  .. */
/*  ..scalar arguments.. */
/*<       real x >*/
/*<       integer n,k,l >*/
/*  ..array arguments.. */
/*<       real t(n),h(6) >*/
/*  ..local scalars.. */
/*<       real f,one >*/
/*<       integer i,j,li,lj >*/
/*  ..local arrays.. */
/*<       real hh(5) >*/
/*  .. */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --t;
    --h__;

    /* Function Body */
    one = 1.;
/*<       h(1) = one >*/
    h__[1] = one;
/*<       do 20 j=1,k >*/
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
/*<         do 10 i=1,j >*/
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           hh(i) = h(i) >*/
	    hh[i__ - 1] = h__[i__];
/*<   10    continue >*/
/* L10: */
	}
/*<         h(1) = 0. >*/
	h__[1] = 0.;
/*<         do 20 i=1,j >*/
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           li = l+i >*/
	    li = *l + i__;
/*<           lj = li-j >*/
	    lj = li - j;
/*<           f = hh(i)/(t(li)-t(lj)) >*/
	    f = hh[i__ - 1] / (t[li] - t[lj]);
/*<           h(i) = h(i)+f*(t(li)-x) >*/
	    h__[i__] += f * (t[li] - *x);
/*<           h(i+1) = f*(x-t(lj)) >*/
	    h__[i__ + 1] = f * (*x - t[lj]);
/*<   20  continue >*/
/* L20: */
	}
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpbspl_ */

