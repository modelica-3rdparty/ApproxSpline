/* fpadpo.f -- translated by f2c (version 20061008).
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

/* Table of constant values */

static integer c__0 = 0;

/*<       subroutine fpadpo(idim,t,n,c,nc,k,cp,np,cc,t1,t2) >*/
/* Subroutine */ int fpadpo_(integer *idim, doublereal *t, integer *n, 
	doublereal *c__, integer *nc, integer *k, doublereal *cp, integer *np,
	 doublereal *cc, doublereal *t1, doublereal *t2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, k1, l1, n1, n2, ii, jj, nk1, nk2;
    extern /* Subroutine */ int fpinst_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *);

/*  given a idim-dimensional spline curve of degree k, in its b-spline */
/*  representation ( knots t(j),j=1,...,n , b-spline coefficients c(j), */
/*  j=1,...,nc) and given also a polynomial curve in its b-spline */
/*  representation ( coefficients cp(j), j=1,...,np), subroutine fpadpo */
/*  calculates the b-spline representation (coefficients c(j),j=1,...,nc) */
/*  of the sum of the two curves. */

/*  other subroutine required : fpinst */

/*  .. */
/*  ..scalar arguments.. */
/*<       integer idim,k,n,nc,np >*/
/*  ..array arguments.. */
/*<       real t(n),c(nc),cp(np),cc(nc),t1(n),t2(n) >*/
/*  ..local scalars.. */
/*<       integer i,ii,j,jj,k1,l,l1,n1,n2,nk1,nk2 >*/
/*  .. */
/*<       k1 = k+1 >*/
    /* Parameter adjustments */
    --t2;
    --t1;
    --t;
    --cc;
    --c__;
    --cp;

    /* Function Body */
    k1 = *k + 1;
/*<       nk1 = n-k1 >*/
    nk1 = *n - k1;
/*  initialization */
/*<       j = 1 >*/
    j = 1;
/*<       l = 1 >*/
    l = 1;
/*<       do 20 jj=1,idim >*/
    i__1 = *idim;
    for (jj = 1; jj <= i__1; ++jj) {
/*<         l1 = j >*/
	l1 = j;
/*<         do 10 ii=1,k1 >*/
	i__2 = k1;
	for (ii = 1; ii <= i__2; ++ii) {
/*<           cc(l1) = cp(l) >*/
	    cc[l1] = cp[l];
/*<           l1 = l1+1 >*/
	    ++l1;
/*<           l = l+1 >*/
	    ++l;
/*<   10    continue >*/
/* L10: */
	}
/*<         j = j+n >*/
	j += *n;
/*<         l = l+k1 >*/
	l += k1;
/*<   20  continue >*/
/* L20: */
    }
/*<       if(nk1.eq.k1) go to 70 >*/
    if (nk1 == k1) {
	goto L70;
    }
/*<       n1 = k1*2 >*/
    n1 = k1 << 1;
/*<       j = n >*/
    j = *n;
/*<       l = n1 >*/
    l = n1;
/*<       do 30 i=1,k1 >*/
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         t1(i) = t(i) >*/
	t1[i__] = t[i__];
/*<         t1(l) = t(j) >*/
	t1[l] = t[j];
/*<         l = l-1 >*/
	--l;
/*<         j = j-1 >*/
	--j;
/*<   30  continue >*/
/* L30: */
    }
/*  find the b-spline representation of the given polynomial curve */
/*  according to the given set of knots. */
/*<       nk2 = nk1-1 >*/
    nk2 = nk1 - 1;
/*<       do 60 l=k1,nk2 >*/
    i__1 = nk2;
    for (l = k1; l <= i__1; ++l) {
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         j = 1 >*/
	j = 1;
/*<         do 40 i=1,idim >*/
	i__2 = *idim;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           call fpinst(0,t1,n1,cc(j),k,t(l1),l,t2,n2,cc(j),n) >*/
	    fpinst_(&c__0, &t1[1], &n1, &cc[j], k, &t[l1], &l, &t2[1], &n2, &
		    cc[j], n);
/*<           j = j+n >*/
	    j += *n;
/*<   40    continue >*/
/* L40: */
	}
/*<         do 50 i=1,n2 >*/
	i__2 = n2;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           t1(i) = t2(i) >*/
	    t1[i__] = t2[i__];
/*<   50    continue >*/
/* L50: */
	}
/*<         n1 = n2 >*/
	n1 = n2;
/*<   60  continue >*/
/* L60: */
    }
/*  find the b-spline representation of the resulting curve. */
/*<   70  j = 1 >*/
L70:
    j = 1;
/*<       do 90 jj=1,idim >*/
    i__1 = *idim;
    for (jj = 1; jj <= i__1; ++jj) {
/*<         l = j >*/
	l = j;
/*<         do 80 i=1,nk1 >*/
	i__2 = nk1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           c(l) = cc(l)+c(l) >*/
	    c__[l] = cc[l] + c__[l];
/*<           l = l+1 >*/
	    ++l;
/*<   80    continue >*/
/* L80: */
	}
/*<         j = j+n >*/
	j += *n;
/*<   90  continue >*/
/* L90: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpadpo_ */

