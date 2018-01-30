/* sproot.f -- translated by f2c (version 20061008).
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

/*<       subroutine sproot(t,n,c,zero,mest,m,ier) >*/
/* Subroutine */ int sproot_(doublereal *t, integer *n, doublereal *c__, 
	doublereal *zero, integer *mest, integer *m, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal y[3], a0, a1, a2, a3, b0, b1, c1, c2, c3, c4;
    static integer j1;
    static doublereal c5, d4, d5, h1, h2;
    static integer n4;
    static doublereal t1, t2, t3, t4, t5;
    static logical z0, z1, z2, z3, z4;
    static doublereal ah, bh, zz;
    static logical nz0, nz1, nz2, nz3, nz4;
    static doublereal two, three;
    extern /* Subroutine */ int fpcuro_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);

/*  subroutine sproot finds the zeros of a cubic spline s(x),which is */
/*  given in its normalized b-spline representation. */

/*  calling sequence: */
/*     call sproot(t,n,c,zero,mest,m,ier) */

/*  input parameters: */
/*    t    : real array,length n, containing the knots of s(x). */
/*    n    : integer, containing the number of knots.  n>=8 */
/*    c    : real array,length n, containing the b-spline coefficients. */
/*    mest : integer, specifying the dimension of array zero. */

/*  output parameters: */
/*    zero : real array,lenth mest, containing the zeros of s(x). */
/*    m    : integer,giving the number of zeros. */
/*    ier  : error flag: */
/*      ier = 0: normal return. */
/*      ier = 1: the number of zeros exceeds mest. */
/*      ier =10: invalid input data (see restrictions). */

/*  other subroutines required: fpcuro */

/*  restrictions: */
/*    1) n>= 8. */
/*    2) t(4) < t(5) < ... < t(n-4) < t(n-3). */
/*       t(1) <= t(2) <= t(3) <= t(4) */
/*       t(n-3) <= t(n-2) <= t(n-1) <= t(n) */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/* .. */
/* ..scalar arguments.. */
/*<       integer n,mest,m,ier >*/
/*  ..array arguments.. */
/*<       real t(n),c(n),zero(mest) >*/
/*  ..local scalars.. */
/*<       integer i,j,j1,l,n4 >*/
/*<    >*/
/*<       logical z0,z1,z2,z3,z4,nz0,nz1,nz2,nz3,nz4 >*/
/*  ..local array.. */
/*<       real y(3) >*/
/*  .. */
/*  set some constants */
/*<       two = 0.2e+01 >*/
    /* Parameter adjustments */
    --c__;
    --t;
    --zero;

    /* Function Body */
    two = 2.;
/*<       three = 0.3e+01 >*/
    three = 3.;
/*  before starting computations a data check is made. if the input data */
/*  are invalid, control is immediately repassed to the calling program. */
/*<       n4 = n-4 >*/
    n4 = *n - 4;
/*<       ier = 10 >*/
    *ier = 10;
/*<       if(n.lt.8) go to 800 >*/
    if (*n < 8) {
	goto L800;
    }
/*<       j = n >*/
    j = *n;
/*<       do 10 i=1,3 >*/
    for (i__ = 1; i__ <= 3; ++i__) {
/*<         if(t(i).gt.t(i+1)) go to 800 >*/
	if (t[i__] > t[i__ + 1]) {
	    goto L800;
	}
/*<         if(t(j).lt.t(j-1)) go to 800 >*/
	if (t[j] < t[j - 1]) {
	    goto L800;
	}
/*<         j = j-1 >*/
	--j;
/*<   10  continue >*/
/* L10: */
    }
/*<       do 20 i=4,n4 >*/
    i__1 = n4;
    for (i__ = 4; i__ <= i__1; ++i__) {
/*<         if(t(i).ge.t(i+1)) go to 800 >*/
	if (t[i__] >= t[i__ + 1]) {
	    goto L800;
	}
/*<   20  continue >*/
/* L20: */
    }
/*  the problem considered reduces to finding the zeros of the cubic */
/*  polynomials pl(x) which define the cubic spline in each knot */
/*  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on */
/*  the condition that it belongs to the knot interval. */
/*  the cubic polynomial pl(x) is determined by computing s(t(l)), */
/*  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute */
/*  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of */
/*  splines and their derivatives, the value of s(t(l)) and s'(t(l)) */
/*  is already known from the foregoing knot interval. */
/*<       ier = 0 >*/
    *ier = 0;
/*  evaluate some constants for the first knot interval */
/*<       h1 = t(4)-t(3) >*/
    h1 = t[4] - t[3];
/*<       h2 = t(5)-t(4) >*/
    h2 = t[5] - t[4];
/*<       t1 = t(4)-t(2) >*/
    t1 = t[4] - t[2];
/*<       t2 = t(5)-t(3) >*/
    t2 = t[5] - t[3];
/*<       t3 = t(6)-t(4) >*/
    t3 = t[6] - t[4];
/*<       t4 = t(5)-t(2) >*/
    t4 = t[5] - t[2];
/*<       t5 = t(6)-t(3) >*/
    t5 = t[6] - t[3];
/*  calculate a0 = s(t(4)) and ah = s'(t(4)). */
/*<       c1 = c(1) >*/
    c1 = c__[1];
/*<       c2 = c(2) >*/
    c2 = c__[2];
/*<       c3 = c(3) >*/
    c3 = c__[3];
/*<       c4 = (c2-c1)/t4 >*/
    c4 = (c2 - c1) / t4;
/*<       c5 = (c3-c2)/t5 >*/
    c5 = (c3 - c2) / t5;
/*<       d4 = (h2*c1+t1*c2)/t4 >*/
    d4 = (h2 * c1 + t1 * c2) / t4;
/*<       d5 = (t3*c2+h1*c3)/t5 >*/
    d5 = (t3 * c2 + h1 * c3) / t5;
/*<       a0 = (h2*d4+h1*d5)/t2 >*/
    a0 = (h2 * d4 + h1 * d5) / t2;
/*<       ah = three*(h2*c4+h1*c5)/t2 >*/
    ah = three * (h2 * c4 + h1 * c5) / t2;
/*<       z1 = .true. >*/
    z1 = TRUE_;
/*<       if(ah.lt.0.) z1 = .false. >*/
    if (ah < 0.) {
	z1 = FALSE_;
    }
/*<       nz1 = .not.z1 >*/
    nz1 = ! z1;
/*<       m = 0 >*/
    *m = 0;
/*  main loop for the different knot intervals. */
/*<       do 300 l=4,n4 >*/
    i__1 = n4;
    for (l = 4; l <= i__1; ++l) {
/*  evaluate some constants for the knot interval t(l) <= x <= t(l+1). */
/*<         h1 = h2 >*/
	h1 = h2;
/*<         h2 = t(l+2)-t(l+1) >*/
	h2 = t[l + 2] - t[l + 1];
/*<         t1 = t2 >*/
	t1 = t2;
/*<         t2 = t3 >*/
	t2 = t3;
/*<         t3 = t(l+3)-t(l+1) >*/
	t3 = t[l + 3] - t[l + 1];
/*<         t4 = t5 >*/
	t4 = t5;
/*<         t5 = t(l+3)-t(l) >*/
	t5 = t[l + 3] - t[l];
/*  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)). */
/*<         c1 = c2 >*/
	c1 = c2;
/*<         c2 = c3 >*/
	c2 = c3;
/*<         c3 = c(l) >*/
	c3 = c__[l];
/*<         c4 = c5 >*/
	c4 = c5;
/*<         c5 = (c3-c2)/t5 >*/
	c5 = (c3 - c2) / t5;
/*<         d4 = (h2*c1+t1*c2)/t4 >*/
	d4 = (h2 * c1 + t1 * c2) / t4;
/*<         d5 = (h1*c3+t3*c2)/t5 >*/
	d5 = (h1 * c3 + t3 * c2) / t5;
/*<         b0 = (h2*d4+h1*d5)/t2 >*/
	b0 = (h2 * d4 + h1 * d5) / t2;
/*<         bh = three*(h2*c4+h1*c5)/t2 >*/
	bh = three * (h2 * c4 + h1 * c5) / t2;
/*  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial */
/*  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)). */
/*<         a1 = ah*h1 >*/
	a1 = ah * h1;
/*<         b1 = bh*h1 >*/
	b1 = bh * h1;
/*<         a2 = three*(b0-a0)-b1-two*a1 >*/
	a2 = three * (b0 - a0) - b1 - two * a1;
/*<         a3 = two*(a0-b0)+b1+a1 >*/
	a3 = two * (a0 - b0) + b1 + a1;
/*  test whether or not pl(x) could have a zero in the range */
/*  t(l) <= x <= t(l+1). */
/*<         z3 = .true. >*/
	z3 = TRUE_;
/*<         if(b1.lt.0.) z3 = .false. >*/
	if (b1 < 0.) {
	    z3 = FALSE_;
	}
/*<         nz3 = .not.z3 >*/
	nz3 = ! z3;
/*<         if(a0*b0.le.0.) go to 100 >*/
	if (a0 * b0 <= 0.) {
	    goto L100;
	}
/*<         z0 = .true. >*/
	z0 = TRUE_;
/*<         if(a0.lt.0.) z0 = .false. >*/
	if (a0 < 0.) {
	    z0 = FALSE_;
	}
/*<         nz0 = .not.z0 >*/
	nz0 = ! z0;
/*<         z2 = .true. >*/
	z2 = TRUE_;
/*<         if(a2.lt.0.) z2 = .false. >*/
	if (a2 < 0.) {
	    z2 = FALSE_;
	}
/*<         nz2 = .not.z2 >*/
	nz2 = ! z2;
/*<         z4 = .true. >*/
	z4 = TRUE_;
/*<         if(3.0*a3+a2.lt.0.) z4 = .false. >*/
	if (a3 * 3. + a2 < 0.) {
	    z4 = FALSE_;
	}
/*<         nz4 = .not.z4 >*/
	nz4 = ! z4;
/*<    >*/
	if (! (z0 && (nz1 && (z3 || z2 && nz4) || nz2 && z3 && z4) || nz0 && (
		z1 && (nz3 || nz2 && z4) || z2 && nz3 && nz4))) {
	    goto L200;
	}
/*  find the zeros of ql(y). */
/*<  100    call fpcuro(a3,a2,a1,a0,y,j) >*/
L100:
	fpcuro_(&a3, &a2, &a1, &a0, y, &j);
/*<         if(j.eq.0) go to 200 >*/
	if (j == 0) {
	    goto L200;
	}
/*  find which zeros of pl(x) are zeros of s(x). */
/*<         do 150 i=1,j >*/
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(y(i).lt.0. .or. y(i).gt.1.0) go to 150 >*/
	    if (y[i__ - 1] < 0. || y[i__ - 1] > 1.) {
		goto L150;
	    }
/*  test whether the number of zeros of s(x) exceeds mest. */
/*<           if(m.ge.mest) go to 700 >*/
	    if (*m >= *mest) {
		goto L700;
	    }
/*<           m = m+1 >*/
	    ++(*m);
/*<           zero(m) = t(l)+h1*y(i) >*/
	    zero[*m] = t[l] + h1 * y[i__ - 1];
/*<  150    continue >*/
L150:
	    ;
	}
/*<  200    a0 = b0 >*/
L200:
	a0 = b0;
/*<         ah = bh >*/
	ah = bh;
/*<         z1 = z3 >*/
	z1 = z3;
/*<         nz1 = nz3 >*/
	nz1 = nz3;
/*<  300  continue >*/
/* L300: */
    }
/*  the zeros of s(x) are arranged in increasing order. */
/*<       if(m.lt.2) go to 800 >*/
    if (*m < 2) {
	goto L800;
    }
/*<       do 400 i=2,m >*/
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         j = i >*/
	j = i__;
/*<  350    j1 = j-1 >*/
L350:
	j1 = j - 1;
/*<         if(j1.eq.0) go to 400 >*/
	if (j1 == 0) {
	    goto L400;
	}
/*<         if(zero(j).ge.zero(j1)) go to 400 >*/
	if (zero[j] >= zero[j1]) {
	    goto L400;
	}
/*<         zz = zero(j) >*/
	zz = zero[j];
/*<         zero(j) = zero(j1) >*/
	zero[j] = zero[j1];
/*<         zero(j1) = zz >*/
	zero[j1] = zz;
/*<         j = j1 >*/
	j = j1;
/*<         go to 350 >*/
	goto L350;
/*<  400  continue >*/
L400:
	;
    }
/*<       j = m >*/
    j = *m;
/*<       m = 1 >*/
    *m = 1;
/*<       do 500 i=2,j >*/
    i__1 = j;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         if(zero(i).eq.zero(m)) go to 500 >*/
	if (zero[i__] == zero[*m]) {
	    goto L500;
	}
/*<         m = m+1 >*/
	++(*m);
/*<         zero(m) = zero(i) >*/
	zero[*m] = zero[i__];
/*<  500  continue >*/
L500:
	;
    }
/*<       go to 800 >*/
    goto L800;
/*<  700  ier = 1 >*/
L700:
    *ier = 1;
/*<  800  return >*/
L800:
    return 0;
/*<       end >*/
} /* sproot_ */

