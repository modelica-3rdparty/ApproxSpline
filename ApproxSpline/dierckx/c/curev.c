/* curev.f -- translated by f2c (version 20061008).
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

/*<       subroutine curev(idim,t,n,c,nc,k,u,m,x,mx,ier) >*/
/* Subroutine */ int curev_(integer *idim, doublereal *t, integer *n, 
	doublereal *c__, integer *nc, integer *k, doublereal *u, integer *m, 
	doublereal *x, integer *mx, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, j1, k1, l1, jj;
    static doublereal tb;
    static integer ll;
    static doublereal te;
    static integer mm;
    static doublereal sp;
    static integer nk1;
    static doublereal arg;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  subroutine curev evaluates in a number of points u(i),i=1,2,...,m */
/*  a spline curve s(u) of degree k and dimension idim, given in its */
/*  b-spline representation. */

/*  calling sequence: */
/*     call curev(idim,t,n,c,nc,k,u,m,x,mx,ier) */

/*  input parameters: */
/*    idim : integer, giving the dimension of the spline curve. */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(u). */
/*    c    : array,length nc, which contains the b-spline coefficients. */
/*    nc   : integer, giving the total number of coefficients of s(u). */
/*    k    : integer, giving the degree of s(u). */
/*    u    : array,length m, which contains the points where s(u) must */
/*           be evaluated. */
/*    m    : integer, giving the number of points where s(u) must be */
/*           evaluated. */
/*    mx   : integer, giving the dimension of the array x. mx >= m*idim */

/*  output parameters: */
/*    x    : array,length mx,giving the value of s(u) at the different */
/*           points. x(idim*(i-1)+j) will contain the j-th coordinate */
/*           of the i-th point on the curve. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    m >= 1 */
/*    mx >= m*idim */
/*    t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1. */

/*  other subroutines required: fpbspl. */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths */
/*                applics 10 (1972) 134-149. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/*  ..scalar arguments.. */
/*<       integer idim,n,nc,k,m,mx,ier >*/
/*  ..array arguments.. */
/*<       real t(n),c(nc),u(m),x(mx) >*/
/*  ..local scalars.. */
/*<       integer i,j,jj,j1,k1,l,ll,l1,mm,nk1 >*/
/*<       real arg,sp,tb,te >*/
/*  ..local array.. */
/*<       real h(6) >*/
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
/*<       ier = 10 >*/
    /* Parameter adjustments */
    --t;
    --c__;
    --u;
    --x;

    /* Function Body */
    *ier = 10;
/*<       if(m-1) 100,30,10 >*/
    if ((i__1 = *m - 1) < 0) {
	goto L100;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L10;
    }
/*<   10  do 20 i=2,m >*/
L10:
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         if(u(i).lt.u(i-1)) go to 100 >*/
	if (u[i__] < u[i__ - 1]) {
	    goto L100;
	}
/*<   20  continue >*/
/* L20: */
    }
/*<   30  if(mx.lt.(m*idim)) go to 100 >*/
L30:
    if (*mx < *m * *idim) {
	goto L100;
    }
/*<       ier = 0 >*/
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
/*<       l = k1 >*/
    l = k1;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*  main loop for the different points. */
/*<       mm = 0 >*/
    mm = 0;
/*<       do 80 i=1,m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch a new u-value arg. */
/*<         arg = u(i) >*/
	arg = u[i__];
/*<         if(arg.lt.tb) arg = tb >*/
	if (arg < tb) {
	    arg = tb;
	}
/*<         if(arg.gt.te) arg = te >*/
	if (arg > te) {
	    arg = te;
	}
/*  search for knot interval t(l) <= arg < t(l+1) */
/*<   40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50 >*/
L40:
	if (arg < t[l1] || l == nk1) {
	    goto L50;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 40 >*/
	goto L40;
/*  evaluate the non-zero b-splines at arg. */
/*<   50    call fpbspl(t,n,k,arg,l,h) >*/
L50:
	fpbspl_(&t[1], n, k, &arg, &l, h__);
/*  find the value of s(u) at u=arg. */
/*<         ll = l-k1 >*/
	ll = l - k1;
/*<         do 70 j1=1,idim >*/
	i__2 = *idim;
	for (j1 = 1; j1 <= i__2; ++j1) {
/*<           jj = ll >*/
	    jj = ll;
/*<           sp = 0. >*/
	    sp = 0.;
/*<           do 60 j=1,k1 >*/
	    i__3 = k1;
	    for (j = 1; j <= i__3; ++j) {
/*<             jj = jj+1 >*/
		++jj;
/*<             sp = sp+c(jj)*h(j) >*/
		sp += c__[jj] * h__[j - 1];
/*<   60      continue >*/
/* L60: */
	    }
/*<           mm = mm+1 >*/
	    ++mm;
/*<           x(mm) = sp >*/
	    x[mm] = sp;
/*<           ll = ll+n >*/
	    ll += *n;
/*<   70    continue >*/
/* L70: */
	}
/*<   80  continue >*/
/* L80: */
    }
/*<  100  return >*/
L100:
    return 0;
/*<       end >*/
} /* curev_ */

