/* evapol.f -- translated by f2c (version 20061008).
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

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__8 = 8;
static integer c__2 = 2;

/*<       real function evapol(tu,nu,tv,nv,c,rad,x,y) >*/
doublereal evapol_(doublereal *tu, integer *nu, doublereal *tv, integer *nv, 
	doublereal *c__, D_fp rad, doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double atan2(doublereal, doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal f, r__, u, v;
    static integer ier;
    static doublereal one, wrk[8], dist;
    static integer iwrk[2];
    extern /* Subroutine */ int bispev_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *);

/*  function program evacir evaluates the function f(x,y) = s(u,v), */
/*  defined through the transformation */
/*      x = u*rad(v)*cos(v)    y = u*rad(v)*sin(v) */
/*  and where s(u,v) is a bicubic spline ( 0<=u<=1 , -pi<=v<=pi ), given */
/*  in its standard b-spline representation. */

/*  calling sequence: */
/*     f = evapol(tu,nu,tv,nv,c,rad,x,y) */

/*  input parameters: */
/*   tu    : real array, length nu, which contains the position of the */
/*           knots in the u-direction. */
/*   nu    : integer, giving the total number of knots in the u-direction */
/*   tv    : real array, length nv, which contains the position of the */
/*           knots in the v-direction. */
/*   nv    : integer, giving the total number of knots in the v-direction */
/*   c     : real array, length (nu-4)*(nv-4), which contains the */
/*           b-spline coefficients. */
/*   rad   : real function subprogram, defining the boundary of the */
/*           approximation domain. must be declared external in the */
/*           calling (sub)-program */
/*   x,y   : real values. */
/*           before entry x and y must be set to the co-ordinates of */
/*           the point where f(x,y) must be evaluated. */

/*  output parameter: */
/*   f     : real */
/*           on exit f contains the value of f(x,y) */

/*  other subroutines required: */
/*    bispev,fpbisp,fpbspl */

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

/*  latest update : march 1989 */

/*  ..scalar arguments.. */
/*<       integer nu,nv >*/
/*<       real x,y >*/
/*  ..array arguments.. */
/*<       real tu(nu),tv(nv),c((nu-4)*(nv-4)) >*/
/*  ..user specified function */
/*<       real rad >*/
/*  ..local scalars.. */
/*<       integer ier >*/
/*<       real u,v,r,f,one,dist >*/
/*  ..local arrays */
/*<       real wrk(8) >*/
/*<       integer iwrk(2) >*/
/*  ..function references */
/*<       real atan2,sqrt >*/
/*  .. */
/*  calculate the (u,v)-coordinates of the given point. */
/*<       one = 1 >*/
    /* Parameter adjustments */
    --tu;
    --c__;
    --tv;

    /* Function Body */
    one = 1.;
/*<       u = 0. >*/
    u = 0.;
/*<       v = 0. >*/
    v = 0.;
/*<       dist = x**2+y**2 >*/
/* Computing 2nd power */
    d__1 = *x;
/* Computing 2nd power */
    d__2 = *y;
    dist = d__1 * d__1 + d__2 * d__2;
/*<       if(dist.le.0.) go to 10 >*/
    if (dist <= 0.) {
	goto L10;
    }
/*<       v = atan2(y,x) >*/
    v = atan2(*y, *x);
/*<       r = rad(v) >*/
    r__ = (*rad)(&v);
/*<       if(r.le.0.) go to 10 >*/
    if (r__ <= 0.) {
	goto L10;
    }
/*<       u = sqrt(dist)/r >*/
    u = sqrt(dist) / r__;
/*<       if(u.gt.one) u = one >*/
    if (u > one) {
	u = one;
    }
/*  evaluate s(u,v) */
/*<   10  call bispev(tu,nu,tv,nv,c,3,3,u,1,v,1,f,wrk,8,iwrk,2,ier) >*/
L10:
    bispev_(&tu[1], nu, &tv[1], nv, &c__[1], &c__3, &c__3, &u, &c__1, &v, &
	    c__1, &f, wrk, &c__8, iwrk, &c__2, &ier);
/*<       evapol = f >*/
    ret_val = f;
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* evapol_ */

