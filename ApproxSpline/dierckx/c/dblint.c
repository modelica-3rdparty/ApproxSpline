/* dblint.f -- translated by f2c (version 20061008).
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

/*<       real function dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk) >*/
doublereal dblint_(doublereal *tx, integer *nx, doublereal *ty, integer *ny, 
	doublereal *c__, integer *kx, integer *ky, doublereal *xb, doublereal 
	*xe, doublereal *yb, doublereal *ye, doublereal *wrk)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, l, m;
    static doublereal res;
    static integer nkx1, nky1;
    extern /* Subroutine */ int fpintb_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *);

/*  function dblint calculates the double integral */
/*         / xe  / ye */
/*        |     |      s(x,y) dx dy */
/*    xb /  yb / */
/*  with s(x,y) a bivariate spline of degrees kx and ky, given in the */
/*  b-spline representation. */

/*  calling sequence: */
/*     aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk) */

/*  input parameters: */
/*   tx    : real array, length nx, which contains the position of the */
/*           knots in the x-direction. */
/*   nx    : integer, giving the total number of knots in the x-direction */
/*   ty    : real array, length ny, which contains the position of the */
/*           knots in the y-direction. */
/*   ny    : integer, giving the total number of knots in the y-direction */
/*   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the */
/*           b-spline coefficients. */
/*   kx,ky : integer values, giving the degrees of the spline. */
/*   xb,xe : real values, containing the boundaries of the integration */
/*   yb,ye   domain. s(x,y) is considered to be identically zero out- */
/*           side the rectangle (tx(kx+1),tx(nx-kx))*(ty(ky+1),ty(ny-ky)) */

/*  output parameters: */
/*   aint  : real , containing the double integral of s(x,y). */
/*   wrk   : real array of dimension at least (nx+ny-kx-ky-2). */
/*           used as working space. */
/*           on exit, wrk(i) will contain the integral */
/*                / xe */
/*               | ni,kx+1(x) dx , i=1,2,...,nx-kx-1 */
/*           xb / */
/*           with ni,kx+1(x) the normalized b-spline defined on */
/*           the knots tx(i),...,tx(i+kx+1) */
/*           wrk(j+nx-kx-1) will contain the integral */
/*                / ye */
/*               | nj,ky+1(y) dy , j=1,2,...,ny-ky-1 */
/*           yb / */
/*           with nj,ky+1(y) the normalized b-spline defined on */
/*           the knots ty(j),...,ty(j+ky+1) */

/*  other subroutines required: fpintb */

/*  references : */
/*    gaffney p.w. : the calculation of indefinite integrals of b-splines */
/*                   j. inst. maths applics 17 (1976) 37-41. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1989 */

/*  ..scalar arguments.. */
/*<       integer nx,ny,kx,ky >*/
/*<       real xb,xe,yb,ye >*/
/*  ..array arguments.. */
/*<       real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),wrk(nx+ny-kx-ky-2) >*/
/*  ..local scalars.. */
/*<       integer i,j,l,m,nkx1,nky1 >*/
/*<       real res >*/
/*  .. */
/*<       nkx1 = nx-kx-1 >*/
    /* Parameter adjustments */
    --tx;
    --ty;
    --wrk;
    --c__;

    /* Function Body */
    nkx1 = *nx - *kx - 1;
/*<       nky1 = ny-ky-1 >*/
    nky1 = *ny - *ky - 1;
/*  we calculate the integrals of the normalized b-splines ni,kx+1(x) */
/*<       call fpintb(tx,nx,wrk,nkx1,xb,xe) >*/
    fpintb_(&tx[1], nx, &wrk[1], &nkx1, xb, xe);
/*  we calculate the integrals of the normalized b-splines nj,ky+1(y) */
/*<       call fpintb(ty,ny,wrk(nkx1+1),nky1,yb,ye) >*/
    fpintb_(&ty[1], ny, &wrk[nkx1 + 1], &nky1, yb, ye);
/*  calculate the integral of s(x,y) */
/*<       dblint = 0. >*/
    ret_val = 0.;
/*<       do 200 i=1,nkx1 >*/
    i__1 = nkx1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         res = wrk(i) >*/
	res = wrk[i__];
/*<         if(res.eq.0.) go to 200 >*/
	if (res == 0.) {
	    goto L200;
	}
/*<         m = (i-1)*nky1 >*/
	m = (i__ - 1) * nky1;
/*<         l = nkx1 >*/
	l = nkx1;
/*<         do 100 j=1,nky1 >*/
	i__2 = nky1;
	for (j = 1; j <= i__2; ++j) {
/*<           m = m+1 >*/
	    ++m;
/*<           l = l+1 >*/
	    ++l;
/*<           dblint = dblint+res*wrk(l)*c(m) >*/
	    ret_val += res * wrk[l] * c__[m];
/*<  100    continue >*/
/* L100: */
	}
/*<  200  continue >*/
L200:
	;
    }
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* dblint_ */

