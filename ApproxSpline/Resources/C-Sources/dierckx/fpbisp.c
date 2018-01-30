/* fpbisp.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly) >*/
/* Subroutine */ int fpbisp_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, doublereal *x,
	 integer *mx, doublereal *y, integer *my, doublereal *z__, doublereal 
	*wx, doublereal *wy, integer *lx, integer *ly)
{
    /* System generated locals */
    integer wx_dim1, wx_offset, wy_dim1, wy_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, m, i1, j1, l1, l2;
    static doublereal tb, te, sp;
    static integer kx1, ky1;
    static doublereal arg;
    static integer nkx1, nky1;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  ..scalar arguments.. */
/*<       integer nx,ny,kx,ky,mx,my >*/
/*  ..array arguments.. */
/*<       integer lx(mx),ly(my) >*/
/*<    >*/
/*  ..local scalars.. */
/*<       integer kx1,ky1,l,l1,l2,m,nkx1,nky1 >*/
/*<       real arg,sp,tb,te >*/
/*  ..local arrays.. */
/*<       real h(6) >*/
/*  ..subroutine references.. */
/*    fpbspl */
/*  .. */
/*<       kx1 = kx+1 >*/
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --lx;
    wx_dim1 = *mx;
    wx_offset = 1 + wx_dim1;
    wx -= wx_offset;
    --x;
    --ly;
    wy_dim1 = *my;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    --z__;
    --y;

    /* Function Body */
    kx1 = *kx + 1;
/*<       nkx1 = nx-kx1 >*/
    nkx1 = *nx - kx1;
/*<       tb = tx(kx1) >*/
    tb = tx[kx1];
/*<       te = tx(nkx1+1) >*/
    te = tx[nkx1 + 1];
/*<       l = kx1 >*/
    l = kx1;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*<       do 40 i=1,mx >*/
    i__1 = *mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         arg = x(i) >*/
	arg = x[i__];
/*<         if(arg.lt.tb) arg = tb >*/
	if (arg < tb) {
	    arg = tb;
	}
/*<         if(arg.gt.te) arg = te >*/
	if (arg > te) {
	    arg = te;
	}
/*<   10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20 >*/
L10:
	if (arg < tx[l1] || l == nkx1) {
	    goto L20;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 10 >*/
	goto L10;
/*<   20    call fpbspl(tx,nx,kx,arg,l,h) >*/
L20:
	fpbspl_(&tx[1], nx, kx, &arg, &l, h__);
/*<         lx(i) = l-kx1 >*/
	lx[i__] = l - kx1;
/*<         do 30 j=1,kx1 >*/
	i__2 = kx1;
	for (j = 1; j <= i__2; ++j) {
/*<           wx(i,j) = h(j) >*/
	    wx[i__ + j * wx_dim1] = h__[j - 1];
/*<   30    continue >*/
/* L30: */
	}
/*<   40  continue >*/
/* L40: */
    }
/*<       ky1 = ky+1 >*/
    ky1 = *ky + 1;
/*<       nky1 = ny-ky1 >*/
    nky1 = *ny - ky1;
/*<       tb = ty(ky1) >*/
    tb = ty[ky1];
/*<       te = ty(nky1+1) >*/
    te = ty[nky1 + 1];
/*<       l = ky1 >*/
    l = ky1;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*<       do 80 i=1,my >*/
    i__1 = *my;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         arg = y(i) >*/
	arg = y[i__];
/*<         if(arg.lt.tb) arg = tb >*/
	if (arg < tb) {
	    arg = tb;
	}
/*<         if(arg.gt.te) arg = te >*/
	if (arg > te) {
	    arg = te;
	}
/*<   50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60 >*/
L50:
	if (arg < ty[l1] || l == nky1) {
	    goto L60;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 50 >*/
	goto L50;
/*<   60    call fpbspl(ty,ny,ky,arg,l,h) >*/
L60:
	fpbspl_(&ty[1], ny, ky, &arg, &l, h__);
/*<         ly(i) = l-ky1 >*/
	ly[i__] = l - ky1;
/*<         do 70 j=1,ky1 >*/
	i__2 = ky1;
	for (j = 1; j <= i__2; ++j) {
/*<           wy(i,j) = h(j) >*/
	    wy[i__ + j * wy_dim1] = h__[j - 1];
/*<   70    continue >*/
/* L70: */
	}
/*<   80  continue >*/
/* L80: */
    }
/*<       m = 0 >*/
    m = 0;
/*<       do 130 i=1,mx >*/
    i__1 = *mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         l = lx(i)*nky1 >*/
	l = lx[i__] * nky1;
/*<         do 90 i1=1,kx1 >*/
	i__2 = kx1;
	for (i1 = 1; i1 <= i__2; ++i1) {
/*<           h(i1) = wx(i,i1) >*/
	    h__[i1 - 1] = wx[i__ + i1 * wx_dim1];
/*<   90    continue >*/
/* L90: */
	}
/*<         do 120 j=1,my >*/
	i__2 = *my;
	for (j = 1; j <= i__2; ++j) {
/*<           l1 = l+ly(j) >*/
	    l1 = l + ly[j];
/*<           sp = 0. >*/
	    sp = 0.;
/*<           do 110 i1=1,kx1 >*/
	    i__3 = kx1;
	    for (i1 = 1; i1 <= i__3; ++i1) {
/*<             l2 = l1 >*/
		l2 = l1;
/*<             do 100 j1=1,ky1 >*/
		i__4 = ky1;
		for (j1 = 1; j1 <= i__4; ++j1) {
/*<               l2 = l2+1 >*/
		    ++l2;
/*<               sp = sp+c(l2)*h(i1)*wy(j,j1) >*/
		    sp += c__[l2] * h__[i1 - 1] * wy[j + j1 * wy_dim1];
/*<  100        continue >*/
/* L100: */
		}
/*<             l1 = l1+nky1 >*/
		l1 += nky1;
/*<  110      continue >*/
/* L110: */
	    }
/*<           m = m+1 >*/
	    ++m;
/*<           z(m) = sp >*/
	    z__[m] = sp;
/*<  120    continue >*/
/* L120: */
	}
/*<  130  continue >*/
/* L130: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpbisp_ */

