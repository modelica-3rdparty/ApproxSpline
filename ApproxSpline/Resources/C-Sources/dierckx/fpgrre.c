/* fpgrre.f -- translated by f2c (version 20061008).
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

/*<    >*/
/* Subroutine */ int fpgrre_(integer *ifsx, integer *ifsy, integer *ifbx, 
	integer *ifby, doublereal *x, integer *mx, doublereal *y, integer *my,
	 doublereal *z__, integer *mz, integer *kx, integer *ky, doublereal *
	tx, integer *nx, doublereal *ty, integer *ny, doublereal *p, 
	doublereal *c__, integer *nc, doublereal *fp, doublereal *fpx, 
	doublereal *fpy, integer *mm, integer *mynx, integer *kx1, integer *
	kx2, integer *ky1, integer *ky2, doublereal *spx, doublereal *spy, 
	doublereal *right, doublereal *q, doublereal *ax, doublereal *ay, 
	doublereal *bx, doublereal *by, integer *nrx, integer *nry)
{
    /* System generated locals */
    integer spx_dim1, spx_offset, spy_dim1, spy_offset, ax_dim1, ax_offset, 
	    bx_dim1, bx_offset, ay_dim1, ay_offset, by_dim1, by_offset, i__1, 
	    i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static doublereal h__[7];
    static integer i__, j, k, l, i1, i2, i3, k1, k2, l1, l2, n1, ic, iq, it, 
	    iz;
    static doublereal fac, arg, one, cos__, sin__, piv;
    static integer nk1x, nk1y;
    static doublereal half;
    static integer ncof;
    static doublereal term, pinv;
    static integer irot, numx, numy, numx1, numy1, nrold;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *);
    static integer ibandx, ibandy;
    extern /* Subroutine */ int fpdisc_(doublereal *, integer *, integer *, 
	    doublereal *, integer *), fpbspl_(doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer number;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer nroldx, nroldy;

/*  .. */
/*  ..scalar arguments.. */
/*<       real p,fp >*/
/*<    >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer nrx(mx),nry(my) >*/
/*  ..local scalars.. */
/*<       real arg,cos,fac,pinv,piv,sin,term,one,half >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real h(7) >*/
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpgivs,fpdisc,fprota */
/*  .. */
/*  the b-spline coefficients of the smoothing spline are calculated as */
/*  the least-squares solution of the over-determined linear system of */
/*  equations  (ay) c (ax)' = q       where */

/*               |   (spx)    |            |   (spy)    | */
/*        (ax) = | ---------- |     (ay) = | ---------- | */
/*               | (1/p) (bx) |            | (1/p) (by) | */

/*                                | z  ' 0 | */
/*                            q = | ------ | */
/*                                | 0  ' 0 | */

/*  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the */
/*                b-spline coefficients. */
/*       z      : the my x mx matrix which contains the function values. */
/*       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation */
/*                matrices according to the least-squares problems in */
/*                the x- and y-direction. */
/*       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1) */
/*                matrices which contain the discontinuity jumps of the */
/*                derivatives of the b-splines in the x- and y-direction. */
/*<       one = 1 >*/
    /* Parameter adjustments */
    --nrx;
    --x;
    --nry;
    --y;
    --z__;
    --fpx;
    --tx;
    --fpy;
    --ty;
    --c__;
    --right;
    --q;
    spx_dim1 = *mx;
    spx_offset = 1 + spx_dim1;
    spx -= spx_offset;
    bx_dim1 = *nx;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
    ax_dim1 = *nx;
    ax_offset = 1 + ax_dim1;
    ax -= ax_offset;
    spy_dim1 = *my;
    spy_offset = 1 + spy_dim1;
    spy -= spy_offset;
    by_dim1 = *ny;
    by_offset = 1 + by_dim1;
    by -= by_offset;
    ay_dim1 = *ny;
    ay_offset = 1 + ay_dim1;
    ay -= ay_offset;

    /* Function Body */
    one = 1.;
/*<       half = 0.5 >*/
    half = .5;
/*<       nk1x = nx-kx1 >*/
    nk1x = *nx - *kx1;
/*<       nk1y = ny-ky1 >*/
    nk1y = *ny - *ky1;
/*<       if(p.gt.0.) pinv = one/p >*/
    if (*p > 0.) {
	pinv = one / *p;
    }
/*  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on */
/*  the value of p whether the matrices (spx),(spy),(bx) and (by) still */
/*  must be determined. */
/*<       if(ifsx.ne.0) go to 50 >*/
    if (*ifsx != 0) {
	goto L50;
    }
/*  calculate the non-zero elements of the matrix (spx) which is the */
/*  observation matrix according to the least-squares spline approximat- */
/*  ion problem in the x-direction. */
/*<       l = kx1 >*/
    l = *kx1;
/*<       l1 = kx2 >*/
    l1 = *kx2;
/*<       number = 0 >*/
    number = 0;
/*<       do 40 it=1,mx >*/
    i__1 = *mx;
    for (it = 1; it <= i__1; ++it) {
/*<         arg = x(it) >*/
	arg = x[it];
/*<   10    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20 >*/
L10:
	if (arg < tx[l1] || l == nk1x) {
	    goto L20;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         number = number+1 >*/
	++number;
/*<         go to 10 >*/
	goto L10;
/*<   20    call fpbspl(tx,nx,kx,arg,l,h) >*/
L20:
	fpbspl_(&tx[1], nx, kx, &arg, &l, h__);
/*<         do 30 i=1,kx1 >*/
	i__2 = *kx1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           spx(it,i) = h(i) >*/
	    spx[it + i__ * spx_dim1] = h__[i__ - 1];
/*<   30    continue >*/
/* L30: */
	}
/*<         nrx(it) = number >*/
	nrx[it] = number;
/*<   40  continue >*/
/* L40: */
    }
/*<       ifsx = 1 >*/
    *ifsx = 1;
/*<   50  if(ifsy.ne.0) go to 100 >*/
L50:
    if (*ifsy != 0) {
	goto L100;
    }
/*  calculate the non-zero elements of the matrix (spy) which is the */
/*  observation matrix according to the least-squares spline approximat- */
/*  ion problem in the y-direction. */
/*<       l = ky1 >*/
    l = *ky1;
/*<       l1 = ky2 >*/
    l1 = *ky2;
/*<       number = 0 >*/
    number = 0;
/*<       do 90 it=1,my >*/
    i__1 = *my;
    for (it = 1; it <= i__1; ++it) {
/*<         arg = y(it) >*/
	arg = y[it];
/*<   60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70 >*/
L60:
	if (arg < ty[l1] || l == nk1y) {
	    goto L70;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         number = number+1 >*/
	++number;
/*<         go to 60 >*/
	goto L60;
/*<   70    call fpbspl(ty,ny,ky,arg,l,h) >*/
L70:
	fpbspl_(&ty[1], ny, ky, &arg, &l, h__);
/*<         do 80 i=1,ky1 >*/
	i__2 = *ky1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           spy(it,i) = h(i) >*/
	    spy[it + i__ * spy_dim1] = h__[i__ - 1];
/*<   80    continue >*/
/* L80: */
	}
/*<         nry(it) = number >*/
	nry[it] = number;
/*<   90  continue >*/
/* L90: */
    }
/*<       ifsy = 1 >*/
    *ifsy = 1;
/*<  100  if(p.le.0.) go to 120 >*/
L100:
    if (*p <= 0.) {
	goto L120;
    }
/*  calculate the non-zero elements of the matrix (bx). */
/*<       if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110 >*/
    if (*ifbx != 0 || *nx == *kx1 << 1) {
	goto L110;
    }
/*<       call fpdisc(tx,nx,kx2,bx,nx) >*/
    fpdisc_(&tx[1], nx, kx2, &bx[bx_offset], nx);
/*<       ifbx = 1 >*/
    *ifbx = 1;
/*  calculate the non-zero elements of the matrix (by). */
/*<  110  if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120 >*/
L110:
    if (*ifby != 0 || *ny == *ky1 << 1) {
	goto L120;
    }
/*<       call fpdisc(ty,ny,ky2,by,ny) >*/
    fpdisc_(&ty[1], ny, ky2, &by[by_offset], ny);
/*<       ifby = 1 >*/
    *ifby = 1;
/*  reduce the matrix (ax) to upper triangular form (rx) using givens */
/*  rotations. apply the same transformations to the rows of matrix q */
/*  to obtain the my x (nx-kx-1) matrix g. */
/*  store matrix (rx) into (ax) and g into q. */
/*<  120  l = my*nk1x >*/
L120:
    l = *my * nk1x;
/*  initialization. */
/*<       do 130 i=1,l >*/
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         q(i) = 0. >*/
	q[i__] = 0.;
/*<  130  continue >*/
/* L130: */
    }
/*<       do 140 i=1,nk1x >*/
    i__1 = nk1x;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         do 140 j=1,kx2 >*/
	i__2 = *kx2;
	for (j = 1; j <= i__2; ++j) {
/*<           ax(i,j) = 0. >*/
	    ax[i__ + j * ax_dim1] = 0.;
/*<  140  continue >*/
/* L140: */
	}
    }
/*<       l = 0 >*/
    l = 0;
/*<       nrold = 0 >*/
    nrold = 0;
/*  ibandx denotes the bandwidth of the matrices (ax) and (rx). */
/*<       ibandx = kx1 >*/
    ibandx = *kx1;
/*<       do 270 it=1,mx >*/
    i__2 = *mx;
    for (it = 1; it <= i__2; ++it) {
/*<         number = nrx(it) >*/
	number = nrx[it];
/*<  150    if(nrold.eq.number) go to 180 >*/
L150:
	if (nrold == number) {
	    goto L180;
	}
/*<         if(p.le.0.) go to 260 >*/
	if (*p <= 0.) {
	    goto L260;
	}
/*<         ibandx = kx2 >*/
	ibandx = *kx2;
/*  fetch a new row of matrix (bx). */
/*<         n1 = nrold+1 >*/
	n1 = nrold + 1;
/*<         do 160 j=1,kx2 >*/
	i__1 = *kx2;
	for (j = 1; j <= i__1; ++j) {
/*<           h(j) = bx(n1,j)*pinv >*/
	    h__[j - 1] = bx[n1 + j * bx_dim1] * pinv;
/*<  160    continue >*/
/* L160: */
	}
/*  find the appropriate column of q. */
/*<         do 170 j=1,my >*/
	i__1 = *my;
	for (j = 1; j <= i__1; ++j) {
/*<           right(j) = 0. >*/
	    right[j] = 0.;
/*<  170    continue >*/
/* L170: */
	}
/*<         irot = nrold >*/
	irot = nrold;
/*<         go to 210 >*/
	goto L210;
/*  fetch a new row of matrix (spx). */
/*<  180    h(ibandx) = 0. >*/
L180:
	h__[ibandx - 1] = 0.;
/*<         do 190 j=1,kx1 >*/
	i__1 = *kx1;
	for (j = 1; j <= i__1; ++j) {
/*<           h(j) = spx(it,j) >*/
	    h__[j - 1] = spx[it + j * spx_dim1];
/*<  190    continue >*/
/* L190: */
	}
/*  find the appropriate column of q. */
/*<         do 200 j=1,my >*/
	i__1 = *my;
	for (j = 1; j <= i__1; ++j) {
/*<           l = l+1 >*/
	    ++l;
/*<           right(j) = z(l) >*/
	    right[j] = z__[l];
/*<  200    continue >*/
/* L200: */
	}
/*<         irot = number >*/
	irot = number;
/*  rotate the new row of matrix (ax) into triangle. */
/*<  210    do 240 i=1,ibandx >*/
L210:
	i__1 = ibandx;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<           irot = irot+1 >*/
	    ++irot;
/*<           piv = h(i) >*/
	    piv = h__[i__ - 1];
/*<           if(piv.eq.0.) go to 240 >*/
	    if (piv == 0.) {
		goto L240;
	    }
/*  calculate the parameters of the givens transformation. */
/*<           call fpgivs(piv,ax(irot,1),cos,sin) >*/
	    fpgivs_(&piv, &ax[irot + ax_dim1], &cos__, &sin__);
/*  apply that transformation to the rows of matrix q. */
/*<           iq = (irot-1)*my >*/
	    iq = (irot - 1) * *my;
/*<           do 220 j=1,my >*/
	    i__3 = *my;
	    for (j = 1; j <= i__3; ++j) {
/*<             iq = iq+1 >*/
		++iq;
/*<             call fprota(cos,sin,right(j),q(iq)) >*/
		fprota_(&cos__, &sin__, &right[j], &q[iq]);
/*<  220      continue >*/
/* L220: */
	    }
/*  apply that transformation to the columns of (ax). */
/*<           if(i.eq.ibandx) go to 250 >*/
	    if (i__ == ibandx) {
		goto L250;
	    }
/*<           i2 = 1 >*/
	    i2 = 1;
/*<           i3 = i+1 >*/
	    i3 = i__ + 1;
/*<           do 230 j=i3,ibandx >*/
	    i__3 = ibandx;
	    for (j = i3; j <= i__3; ++j) {
/*<             i2 = i2+1 >*/
		++i2;
/*<             call fprota(cos,sin,h(j),ax(irot,i2)) >*/
		fprota_(&cos__, &sin__, &h__[j - 1], &ax[irot + i2 * ax_dim1])
			;
/*<  230      continue >*/
/* L230: */
	    }
/*<  240    continue >*/
L240:
	    ;
	}
/*<  250    if(nrold.eq.number) go to 270 >*/
L250:
	if (nrold == number) {
	    goto L270;
	}
/*<  260    nrold = nrold+1 >*/
L260:
	++nrold;
/*<         go to 150 >*/
	goto L150;
/*<  270  continue >*/
L270:
	;
    }
/*  reduce the matrix (ay) to upper triangular form (ry) using givens */
/*  rotations. apply the same transformations to the columns of matrix g */
/*  to obtain the (ny-ky-1) x (nx-kx-1) matrix h. */
/*  store matrix (ry) into (ay) and h into c. */
/*<       ncof = nk1x*nk1y >*/
    ncof = nk1x * nk1y;
/*  initialization. */
/*<       do 280 i=1,ncof >*/
    i__2 = ncof;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         c(i) = 0. >*/
	c__[i__] = 0.;
/*<  280  continue >*/
/* L280: */
    }
/*<       do 290 i=1,nk1y >*/
    i__2 = nk1y;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         do 290 j=1,ky2 >*/
	i__1 = *ky2;
	for (j = 1; j <= i__1; ++j) {
/*<           ay(i,j) = 0. >*/
	    ay[i__ + j * ay_dim1] = 0.;
/*<  290  continue >*/
/* L290: */
	}
    }
/*<       nrold = 0 >*/
    nrold = 0;
/*  ibandy denotes the bandwidth of the matrices (ay) and (ry). */
/*<       ibandy = ky1 >*/
    ibandy = *ky1;
/*<       do 420 it=1,my >*/
    i__1 = *my;
    for (it = 1; it <= i__1; ++it) {
/*<         number = nry(it) >*/
	number = nry[it];
/*<  300    if(nrold.eq.number) go to 330 >*/
L300:
	if (nrold == number) {
	    goto L330;
	}
/*<         if(p.le.0.) go to 410 >*/
	if (*p <= 0.) {
	    goto L410;
	}
/*<         ibandy = ky2 >*/
	ibandy = *ky2;
/*  fetch a new row of matrix (by). */
/*<         n1 = nrold+1 >*/
	n1 = nrold + 1;
/*<         do 310 j=1,ky2 >*/
	i__2 = *ky2;
	for (j = 1; j <= i__2; ++j) {
/*<           h(j) = by(n1,j)*pinv >*/
	    h__[j - 1] = by[n1 + j * by_dim1] * pinv;
/*<  310    continue >*/
/* L310: */
	}
/*  find the appropiate row of g. */
/*<         do 320 j=1,nk1x >*/
	i__2 = nk1x;
	for (j = 1; j <= i__2; ++j) {
/*<           right(j) = 0. >*/
	    right[j] = 0.;
/*<  320    continue >*/
/* L320: */
	}
/*<         irot = nrold >*/
	irot = nrold;
/*<         go to 360 >*/
	goto L360;
/*  fetch a new row of matrix (spy) */
/*<  330    h(ibandy) = 0. >*/
L330:
	h__[ibandy - 1] = 0.;
/*<         do 340 j=1,ky1 >*/
	i__2 = *ky1;
	for (j = 1; j <= i__2; ++j) {
/*<           h(j) = spy(it,j) >*/
	    h__[j - 1] = spy[it + j * spy_dim1];
/*<  340    continue >*/
/* L340: */
	}
/*  find the appropiate row of g. */
/*<         l = it >*/
	l = it;
/*<         do 350 j=1,nk1x >*/
	i__2 = nk1x;
	for (j = 1; j <= i__2; ++j) {
/*<           right(j) = q(l) >*/
	    right[j] = q[l];
/*<           l = l+my >*/
	    l += *my;
/*<  350    continue >*/
/* L350: */
	}
/*<         irot = number >*/
	irot = number;
/*  rotate the new row of matrix (ay) into triangle. */
/*<  360    do 390 i=1,ibandy >*/
L360:
	i__2 = ibandy;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           irot = irot+1 >*/
	    ++irot;
/*<           piv = h(i) >*/
	    piv = h__[i__ - 1];
/*<           if(piv.eq.0.) go to 390 >*/
	    if (piv == 0.) {
		goto L390;
	    }
/*  calculate the parameters of the givens transformation. */
/*<           call fpgivs(piv,ay(irot,1),cos,sin) >*/
	    fpgivs_(&piv, &ay[irot + ay_dim1], &cos__, &sin__);
/*  apply that transformation to the colums of matrix g. */
/*<           ic = irot >*/
	    ic = irot;
/*<           do 370 j=1,nk1x >*/
	    i__3 = nk1x;
	    for (j = 1; j <= i__3; ++j) {
/*<             call fprota(cos,sin,right(j),c(ic)) >*/
		fprota_(&cos__, &sin__, &right[j], &c__[ic]);
/*<             ic = ic+nk1y >*/
		ic += nk1y;
/*<  370      continue >*/
/* L370: */
	    }
/*  apply that transformation to the columns of matrix (ay). */
/*<           if(i.eq.ibandy) go to 400 >*/
	    if (i__ == ibandy) {
		goto L400;
	    }
/*<           i2 = 1 >*/
	    i2 = 1;
/*<           i3 = i+1 >*/
	    i3 = i__ + 1;
/*<           do 380 j=i3,ibandy >*/
	    i__3 = ibandy;
	    for (j = i3; j <= i__3; ++j) {
/*<             i2 = i2+1 >*/
		++i2;
/*<             call fprota(cos,sin,h(j),ay(irot,i2)) >*/
		fprota_(&cos__, &sin__, &h__[j - 1], &ay[irot + i2 * ay_dim1])
			;
/*<  380      continue >*/
/* L380: */
	    }
/*<  390    continue >*/
L390:
	    ;
	}
/*<  400    if(nrold.eq.number) go to 420 >*/
L400:
	if (nrold == number) {
	    goto L420;
	}
/*<  410    nrold = nrold+1 >*/
L410:
	++nrold;
/*<         go to 300 >*/
	goto L300;
/*<  420  continue >*/
L420:
	;
    }
/*  backward substitution to obtain the b-spline coefficients as the */
/*  solution of the linear system    (ry) c (rx)' = h. */
/*  first step: solve the system  (ry) (c1) = h. */
/*<       k = 1 >*/
    k = 1;
/*<       do 450 i=1,nk1x >*/
    i__1 = nk1x;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         call fpback(ay,c(k),nk1y,ibandy,c(k),ny) >*/
	fpback_(&ay[ay_offset], &c__[k], &nk1y, &ibandy, &c__[k], ny);
/*<         k = k+nk1y >*/
	k += nk1y;
/*<  450  continue >*/
/* L450: */
    }
/*  second step: solve the system  c (rx)' = (c1). */
/*<       k = 0 >*/
    k = 0;
/*<       do 480 j=1,nk1y >*/
    i__1 = nk1y;
    for (j = 1; j <= i__1; ++j) {
/*<         k = k+1 >*/
	++k;
/*<         l = k >*/
	l = k;
/*<         do 460 i=1,nk1x >*/
	i__2 = nk1x;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           right(i) = c(l) >*/
	    right[i__] = c__[l];
/*<           l = l+nk1y >*/
	    l += nk1y;
/*<  460    continue >*/
/* L460: */
	}
/*<         call fpback(ax,right,nk1x,ibandx,right,nx) >*/
	fpback_(&ax[ax_offset], &right[1], &nk1x, &ibandx, &right[1], nx);
/*<         l = k >*/
	l = k;
/*<         do 470 i=1,nk1x >*/
	i__2 = nk1x;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           c(l) = right(i) >*/
	    c__[l] = right[i__];
/*<           l = l+nk1y >*/
	    l += nk1y;
/*<  470    continue >*/
/* L470: */
	}
/*<  480  continue >*/
/* L480: */
    }
/*  calculate the quantities */
/*    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my */
/*    fp = sumi=1,mx(sumj=1,my(res(i,j))) */
/*    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1 */
/*                  tx(r+kx) <= x(i) <= tx(r+kx+1) */
/*    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1 */
/*                  ty(r+ky) <= y(j) <= ty(r+ky+1) */
/*<       fp = 0. >*/
    *fp = 0.;
/*<       do 490 i=1,nx >*/
    i__1 = *nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         fpx(i) = 0. >*/
	fpx[i__] = 0.;
/*<  490  continue >*/
/* L490: */
    }
/*<       do 500 i=1,ny >*/
    i__1 = *ny;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         fpy(i) = 0. >*/
	fpy[i__] = 0.;
/*<  500  continue >*/
/* L500: */
    }
/*<       nk1y = ny-ky1 >*/
    nk1y = *ny - *ky1;
/*<       iz = 0 >*/
    iz = 0;
/*<       nroldx = 0 >*/
    nroldx = 0;
/*  main loop for the different grid points. */
/*<       do 550 i1=1,mx >*/
    i__1 = *mx;
    for (i1 = 1; i1 <= i__1; ++i1) {
/*<         numx = nrx(i1) >*/
	numx = nrx[i1];
/*<         numx1 = numx+1 >*/
	numx1 = numx + 1;
/*<         nroldy = 0 >*/
	nroldy = 0;
/*<         do 540 i2=1,my >*/
	i__2 = *my;
	for (i2 = 1; i2 <= i__2; ++i2) {
/*<           numy = nry(i2) >*/
	    numy = nry[i2];
/*<           numy1 = numy+1 >*/
	    numy1 = numy + 1;
/*<           iz = iz+1 >*/
	    ++iz;
/*  evaluate s(x,y) at the current grid point by making the sum of the */
/*  cross products of the non-zero b-splines at (x,y), multiplied with */
/*  the appropiate b-spline coefficients. */
/*<           term = 0. >*/
	    term = 0.;
/*<           k1 = numx*nk1y+numy >*/
	    k1 = numx * nk1y + numy;
/*<           do 520 l1=1,kx1 >*/
	    i__3 = *kx1;
	    for (l1 = 1; l1 <= i__3; ++l1) {
/*<             k2 = k1 >*/
		k2 = k1;
/*<             fac = spx(i1,l1) >*/
		fac = spx[i1 + l1 * spx_dim1];
/*<             do 510 l2=1,ky1 >*/
		i__4 = *ky1;
		for (l2 = 1; l2 <= i__4; ++l2) {
/*<               k2 = k2+1 >*/
		    ++k2;
/*<               term = term+fac*spy(i2,l2)*c(k2) >*/
		    term += fac * spy[i2 + l2 * spy_dim1] * c__[k2];
/*<  510        continue >*/
/* L510: */
		}
/*<             k1 = k1+nk1y >*/
		k1 += nk1y;
/*<  520      continue >*/
/* L520: */
	    }
/*  calculate the squared residual at the current grid point. */
/*<           term = (z(iz)-term)**2 >*/
/* Computing 2nd power */
	    d__1 = z__[iz] - term;
	    term = d__1 * d__1;
/*  adjust the different parameters. */
/*<           fp = fp+term >*/
	    *fp += term;
/*<           fpx(numx1) = fpx(numx1)+term >*/
	    fpx[numx1] += term;
/*<           fpy(numy1) = fpy(numy1)+term >*/
	    fpy[numy1] += term;
/*<           fac = term*half >*/
	    fac = term * half;
/*<           if(numy.eq.nroldy) go to 530 >*/
	    if (numy == nroldy) {
		goto L530;
	    }
/*<           fpy(numy1) = fpy(numy1)-fac >*/
	    fpy[numy1] -= fac;
/*<           fpy(numy) = fpy(numy)+fac >*/
	    fpy[numy] += fac;
/*<  530      nroldy = numy >*/
L530:
	    nroldy = numy;
/*<           if(numx.eq.nroldx) go to 540 >*/
	    if (numx == nroldx) {
		goto L540;
	    }
/*<           fpx(numx1) = fpx(numx1)-fac >*/
	    fpx[numx1] -= fac;
/*<           fpx(numx) = fpx(numx)+fac >*/
	    fpx[numx] += fac;
/*<  540    continue >*/
L540:
	    ;
	}
/*<         nroldx = numx >*/
	nroldx = numx;
/*<  550  continue >*/
/* L550: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpgrre_ */

