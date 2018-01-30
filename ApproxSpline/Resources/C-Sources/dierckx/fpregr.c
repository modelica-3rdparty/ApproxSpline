/* fpregr.f -- translated by f2c (version 20061008).
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

static integer c__1 = 1;

/*<    >*/
/* Subroutine */ int fpregr_(integer *iopt, doublereal *x, integer *mx, 
	doublereal *y, integer *my, doublereal *z__, integer *mz, doublereal *
	xb, doublereal *xe, doublereal *yb, doublereal *ye, integer *kx, 
	integer *ky, doublereal *s, integer *nxest, integer *nyest, 
	doublereal *tol, integer *maxit, integer *nc, integer *nx, doublereal 
	*tx, integer *ny, doublereal *ty, doublereal *c__, doublereal *fp, 
	doublereal *fp0, doublereal *fpold, doublereal *reducx, doublereal *
	reducy, doublereal *fpintx, doublereal *fpinty, integer *lastdi, 
	integer *nplusx, integer *nplusy, integer *nrx, integer *nry, integer 
	*nrdatx, integer *nrdaty, doublereal *wrk, integer *lwrk, integer *
	ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, l;
    static doublereal p, f1, f2, f3;
    static integer k3;
    static doublereal p1, p2, p3;
    static integer mm, lq;
    static doublereal rn;
    static integer mk1, kx1, kx2, ky1, ky2;
    static doublereal acc, one;
    static integer lax, lay, lbx, lby, lri, mpm, nxe, nye, nxk, lsx, lsy, 
	    ich1, ich3;
    static doublereal con1, con4, con9;
    static integer npl1, nk1x, nk1y;
    static doublereal half;
    static integer ncof, ifbx, ifby, iter;
    static doublereal fpms;
    static integer ifsx, ifsy, nplx, nply, mynx, nminx, nminy, nmaxx, nmaxy;
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpgrre_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fpknot_(doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static integer nrintx, nrinty;

/*  .. */
/*  ..scalar arguments.. */
/*<       real xb,xe,yb,ye,s,tol,fp,fp0,fpold,reducx,reducy >*/
/*<    >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my) >*/
/*  ..local scalars */
/*<       real acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,half,con1,con9,con4 >*/
/*<    >*/
/*  ..function references.. */
/*<       real abs,fprati >*/
/*<       integer max0,min0 >*/
/*  ..subroutine references.. */
/*    fpgrre,fpknot */
/*  .. */
/*   set constants */
/*<       one = 1 >*/
    /* Parameter adjustments */
    --nrx;
    --x;
    --nry;
    --y;
    --z__;
    --nrdatx;
    --fpintx;
    --tx;
    --nrdaty;
    --fpinty;
    --ty;
    --c__;
    --wrk;

    /* Function Body */
    one = 1.;
/*<       half = 0.5e0 >*/
    half = .5;
/*<       con1 = 0.1e0 >*/
    con1 = .1;
/*<       con9 = 0.9e0 >*/
    con9 = .9;
/*<       con4 = 0.4e-01 >*/
    con4 = .04;
/*  we partition the working space. */
/*<       kx1 = kx+1 >*/
    kx1 = *kx + 1;
/*<       ky1 = ky+1 >*/
    ky1 = *ky + 1;
/*<       kx2 = kx1+1 >*/
    kx2 = kx1 + 1;
/*<       ky2 = ky1+1 >*/
    ky2 = ky1 + 1;
/*<       lsx = 1 >*/
    lsx = 1;
/*<       lsy = lsx+mx*kx1 >*/
    lsy = lsx + *mx * kx1;
/*<       lri = lsy+my*ky1 >*/
    lri = lsy + *my * ky1;
/*<       mm = max0(nxest,my) >*/
    mm = max(*nxest,*my);
/*<       lq = lri+mm >*/
    lq = lri + mm;
/*<       mynx = nxest*my >*/
    mynx = *nxest * *my;
/*<       lax = lq+mynx >*/
    lax = lq + mynx;
/*<       nxk = nxest*kx2 >*/
    nxk = *nxest * kx2;
/*<       lbx = lax+nxk >*/
    lbx = lax + nxk;
/*<       lay = lbx+nxk >*/
    lay = lbx + nxk;
/*<       lby = lay+nyest*ky2 >*/
    lby = lay + *nyest * ky2;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 1: determination of the number of knots and their position.     c */
/* ****************************************************************     c */
/*  given a set of knots we compute the least-squares spline sinf(x,y), c */
/*  and the corresponding sum of squared residuals fp=f(p=inf).         c */
/*  if iopt=-1  sinf(x,y) is the requested approximation.               c */
/*  if iopt=0 or iopt=1 we check whether we can accept the knots:       c */
/*    if fp <=s we will continue with the current set of knots.         c */
/*    if fp > s we will increase the number of knots and compute the    c */
/*       corresponding least-squares spline until finally fp<=s.        c */
/*    the initial choice of knots depends on the value of s and iopt.   c */
/*    if s=0 we have spline interpolation; in that case the number of   c */
/*    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c */
/*    if s>0 and                                                        c */
/*     *iopt=0 we first compute the least-squares polynomial of degree  c */
/*      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c */
/*     *iopt=1 we start with the knots found at the last call of the    c */
/*      routine, except for the case that s > fp0; then we can compute  c */
/*      the least-squares polynomial directly.                          c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  determine the number of knots for polynomial approximation. */
/*<       nminx = 2*kx1 >*/
    nminx = kx1 << 1;
/*<       nminy = 2*ky1 >*/
    nminy = ky1 << 1;
/*<       if(iopt.lt.0) go to 120 >*/
    if (*iopt < 0) {
	goto L120;
    }
/*  acc denotes the absolute tolerance for the root of f(p)=s. */
/*<       acc = tol*s >*/
    acc = *tol * *s;
/*  find nmaxx and nmaxy which denote the number of knots in x- and y- */
/*  direction in case of spline interpolation. */
/*<       nmaxx = mx+kx1 >*/
    nmaxx = *mx + kx1;
/*<       nmaxy = my+ky1 >*/
    nmaxy = *my + ky1;
/*  find nxe and nye which denote the maximum number of knots */
/*  allowed in each direction */
/*<       nxe = min0(nmaxx,nxest) >*/
    nxe = min(nmaxx,*nxest);
/*<       nye = min0(nmaxy,nyest) >*/
    nye = min(nmaxy,*nyest);
/*<       if(s.gt.0.) go to 100 >*/
    if (*s > 0.) {
	goto L100;
    }
/*  if s = 0, s(x,y) is an interpolating spline. */
/*<       nx = nmaxx >*/
    *nx = nmaxx;
/*<       ny = nmaxy >*/
    *ny = nmaxy;
/*  test whether the required storage space exceeds the available one. */
/*<       if(ny.gt.nyest .or. nx.gt.nxest) go to 420 >*/
    if (*ny > *nyest || *nx > *nxest) {
	goto L420;
    }
/*  find the position of the interior knots in case of interpolation. */
/*  the knots in the x-direction. */
/*<       mk1 = mx-kx1 >*/
    mk1 = *mx - kx1;
/*<       if(mk1.eq.0) go to 60 >*/
    if (mk1 == 0) {
	goto L60;
    }
/*<       k3 = kx/2 >*/
    k3 = *kx / 2;
/*<       i = kx1+1 >*/
    i__ = kx1 + 1;
/*<       j = k3+2 >*/
    j = k3 + 2;
/*<       if(k3*2.eq.kx) go to 40 >*/
    if (k3 << 1 == *kx) {
	goto L40;
    }
/*<       do 30 l=1,mk1 >*/
    i__1 = mk1;
    for (l = 1; l <= i__1; ++l) {
/*<         tx(i) = x(j) >*/
	tx[i__] = x[j];
/*<         i = i+1 >*/
	++i__;
/*<         j = j+1 >*/
	++j;
/*<   30  continue >*/
/* L30: */
    }
/*<       go to 60 >*/
    goto L60;
/*<   40  do 50 l=1,mk1 >*/
L40:
    i__1 = mk1;
    for (l = 1; l <= i__1; ++l) {
/*<         tx(i) = (x(j)+x(j-1))*half >*/
	tx[i__] = (x[j] + x[j - 1]) * half;
/*<         i = i+1 >*/
	++i__;
/*<         j = j+1 >*/
	++j;
/*<   50  continue >*/
/* L50: */
    }
/*  the knots in the y-direction. */
/*<   60  mk1 = my-ky1 >*/
L60:
    mk1 = *my - ky1;
/*<       if(mk1.eq.0) go to 120 >*/
    if (mk1 == 0) {
	goto L120;
    }
/*<       k3 = ky/2 >*/
    k3 = *ky / 2;
/*<       i = ky1+1 >*/
    i__ = ky1 + 1;
/*<       j = k3+2 >*/
    j = k3 + 2;
/*<       if(k3*2.eq.ky) go to 80 >*/
    if (k3 << 1 == *ky) {
	goto L80;
    }
/*<       do 70 l=1,mk1 >*/
    i__1 = mk1;
    for (l = 1; l <= i__1; ++l) {
/*<         ty(i) = y(j) >*/
	ty[i__] = y[j];
/*<         i = i+1 >*/
	++i__;
/*<         j = j+1 >*/
	++j;
/*<   70  continue >*/
/* L70: */
    }
/*<       go to 120 >*/
    goto L120;
/*<   80  do 90 l=1,mk1 >*/
L80:
    i__1 = mk1;
    for (l = 1; l <= i__1; ++l) {
/*<         ty(i) = (y(j)+y(j-1))*half >*/
	ty[i__] = (y[j] + y[j - 1]) * half;
/*<         i = i+1 >*/
	++i__;
/*<         j = j+1 >*/
	++j;
/*<   90  continue >*/
/* L90: */
    }
/*<       go to 120 >*/
    goto L120;
/*  if s > 0 our initial choice of knots depends on the value of iopt. */
/*<  100  if(iopt.eq.0) go to 115 >*/
L100:
    if (*iopt == 0) {
	goto L115;
    }
/*<       if(fp0.le.s) go to 115 >*/
    if (*fp0 <= *s) {
	goto L115;
    }
/*  if iopt=1 and fp0 > s we start computing the least- squares spline */
/*  according to the set of knots found at the last call of the routine. */
/*  we determine the number of grid coordinates x(i) inside each knot */
/*  interval (tx(l),tx(l+1)). */
/*<       l = kx2 >*/
    l = kx2;
/*<       j = 1 >*/
    j = 1;
/*<       nrdatx(1) = 0 >*/
    nrdatx[1] = 0;
/*<       mpm = mx-1 >*/
    mpm = *mx - 1;
/*<       do 105 i=2,mpm >*/
    i__1 = mpm;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         nrdatx(j) = nrdatx(j)+1 >*/
	++nrdatx[j];
/*<         if(x(i).lt.tx(l)) go to 105 >*/
	if (x[i__] < tx[l]) {
	    goto L105;
	}
/*<         nrdatx(j) = nrdatx(j)-1 >*/
	--nrdatx[j];
/*<         l = l+1 >*/
	++l;
/*<         j = j+1 >*/
	++j;
/*<         nrdatx(j) = 0 >*/
	nrdatx[j] = 0;
/*<  105  continue >*/
L105:
	;
    }
/*  we determine the number of grid coordinates y(i) inside each knot */
/*  interval (ty(l),ty(l+1)). */
/*<       l = ky2 >*/
    l = ky2;
/*<       j = 1 >*/
    j = 1;
/*<       nrdaty(1) = 0 >*/
    nrdaty[1] = 0;
/*<       mpm = my-1 >*/
    mpm = *my - 1;
/*<       do 110 i=2,mpm >*/
    i__1 = mpm;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         nrdaty(j) = nrdaty(j)+1 >*/
	++nrdaty[j];
/*<         if(y(i).lt.ty(l)) go to 110 >*/
	if (y[i__] < ty[l]) {
	    goto L110;
	}
/*<         nrdaty(j) = nrdaty(j)-1 >*/
	--nrdaty[j];
/*<         l = l+1 >*/
	++l;
/*<         j = j+1 >*/
	++j;
/*<         nrdaty(j) = 0 >*/
	nrdaty[j] = 0;
/*<  110  continue >*/
L110:
	;
    }
/*<       go to 120 >*/
    goto L120;
/*  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares */
/*  polynomial of degree kx in x and ky in y (which is a spline without */
/*  interior knots). */
/*<  115  nx = nminx >*/
L115:
    *nx = nminx;
/*<       ny = nminy >*/
    *ny = nminy;
/*<       nrdatx(1) = mx-2 >*/
    nrdatx[1] = *mx - 2;
/*<       nrdaty(1) = my-2 >*/
    nrdaty[1] = *my - 2;
/*<       lastdi = 0 >*/
    *lastdi = 0;
/*<       nplusx = 0 >*/
    *nplusx = 0;
/*<       nplusy = 0 >*/
    *nplusy = 0;
/*<       fp0 = 0. >*/
    *fp0 = 0.;
/*<       fpold = 0. >*/
    *fpold = 0.;
/*<       reducx = 0. >*/
    *reducx = 0.;
/*<       reducy = 0. >*/
    *reducy = 0.;
/*<  120  mpm = mx+my >*/
L120:
    mpm = *mx + *my;
/*<       ifsx = 0 >*/
    ifsx = 0;
/*<       ifsy = 0 >*/
    ifsy = 0;
/*<       ifbx = 0 >*/
    ifbx = 0;
/*<       ifby = 0 >*/
    ifby = 0;
/*<       p = -one >*/
    p = -one;
/*  main loop for the different sets of knots.mpm=mx+my is a save upper */
/*  bound for the number of trials. */
/*<       do 250 iter=1,mpm >*/
    i__1 = mpm;
    for (iter = 1; iter <= i__1; ++iter) {
/*<         if(nx.eq.nminx .and. ny.eq.nminy) ier = -2 >*/
	if (*nx == nminx && *ny == nminy) {
	    *ier = -2;
	}
/*  find nrintx (nrinty) which is the number of knot intervals in the */
/*  x-direction (y-direction). */
/*<         nrintx = nx-nminx+1 >*/
	nrintx = *nx - nminx + 1;
/*<         nrinty = ny-nminy+1 >*/
	nrinty = *ny - nminy + 1;
/*  find ncof, the number of b-spline coefficients for the current set */
/*  of knots. */
/*<         nk1x = nx-kx1 >*/
	nk1x = *nx - kx1;
/*<         nk1y = ny-ky1 >*/
	nk1y = *ny - ky1;
/*<         ncof = nk1x*nk1y >*/
	ncof = nk1x * nk1y;
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(x,y). */
/*<         i = nx >*/
	i__ = *nx;
/*<         do 130 j=1,kx1 >*/
	i__2 = kx1;
	for (j = 1; j <= i__2; ++j) {
/*<           tx(j) = xb >*/
	    tx[j] = *xb;
/*<           tx(i) = xe >*/
	    tx[i__] = *xe;
/*<           i = i-1 >*/
	    --i__;
/*<  130    continue >*/
/* L130: */
	}
/*<         i = ny >*/
	i__ = *ny;
/*<         do 140 j=1,ky1 >*/
	i__2 = ky1;
	for (j = 1; j <= i__2; ++j) {
/*<           ty(j) = yb >*/
	    ty[j] = *yb;
/*<           ty(i) = ye >*/
	    ty[i__] = *ye;
/*<           i = i-1 >*/
	    --i__;
/*<  140    continue >*/
/* L140: */
	}
/*  find the least-squares spline sinf(x,y) and calculate for each knot */
/*  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum */
/*  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2, */
/*  ...,ny-2*ky-1) for the data points having their absciss (ordinate)- */
/*  value belonging to that interval. */
/*  fp gives the total sum of squared residuals. */
/*<    >*/
	fpgrre_(&ifsx, &ifsy, &ifbx, &ifby, &x[1], mx, &y[1], my, &z__[1], mz,
		 kx, ky, &tx[1], nx, &ty[1], ny, &p, &c__[1], nc, fp, &fpintx[
		1], &fpinty[1], &mm, &mynx, &kx1, &kx2, &ky1, &ky2, &wrk[lsx],
		 &wrk[lsy], &wrk[lri], &wrk[lq], &wrk[lax], &wrk[lay], &wrk[
		lbx], &wrk[lby], &nrx[1], &nry[1]);
/*<         if(ier.eq.(-2)) fp0 = fp >*/
	if (*ier == -2) {
	    *fp0 = *fp;
	}
/*  test whether the least-squares spline is an acceptable solution. */
/*<         if(iopt.lt.0) go to 440 >*/
	if (*iopt < 0) {
	    goto L440;
	}
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms) .lt. acc) go to 440 >*/
	if (abs(fpms) < acc) {
	    goto L440;
	}
/*  if f(p=inf) < s, we accept the choice of knots. */
/*<         if(fpms.lt.0.) go to 300 >*/
	if (fpms < 0.) {
	    goto L300;
	}
/*  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline. */
/*<         if(nx.eq.nmaxx .and. ny.eq.nmaxy) go to 430 >*/
	if (*nx == nmaxx && *ny == nmaxy) {
	    goto L430;
	}
/*  increase the number of knots. */
/*  if nx=nxe and ny=nye we cannot further increase the number of knots */
/*  because of the storage capacity limitation. */
/*<         if(nx.eq.nxe .and. ny.eq.nye) go to 420 >*/
	if (*nx == nxe && *ny == nye) {
	    goto L420;
	}
/*<         ier = 0 >*/
	*ier = 0;
/*  adjust the parameter reducx or reducy according to the direction */
/*  in which the last added knots were located. */
/*<         if(lastdi) 150,170,160 >*/
	if (*lastdi < 0) {
	    goto L150;
	} else if (*lastdi == 0) {
	    goto L170;
	} else {
	    goto L160;
	}
/*<  150    reducx = fpold-fp >*/
L150:
	*reducx = *fpold - *fp;
/*<         go to 170 >*/
	goto L170;
/*<  160    reducy = fpold-fp >*/
L160:
	*reducy = *fpold - *fp;
/*  store the sum of squared residuals for the current set of knots. */
/*<  170    fpold = fp >*/
L170:
	*fpold = *fp;
/*  find nplx, the number of knots we should add in the x-direction. */
/*<         nplx = 1 >*/
	nplx = 1;
/*<         if(nx.eq.nminx) go to 180 >*/
	if (*nx == nminx) {
	    goto L180;
	}
/*<         npl1 = nplusx*2 >*/
	npl1 = *nplusx << 1;
/*<         rn = nplusx >*/
	rn = (doublereal) (*nplusx);
/*<         if(reducx.gt.acc) npl1 = rn*fpms/reducx >*/
	if (*reducx > acc) {
	    npl1 = (integer) (rn * fpms / *reducx);
	}
/*<         nplx = min0(nplusx*2,max0(npl1,nplusx/2,1)) >*/
/* Computing MIN */
/* Computing MAX */
	i__4 = npl1, i__5 = *nplusx / 2, i__4 = max(i__4,i__5);
	i__2 = *nplusx << 1, i__3 = max(i__4,1);
	nplx = min(i__2,i__3);
/*  find nply, the number of knots we should add in the y-direction. */
/*<  180    nply = 1 >*/
L180:
	nply = 1;
/*<         if(ny.eq.nminy) go to 190 >*/
	if (*ny == nminy) {
	    goto L190;
	}
/*<         npl1 = nplusy*2 >*/
	npl1 = *nplusy << 1;
/*<         rn = nplusy >*/
	rn = (doublereal) (*nplusy);
/*<         if(reducy.gt.acc) npl1 = rn*fpms/reducy >*/
	if (*reducy > acc) {
	    npl1 = (integer) (rn * fpms / *reducy);
	}
/*<         nply = min0(nplusy*2,max0(npl1,nplusy/2,1)) >*/
/* Computing MIN */
/* Computing MAX */
	i__4 = npl1, i__5 = *nplusy / 2, i__4 = max(i__4,i__5);
	i__2 = *nplusy << 1, i__3 = max(i__4,1);
	nply = min(i__2,i__3);
/*<  190    if(nplx-nply) 210,200,230 >*/
L190:
	if ((i__2 = nplx - nply) < 0) {
	    goto L210;
	} else if (i__2 == 0) {
	    goto L200;
	} else {
	    goto L230;
	}
/*<  200    if(lastdi.lt.0) go to 230 >*/
L200:
	if (*lastdi < 0) {
	    goto L230;
	}
/*<  210    if(nx.eq.nxe) go to 230 >*/
L210:
	if (*nx == nxe) {
	    goto L230;
	}
/*  addition in the x-direction. */
/*<         lastdi = -1 >*/
	*lastdi = -1;
/*<         nplusx = nplx >*/
	*nplusx = nplx;
/*<         ifsx = 0 >*/
	ifsx = 0;
/*<         do 220 l=1,nplusx >*/
	i__2 = *nplusx;
	for (l = 1; l <= i__2; ++l) {
/*  add a new knot in the x-direction */
/*<           call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1) >*/
	    fpknot_(&x[1], mx, &tx[1], nx, &fpintx[1], &nrdatx[1], &nrintx, 
		    nxest, &c__1);
/*  test whether we cannot further increase the number of knots in the */
/*  x-direction. */
/*<           if(nx.eq.nxe) go to 250 >*/
	    if (*nx == nxe) {
		goto L250;
	    }
/*<  220    continue >*/
/* L220: */
	}
/*<         go to 250 >*/
	goto L250;
/*<  230    if(ny.eq.nye) go to 210 >*/
L230:
	if (*ny == nye) {
	    goto L210;
	}
/*  addition in the y-direction. */
/*<         lastdi = 1 >*/
	*lastdi = 1;
/*<         nplusy = nply >*/
	*nplusy = nply;
/*<         ifsy = 0 >*/
	ifsy = 0;
/*<         do 240 l=1,nplusy >*/
	i__2 = *nplusy;
	for (l = 1; l <= i__2; ++l) {
/*  add a new knot in the y-direction. */
/*<           call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1) >*/
	    fpknot_(&y[1], my, &ty[1], ny, &fpinty[1], &nrdaty[1], &nrinty, 
		    nyest, &c__1);
/*  test whether we cannot further increase the number of knots in the */
/*  y-direction. */
/*<           if(ny.eq.nye) go to 250 >*/
	    if (*ny == nye) {
		goto L250;
	    }
/*<  240    continue >*/
/* L240: */
	}
/*  restart the computations with the new set of knots. */
/*<  250  continue >*/
L250:
	;
    }
/*  test whether the least-squares polynomial is a solution of our */
/*  approximation problem. */
/*<  300  if(ier.eq.(-2)) go to 440 >*/
L300:
    if (*ier == -2) {
	goto L440;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing spline sp(x,y)                c */
/* *****************************************************                c */
/*  we have determined the number of knots and their position. we now   c */
/*  compute the b-spline coefficients of the smoothing spline sp(x,y).  c */
/*  this smoothing spline varies with the parameter p in such a way thatc */
/*    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c */
/*  is a continuous, strictly decreasing function of p. moreover the    c */
/*  least-squares polynomial corresponds to p=0 and the least-squares   c */
/*  spline to p=infinity. iteratively we then have to determine the     c */
/*  positive value of p such that f(p)=s. the process which is proposed c */
/*  here makes use of rational interpolation. f(p) is approximated by a c */
/*  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c */
/*  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c */
/*  are used to calculate the new value of p such that r(p)=s.          c */
/*  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  initial value for p. */
/*<       p1 = 0. >*/
    p1 = 0.;
/*<       f1 = fp0-s >*/
    f1 = *fp0 - *s;
/*<       p3 = -one >*/
    p3 = -one;
/*<       f3 = fpms >*/
    f3 = fpms;
/*<       p = one >*/
    p = one;
/*<       ich1 = 0 >*/
    ich1 = 0;
/*<       ich3 = 0 >*/
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
/*<       do 350 iter = 1,maxit >*/
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the smoothing spline sp(x,y) and the corresponding sum of */
/*  squared residuals fp. */
/*<    >*/
	fpgrre_(&ifsx, &ifsy, &ifbx, &ifby, &x[1], mx, &y[1], my, &z__[1], mz,
		 kx, ky, &tx[1], nx, &ty[1], ny, &p, &c__[1], nc, fp, &fpintx[
		1], &fpinty[1], &mm, &mynx, &kx1, &kx2, &ky1, &ky2, &wrk[lsx],
		 &wrk[lsy], &wrk[lri], &wrk[lq], &wrk[lax], &wrk[lay], &wrk[
		lbx], &wrk[lby], &nrx[1], &nry[1]);
/*  test whether the approximation sp(x,y) is an acceptable solution. */
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).lt.acc) go to 440 >*/
	if (abs(fpms) < acc) {
	    goto L440;
	}
/*  test whether the maximum allowable number of iterations has been */
/*  reached. */
/*<         if(iter.eq.maxit) go to 400 >*/
	if (iter == *maxit) {
	    goto L400;
	}
/*  carry out one more step of the iteration process. */
/*<         p2 = p >*/
	p2 = p;
/*<         f2 = fpms >*/
	f2 = fpms;
/*<         if(ich3.ne.0) go to 320 >*/
	if (ich3 != 0) {
	    goto L320;
	}
/*<         if((f2-f3).gt.acc) go to 310 >*/
	if (f2 - f3 > acc) {
	    goto L310;
	}
/*  our initial choice of p is too large. */
/*<         p3 = p2 >*/
	p3 = p2;
/*<         f3 = f2 >*/
	f3 = f2;
/*<         p = p*con4 >*/
	p *= con4;
/*<         if(p.le.p1) p = p1*con9 + p2*con1 >*/
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
/*<         go to 350 >*/
	goto L350;
/*<  310    if(f2.lt.0.) ich3 = 1 >*/
L310:
	if (f2 < 0.) {
	    ich3 = 1;
	}
/*<  320    if(ich1.ne.0) go to 340 >*/
L320:
	if (ich1 != 0) {
	    goto L340;
	}
/*<         if((f1-f2).gt.acc) go to 330 >*/
	if (f1 - f2 > acc) {
	    goto L330;
	}
/*  our initial choice of p is too small */
/*<         p1 = p2 >*/
	p1 = p2;
/*<         f1 = f2 >*/
	f1 = f2;
/*<         p = p/con4 >*/
	p /= con4;
/*<         if(p3.lt.0.) go to 350 >*/
	if (p3 < 0.) {
	    goto L350;
	}
/*<         if(p.ge.p3) p = p2*con1 + p3*con9 >*/
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
/*<         go to 350 >*/
	goto L350;
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
/*<  330    if(f2.gt.0.) ich1 = 1 >*/
L330:
	if (f2 > 0.) {
	    ich1 = 1;
	}
/*<  340    if(f2.ge.f1 .or. f2.le.f3) go to 410 >*/
L340:
	if (f2 >= f1 || f2 <= f3) {
	    goto L410;
	}
/*  find the new value of p. */
/*<         p = fprati(p1,f1,p2,f2,p3,f3) >*/
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
/*<  350  continue >*/
L350:
	;
    }
/*  error codes and messages. */
/*<  400  ier = 3 >*/
L400:
    *ier = 3;
/*<       go to 440 >*/
    goto L440;
/*<  410  ier = 2 >*/
L410:
    *ier = 2;
/*<       go to 440 >*/
    goto L440;
/*<  420  ier = 1 >*/
L420:
    *ier = 1;
/*<       go to 440 >*/
    goto L440;
/*<  430  ier = -1 >*/
L430:
    *ier = -1;
/*<       fp = 0. >*/
    *fp = 0.;
/*<  440  return >*/
L440:
    return 0;
/*<       end >*/
} /* fpregr_ */

