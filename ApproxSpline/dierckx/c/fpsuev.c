/* fpsuev.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wu,wv,lu,lv) >*/
/* Subroutine */ int fpsuev_(integer *idim, doublereal *tu, integer *nu, 
	doublereal *tv, integer *nv, doublereal *c__, doublereal *u, integer *
	mu, doublereal *v, integer *mv, doublereal *f, doublereal *wu, 
	doublereal *wv, integer *lu, integer *lv)
{
    /* System generated locals */
    integer wu_dim1, wu_offset, wv_dim1, wv_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[4];
    static integer i__, j, k, l, m, i1, j1, l1, l2, l3;
    static doublereal tb, te, sp;
    static integer nu4, nv4;
    static doublereal arg;
    static integer nuv;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  ..scalar arguments.. */
/*<       integer idim,nu,nv,mu,mv >*/
/*  ..array arguments.. */
/*<       integer lu(mu),lv(mv) >*/
/*<    >*/
/*  ..local scalars.. */
/*<       integer i,i1,j,j1,k,l,l1,l2,l3,m,nuv,nu4,nv4 >*/
/*<       real arg,sp,tb,te >*/
/*  ..local arrays.. */
/*<       real h(4) >*/
/*  ..subroutine references.. */
/*    fpbspl */
/*  .. */
/*<       nu4 = nu-4 >*/
    /* Parameter adjustments */
    --tu;
    --c__;
    --tv;
    --lu;
    wu_dim1 = *mu;
    wu_offset = 1 + wu_dim1;
    wu -= wu_offset;
    --u;
    --lv;
    wv_dim1 = *mv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    --f;
    --v;

    /* Function Body */
    nu4 = *nu - 4;
/*<       tb = tu(4) >*/
    tb = tu[4];
/*<       te = tu(nu4+1) >*/
    te = tu[nu4 + 1];
/*<       l = 4 >*/
    l = 4;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*<       do 40 i=1,mu >*/
    i__1 = *mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
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
/*<   10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20 >*/
L10:
	if (arg < tu[l1] || l == nu4) {
	    goto L20;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 10 >*/
	goto L10;
/*<   20    call fpbspl(tu,nu,3,arg,l,h) >*/
L20:
	fpbspl_(&tu[1], nu, &c__3, &arg, &l, h__);
/*<         lu(i) = l-4 >*/
	lu[i__] = l - 4;
/*<         do 30 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           wu(i,j) = h(j) >*/
	    wu[i__ + j * wu_dim1] = h__[j - 1];
/*<   30    continue >*/
/* L30: */
	}
/*<   40  continue >*/
/* L40: */
    }
/*<       nv4 = nv-4 >*/
    nv4 = *nv - 4;
/*<       tb = tv(4) >*/
    tb = tv[4];
/*<       te = tv(nv4+1) >*/
    te = tv[nv4 + 1];
/*<       l = 4 >*/
    l = 4;
/*<       l1 = l+1 >*/
    l1 = l + 1;
/*<       do 80 i=1,mv >*/
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         arg = v(i) >*/
	arg = v[i__];
/*<         if(arg.lt.tb) arg = tb >*/
	if (arg < tb) {
	    arg = tb;
	}
/*<         if(arg.gt.te) arg = te >*/
	if (arg > te) {
	    arg = te;
	}
/*<   50    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 60 >*/
L50:
	if (arg < tv[l1] || l == nv4) {
	    goto L60;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 50 >*/
	goto L50;
/*<   60    call fpbspl(tv,nv,3,arg,l,h) >*/
L60:
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
/*<         lv(i) = l-4 >*/
	lv[i__] = l - 4;
/*<         do 70 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           wv(i,j) = h(j) >*/
	    wv[i__ + j * wv_dim1] = h__[j - 1];
/*<   70    continue >*/
/* L70: */
	}
/*<   80  continue >*/
/* L80: */
    }
/*<       m = 0 >*/
    m = 0;
/*<       nuv = nu4*nv4 >*/
    nuv = nu4 * nv4;
/*<       do 140 k=1,idim >*/
    i__1 = *idim;
    for (k = 1; k <= i__1; ++k) {
/*<         l3 = (k-1)*nuv >*/
	l3 = (k - 1) * nuv;
/*<         do 130 i=1,mu >*/
	i__2 = *mu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           l = lu(i)*nv4+l3 >*/
	    l = lu[i__] * nv4 + l3;
/*<           do 90 i1=1,4 >*/
	    for (i1 = 1; i1 <= 4; ++i1) {
/*<             h(i1) = wu(i,i1) >*/
		h__[i1 - 1] = wu[i__ + i1 * wu_dim1];
/*<   90      continue >*/
/* L90: */
	    }
/*<           do 120 j=1,mv >*/
	    i__3 = *mv;
	    for (j = 1; j <= i__3; ++j) {
/*<             l1 = l+lv(j) >*/
		l1 = l + lv[j];
/*<             sp = 0. >*/
		sp = 0.;
/*<             do 110 i1=1,4 >*/
		for (i1 = 1; i1 <= 4; ++i1) {
/*<               l2 = l1 >*/
		    l2 = l1;
/*<               do 100 j1=1,4 >*/
		    for (j1 = 1; j1 <= 4; ++j1) {
/*<                 l2 = l2+1 >*/
			++l2;
/*<                 sp = sp+c(l2)*h(i1)*wv(j,j1) >*/
			sp += c__[l2] * h__[i1 - 1] * wv[j + j1 * wv_dim1];
/*<  100          continue >*/
/* L100: */
		    }
/*<               l1 = l1+nv4 >*/
		    l1 += nv4;
/*<  110        continue >*/
/* L110: */
		}
/*<             m = m+1 >*/
		++m;
/*<             f(m) = sp >*/
		f[m] = sp;
/*<  120      continue >*/
/* L120: */
	    }
/*<  130    continue >*/
/* L130: */
	}
/*<  140  continue >*/
/* L140: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpsuev_ */

