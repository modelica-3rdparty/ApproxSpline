/* fpbfout.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpbfou(t,n,par,ress,resc) >*/
/* Subroutine */ int fpbfou_(doublereal *t, integer *n, doublereal *par, 
	doublereal *ress, doublereal *resc)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal c1, c2, f1, f2, f3, s1, s2, hc[5];
    static integer ic;
    static doublereal ak, co[5];
    static integer jj, li, lj;
    static doublereal rc[3];
    static integer ll;
    static doublereal hs[5];
    static integer is;
    static doublereal si[5], rs[3];
    static integer jp1, jp4, nm3, nm7;
    static doublereal fac, one;
    static integer ipj, nmj;
    static doublereal eps, six, con1, con2, beta, sign, term, delta, quart;
    extern /* Subroutine */ int fpcsin_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/*  subroutine fpbfou calculates the integrals */
/*                    /t(n-3) */
/*    ress(j) =      !        nj,4(x)*sin(par*x) dx    and */
/*              t(4)/ */
/*                    /t(n-3) */
/*    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4 */
/*              t(4)/ */
/*  where nj,4(x) denotes the cubic b-spline defined on the knots */
/*  t(j),t(j+1),...,t(j+4). */

/*  calling sequence: */
/*     call fpbfou(t,n,par,ress,resc) */

/*  input parameters: */
/*    t    : real array,length n, containing the knots. */
/*    n    : integer, containing the number of knots. */
/*    par  : real, containing the value of the parameter par. */

/*  output parameters: */
/*    ress  : real array,length n, containing the integrals ress(j). */
/*    resc  : real array,length n, containing the integrals resc(j). */

/*  restrictions: */
/*    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3). */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer n >*/
/*<       real par >*/
/*  ..array arguments.. */
/*<       real t(n),ress(n),resc(n) >*/
/*  ..local scalars.. */
/*<       integer i,ic,ipj,is,j,jj,jp1,jp4,k,li,lj,ll,nmj,nm3,nm7 >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real co(5),si(5),hs(5),hc(5),rs(3),rc(3) >*/
/*  ..function references.. */
/*<       real cos,sin,abs >*/
/*  .. */
/*  initialization. */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --resc;
    --ress;
    --t;

    /* Function Body */
    one = 1.;
/*<       six = 0.6e+01 >*/
    six = 6.;
/*<       eps = 0.1e-07 >*/
    eps = 1e-8;
/*<       quart = 0.25e0 >*/
    quart = .25;
/*<       con1 = 0.5e-01 >*/
    con1 = .05;
/*<       con2 = 0.12e+03 >*/
    con2 = 120.;
/*<       nm3 = n-3 >*/
    nm3 = *n - 3;
/*<       nm7 = n-7 >*/
    nm7 = *n - 7;
/*<       if(par.ne.0.) term = six/par >*/
    if (*par != 0.) {
	term = six / *par;
    }
/*<       beta = par*t(4) >*/
    beta = *par * t[4];
/*<       co(1) = cos(beta) >*/
    co[0] = cos(beta);
/*<       si(1) = sin(beta) >*/
    si[0] = sin(beta);
/*  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up */
/*  a divided difference table. */
/*<       do 30 j=1,3 >*/
    for (j = 1; j <= 3; ++j) {
/*<         jp1 = j+1 >*/
	jp1 = j + 1;
/*<         jp4 = j+4 >*/
	jp4 = j + 4;
/*<         beta = par*t(jp4) >*/
	beta = *par * t[jp4];
/*<         co(jp1) = cos(beta) >*/
	co[jp1 - 1] = cos(beta);
/*<         si(jp1) = sin(beta) >*/
	si[jp1 - 1] = sin(beta);
/*<    >*/
	fpcsin_(&t[4], &t[jp4], par, si, co, &si[jp1 - 1], &co[jp1 - 1], &rs[
		j - 1], &rc[j - 1]);
/*<         i = 5-j >*/
	i__ = 5 - j;
/*<         hs(i) = 0. >*/
	hs[i__ - 1] = 0.;
/*<         hc(i) = 0. >*/
	hc[i__ - 1] = 0.;
/*<         do 10 jj=1,j >*/
	i__1 = j;
	for (jj = 1; jj <= i__1; ++jj) {
/*<           ipj = i+jj >*/
	    ipj = i__ + jj;
/*<           hs(ipj) = rs(jj) >*/
	    hs[ipj - 1] = rs[jj - 1];
/*<           hc(ipj) = rc(jj) >*/
	    hc[ipj - 1] = rc[jj - 1];
/*<   10    continue >*/
/* L10: */
	}
/*<         do 20 jj=1,3 >*/
	for (jj = 1; jj <= 3; ++jj) {
/*<           if(i.lt.jj) i = jj >*/
	    if (i__ < jj) {
		i__ = jj;
	    }
/*<           k = 5 >*/
	    k = 5;
/*<           li = jp4 >*/
	    li = jp4;
/*<           do 20 ll=i,4 >*/
	    for (ll = i__; ll <= 4; ++ll) {
/*<             lj = li-jj >*/
		lj = li - jj;
/*<             fac = t(li)-t(lj) >*/
		fac = t[li] - t[lj];
/*<             hs(k) = (hs(k)-hs(k-1))/fac >*/
		hs[k - 1] = (hs[k - 1] - hs[k - 2]) / fac;
/*<             hc(k) = (hc(k)-hc(k-1))/fac >*/
		hc[k - 1] = (hc[k - 1] - hc[k - 2]) / fac;
/*<             k = k-1 >*/
		--k;
/*<             li = li-1 >*/
		--li;
/*<   20    continue >*/
/* L20: */
	    }
	}
/*<         ress(j) = hs(5)-hs(4) >*/
	ress[j] = hs[4] - hs[3];
/*<         resc(j) = hc(5)-hc(4) >*/
	resc[j] = hc[4] - hc[3];
/*<   30  continue >*/
/* L30: */
    }
/*<       if(nm7.lt.4) go to 160 >*/
    if (nm7 < 4) {
	goto L160;
    }
/*  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7. */
/*<       do 150 j=4,nm7 >*/
    i__1 = nm7;
    for (j = 4; j <= i__1; ++j) {
/*<         jp4 = j+4 >*/
	jp4 = j + 4;
/*<         beta = par*t(jp4) >*/
	beta = *par * t[jp4];
/*<         co(5) = cos(beta) >*/
	co[4] = cos(beta);
/*<         si(5) = sin(beta) >*/
	si[4] = sin(beta);
/*<         delta = t(jp4)-t(j) >*/
	delta = t[jp4] - t[j];
/*  the way of computing ress(j) and resc(j) depends on the value of */
/*  beta = par*(t(j+4)-t(j)). */
/*<         beta = delta*par >*/
	beta = delta * *par;
/*<         if(abs(beta).le.one) go to 60 >*/
	if (abs(beta) <= one) {
	    goto L60;
	}
/*  if !beta! > 1 the integrals are calculated by setting up a divided */
/*  difference table. */
/*<         do 40 k=1,5 >*/
	for (k = 1; k <= 5; ++k) {
/*<           hs(k) = si(k) >*/
	    hs[k - 1] = si[k - 1];
/*<           hc(k) = co(k) >*/
	    hc[k - 1] = co[k - 1];
/*<   40    continue >*/
/* L40: */
	}
/*<         do 50 jj=1,3 >*/
	for (jj = 1; jj <= 3; ++jj) {
/*<           k = 5 >*/
	    k = 5;
/*<           li = jp4 >*/
	    li = jp4;
/*<           do 50 ll=jj,4 >*/
	    for (ll = jj; ll <= 4; ++ll) {
/*<             lj = li-jj >*/
		lj = li - jj;
/*<             fac = par*(t(li)-t(lj)) >*/
		fac = *par * (t[li] - t[lj]);
/*<             hs(k) = (hs(k)-hs(k-1))/fac >*/
		hs[k - 1] = (hs[k - 1] - hs[k - 2]) / fac;
/*<             hc(k) = (hc(k)-hc(k-1))/fac >*/
		hc[k - 1] = (hc[k - 1] - hc[k - 2]) / fac;
/*<             k = k-1 >*/
		--k;
/*<             li = li-1 >*/
		--li;
/*<   50    continue >*/
/* L50: */
	    }
	}
/*<         s2 = (hs(5)-hs(4))*term >*/
	s2 = (hs[4] - hs[3]) * term;
/*<         c2 = (hc(5)-hc(4))*term >*/
	c2 = (hc[4] - hc[3]) * term;
/*<         go to 130 >*/
	goto L130;
/*  if !beta! <= 1 the integrals are calculated by evaluating a series */
/*  expansion. */
/*<   60    f3 = 0. >*/
L60:
	f3 = 0.;
/*<         do 70 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           ipj = i+j >*/
	    ipj = i__ + j;
/*<           hs(i) = par*(t(ipj)-t(j)) >*/
	    hs[i__ - 1] = *par * (t[ipj] - t[j]);
/*<           hc(i) = hs(i) >*/
	    hc[i__ - 1] = hs[i__ - 1];
/*<           f3 = f3+hs(i) >*/
	    f3 += hs[i__ - 1];
/*<   70    continue >*/
/* L70: */
	}
/*<         f3 = f3*con1 >*/
	f3 *= con1;
/*<         c1 = quart >*/
	c1 = quart;
/*<         s1 = f3 >*/
	s1 = f3;
/*<         if(abs(f3).le.eps) go to 120 >*/
	if (abs(f3) <= eps) {
	    goto L120;
	}
/*<         sign = one >*/
	sign = one;
/*<         fac = con2 >*/
	fac = con2;
/*<         k = 5 >*/
	k = 5;
/*<         is = 0 >*/
	is = 0;
/*<         do 110 ic=1,20 >*/
	for (ic = 1; ic <= 20; ++ic) {
/*<           k = k+1 >*/
	    ++k;
/*<           ak = k >*/
	    ak = (doublereal) k;
/*<           fac = fac*ak >*/
	    fac *= ak;
/*<           f1 = 0. >*/
	    f1 = 0.;
/*<           f3 = 0. >*/
	    f3 = 0.;
/*<           do 80 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<             f1 = f1+hc(i) >*/
		f1 += hc[i__ - 1];
/*<             f2 = f1*hs(i) >*/
		f2 = f1 * hs[i__ - 1];
/*<             hc(i) = f2 >*/
		hc[i__ - 1] = f2;
/*<             f3 = f3+f2 >*/
		f3 += f2;
/*<   80      continue >*/
/* L80: */
	    }
/*<           f3 = f3*six/fac >*/
	    f3 = f3 * six / fac;
/*<           if(is.eq.0) go to 90 >*/
	    if (is == 0) {
		goto L90;
	    }
/*<           is = 0 >*/
	    is = 0;
/*<           s1 = s1+f3*sign >*/
	    s1 += f3 * sign;
/*<           go to 100 >*/
	    goto L100;
/*<   90      sign = -sign >*/
L90:
	    sign = -sign;
/*<           is = 1 >*/
	    is = 1;
/*<           c1 = c1+f3*sign >*/
	    c1 += f3 * sign;
/*<  100      if(abs(f3).le.eps) go to 120 >*/
L100:
	    if (abs(f3) <= eps) {
		goto L120;
	    }
/*<  110    continue >*/
/* L110: */
	}
/*<  120    s2 = delta*(co(1)*s1+si(1)*c1) >*/
L120:
	s2 = delta * (co[0] * s1 + si[0] * c1);
/*<         c2 = delta*(co(1)*c1-si(1)*s1) >*/
	c2 = delta * (co[0] * c1 - si[0] * s1);
/*<  130    ress(j) = s2 >*/
L130:
	ress[j] = s2;
/*<         resc(j) = c2 >*/
	resc[j] = c2;
/*<         do 140 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           co(i) = co(i+1) >*/
	    co[i__ - 1] = co[i__];
/*<           si(i) = si(i+1) >*/
	    si[i__ - 1] = si[i__];
/*<  140    continue >*/
/* L140: */
	}
/*<  150  continue >*/
/* L150: */
    }
/*  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting */
/*  up a divided difference table. */
/*<  160  do 190 j=1,3 >*/
L160:
    for (j = 1; j <= 3; ++j) {
/*<         nmj = nm3-j >*/
	nmj = nm3 - j;
/*<         i = 5-j >*/
	i__ = 5 - j;
/*<    >*/
	fpcsin_(&t[nm3], &t[nmj], par, &si[3], &co[3], &si[i__ - 2], &co[i__ 
		- 2], &rs[j - 1], &rc[j - 1]);
/*<         hs(i) = 0. >*/
	hs[i__ - 1] = 0.;
/*<         hc(i) = 0. >*/
	hc[i__ - 1] = 0.;
/*<         do 170 jj=1,j >*/
	i__1 = j;
	for (jj = 1; jj <= i__1; ++jj) {
/*<           ipj = i+jj >*/
	    ipj = i__ + jj;
/*<           hc(ipj) = rc(jj) >*/
	    hc[ipj - 1] = rc[jj - 1];
/*<           hs(ipj) = rs(jj) >*/
	    hs[ipj - 1] = rs[jj - 1];
/*<  170    continue >*/
/* L170: */
	}
/*<         do 180 jj=1,3 >*/
	for (jj = 1; jj <= 3; ++jj) {
/*<           if(i.lt.jj) i = jj >*/
	    if (i__ < jj) {
		i__ = jj;
	    }
/*<           k = 5 >*/
	    k = 5;
/*<           li = nmj >*/
	    li = nmj;
/*<           do 180 ll=i,4 >*/
	    for (ll = i__; ll <= 4; ++ll) {
/*<             lj = li+jj >*/
		lj = li + jj;
/*<             fac = t(lj)-t(li) >*/
		fac = t[lj] - t[li];
/*<             hs(k) = (hs(k-1)-hs(k))/fac >*/
		hs[k - 1] = (hs[k - 2] - hs[k - 1]) / fac;
/*<             hc(k) = (hc(k-1)-hc(k))/fac >*/
		hc[k - 1] = (hc[k - 2] - hc[k - 1]) / fac;
/*<             k = k-1 >*/
		--k;
/*<             li = li+1 >*/
		++li;
/*<  180    continue >*/
/* L180: */
	    }
	}
/*<         ress(nmj) = hs(4)-hs(5) >*/
	ress[nmj] = hs[3] - hs[4];
/*<         resc(nmj) = hc(4)-hc(5) >*/
	resc[nmj] = hc[3] - hc[4];
/*<  190  continue >*/
/* L190: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpbfou_ */

