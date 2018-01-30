/* fpperi.f -- translated by f2c (version 20061008).
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
/* Subroutine */ int fpperi_(integer *iopt, doublereal *x, doublereal *y, 
	doublereal *w, integer *m, integer *k, doublereal *s, integer *nest, 
	doublereal *tol, integer *maxit, integer *k1, integer *k2, integer *n,
	 doublereal *t, doublereal *c__, doublereal *fp, doublereal *fpint, 
	doublereal *z__, doublereal *a1, doublereal *a2, doublereal *b, 
	doublereal *g1, doublereal *g2, doublereal *q, integer *nrdata, 
	integer *ier)
{
    /* System generated locals */
    integer a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, b_offset, g1_dim1,
	     g1_offset, g2_dim1, g2_offset, q_dim1, q_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l;
    static doublereal p, c1, d1, f1, f2, f3;
    static integer i1, i2, i3;
    static doublereal p1, p2, p3;
    static integer j1, j2, k3, l0, l1, l5, m1, n7, n8;
    static doublereal h1[7], h2[6];
    static integer n10, n11, ij, ik, jk, kk, mm, it;
    static doublereal wi, xi, yi, rn, fp0;
    static integer kk1, nk1, nk2;
    static doublereal acc, one, cos__, per, sin__;
    static integer new__;
    static doublereal piv;
    static integer ich1, ich3;
    static doublereal con1, con4, con9;
    static integer npl1;
    static doublereal half;
    static integer jper, nmin, iter, nmax;
    static doublereal fpms, term, pinv, fpold, fpart;
    static integer nrint;
    static doublereal store;
    static integer nplus;
    extern /* Subroutine */ int fpbacp_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), fpdisc_(doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), fprota_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpgivs_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpknot_(doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);

/*  .. */
/*  ..scalar arguments.. */
/*<       real s,tol,fp >*/
/*<       integer iopt,m,k,nest,maxit,k1,k2,n,ier >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer nrdata(nest) >*/
/*  ..local scalars.. */
/*<    >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real h(6),h1(7),h2(6) >*/
/*  ..function references.. */
/*<       real abs,fprati >*/
/*<       integer max0,min0 >*/
/*  ..subroutine references.. */
/*    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota */
/*  .. */
/*  set constants */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --w;
    --y;
    --x;
    --nrdata;
    a2_dim1 = *nest;
    a2_offset = 1 + a2_dim1;
    a2 -= a2_offset;
    --z__;
    --fpint;
    --c__;
    --t;
    q_dim1 = *m;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    g2_dim1 = *nest;
    g2_offset = 1 + g2_dim1;
    g2 -= g2_offset;
    a1_dim1 = *nest;
    a1_offset = 1 + a1_dim1;
    a1 -= a1_offset;
    g1_dim1 = *nest;
    g1_offset = 1 + g1_dim1;
    g1 -= g1_offset;
    b_dim1 = *nest;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    one = 1.;
/*<       con1 = 0.1e0 >*/
    con1 = .1;
/*<       con9 = 0.9e0 >*/
    con9 = .9;
/*<       con4 = 0.4e-01 >*/
    con4 = .04;
/*<       half = 0.5e0 >*/
    half = .5;
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 1: determination of the number of knots and their position     c */
/*  **************************************************************      c */
/*  given a set of knots we compute the least-squares periodic spline   c */
/*  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c */
/*  the initial choice of knots depends on the value of s and iopt.     c */
/*    if s=0 we have spline interpolation; in that case the number of   c */
/*    knots equals nmax = m+2*k.                                        c */
/*    if s > 0 and                                                      c */
/*      iopt=0 we first compute the least-squares polynomial of         c */
/*      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c */
/*      find that s(x) is a constant function.                          c */
/*      iopt=1 we start with the set of knots found at the last         c */
/*      call of the routine, except for the case that s > fp0; then     c */
/*      we compute directly the least-squares periodic polynomial.      c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*<       m1 = m-1 >*/
    m1 = *m - 1;
/*<       kk = k >*/
    kk = *k;
/*<       kk1 = k1 >*/
    kk1 = *k1;
/*<       k3 = 3*k+1 >*/
    k3 = *k * 3 + 1;
/*<       nmin = 2*k1 >*/
    nmin = *k1 << 1;
/*  determine the length of the period of s(x). */
/*<       per = x(m)-x(1) >*/
    per = x[*m] - x[1];
/*<       if(iopt.lt.0) go to 50 >*/
    if (*iopt < 0) {
	goto L50;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
/*<       acc = tol*s >*/
    acc = *tol * *s;
/*  determine nmax, the number of knots for periodic spline interpolation */
/*<       nmax = m+2*k >*/
    nmax = *m + (*k << 1);
/*<       if(s.gt.0. .or. nmax.eq.nmin) go to 30 >*/
    if (*s > 0. || nmax == nmin) {
	goto L30;
    }
/*  if s=0, s(x) is an interpolating spline. */
/*<       n = nmax >*/
    *n = nmax;
/*  test whether the required storage space exceeds the available one. */
/*<       if(n.gt.nest) go to 620 >*/
    if (*n > *nest) {
	goto L620;
    }
/*  find the position of the interior knots in case of interpolation. */
/*<    5  if((k/2)*2 .eq. k) go to 20 >*/
L5:
    if (*k / 2 << 1 == *k) {
	goto L20;
    }
/*<       do 10 i=2,m1 >*/
    i__1 = m1;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         j = i+k >*/
	j = i__ + *k;
/*<         t(j) = x(i) >*/
	t[j] = x[i__];
/*<   10  continue >*/
/* L10: */
    }
/*<       if(s.gt.0.) go to 50 >*/
    if (*s > 0.) {
	goto L50;
    }
/*<       kk = k-1 >*/
    kk = *k - 1;
/*<       kk1 = k >*/
    kk1 = *k;
/*<       if(kk.gt.0) go to 50 >*/
    if (kk > 0) {
	goto L50;
    }
/*<       t(1) = t(m)-per >*/
    t[1] = t[*m] - per;
/*<       t(2) = x(1) >*/
    t[2] = x[1];
/*<       t(m+1) = x(m) >*/
    t[*m + 1] = x[*m];
/*<       t(m+2) = t(3)+per >*/
    t[*m + 2] = t[3] + per;
/*<       do 15 i=1,m1 >*/
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         c(i) = y(i) >*/
	c__[i__] = y[i__];
/*<   15  continue >*/
/* L15: */
    }
/*<       c(m) = c(1) >*/
    c__[*m] = c__[1];
/*<       fp = 0. >*/
    *fp = 0.;
/*<       fpint(n) = fp0 >*/
    fpint[*n] = fp0;
/*<       fpint(n-1) = 0. >*/
    fpint[*n - 1] = 0.;
/*<       nrdata(n) = 0 >*/
    nrdata[*n] = 0;
/*<       go to 630 >*/
    goto L630;
/*<   20  do 25 i=2,m1 >*/
L20:
    i__1 = m1;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<         j = i+k >*/
	j = i__ + *k;
/*<         t(j) = (x(i)+x(i-1))*half >*/
	t[j] = (x[i__] + x[i__ - 1]) * half;
/*<   25  continue >*/
/* L25: */
    }
/*<       go to 50 >*/
    goto L50;
/*  if s > 0 our initial choice depends on the value of iopt. */
/*  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares */
/*  periodic polynomial. (i.e. a constant function). */
/*  if iopt=1 and fp0>s we start computing the least-squares periodic */
/*  spline according the set of knots found at the last call of the */
/*  routine. */
/*<   30  if(iopt.eq.0) go to 35 >*/
L30:
    if (*iopt == 0) {
	goto L35;
    }
/*<       if(n.eq.nmin) go to 35 >*/
    if (*n == nmin) {
	goto L35;
    }
/*<       fp0 = fpint(n) >*/
    fp0 = fpint[*n];
/*<       fpold = fpint(n-1) >*/
    fpold = fpint[*n - 1];
/*<       nplus = nrdata(n) >*/
    nplus = nrdata[*n];
/*<       if(fp0.gt.s) go to 50 >*/
    if (fp0 > *s) {
	goto L50;
    }
/*  the case that s(x) is a constant function is treated separetely. */
/*  find the least-squares constant c1 and compute fp0 at the same time. */
/*<   35  fp0 = 0. >*/
L35:
    fp0 = 0.;
/*<       d1 = 0. >*/
    d1 = 0.;
/*<       c1 = 0. >*/
    c1 = 0.;
/*<       do 40 it=1,m1 >*/
    i__1 = m1;
    for (it = 1; it <= i__1; ++it) {
/*<         wi = w(it) >*/
	wi = w[it];
/*<         yi = y(it)*wi >*/
	yi = y[it] * wi;
/*<         call fpgivs(wi,d1,cos,sin) >*/
	fpgivs_(&wi, &d1, &cos__, &sin__);
/*<         call fprota(cos,sin,yi,c1) >*/
	fprota_(&cos__, &sin__, &yi, &c1);
/*<         fp0 = fp0+yi**2 >*/
/* Computing 2nd power */
	d__1 = yi;
	fp0 += d__1 * d__1;
/*<   40  continue >*/
/* L40: */
    }
/*<       c1 = c1/d1 >*/
    c1 /= d1;
/*  test whether that constant function is a solution of our problem. */
/*<       fpms = fp0-s >*/
    fpms = fp0 - *s;
/*<       if(fpms.lt.acc .or. nmax.eq.nmin) go to 640 >*/
    if (fpms < acc || nmax == nmin) {
	goto L640;
    }
/*<       fpold = fp0 >*/
    fpold = fp0;
/*  test whether the required storage space exceeds the available one. */
/*<       if(nmin.ge.nest) go to 620 >*/
    if (nmin >= *nest) {
	goto L620;
    }
/*  start computing the least-squares periodic spline with one */
/*  interior knot. */
/*<       nplus = 1 >*/
    nplus = 1;
/*<       n = nmin+1 >*/
    *n = nmin + 1;
/*<       mm = (m+1)/2 >*/
    mm = (*m + 1) / 2;
/*<       t(k2) = x(mm) >*/
    t[*k2] = x[mm];
/*<       nrdata(1) = mm-2 >*/
    nrdata[1] = mm - 2;
/*<       nrdata(2) = m1-mm >*/
    nrdata[2] = m1 - mm;
/*  main loop for the different sets of knots. m is a save upper */
/*  bound for the number of trials. */
/*<   50  do 340 iter=1,m >*/
L50:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find nrint, the number of knot intervals. */
/*<         nrint = n-nmin+1 >*/
	nrint = *n - nmin + 1;
/*  find the position of the additional knots which are needed for */
/*  the b-spline representation of s(x). if we take */
/*      t(k+1) = x(1), t(n-k) = x(m) */
/*      t(k+1-j) = t(n-k-j) - per, j=1,2,...k */
/*      t(n-k+j) = t(k+1+j) + per, j=1,2,...k */
/*  then s(x) is a periodic spline with period per if the b-spline */
/*  coefficients satisfy the following conditions */
/*      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1. */
/*<         t(k1) = x(1) >*/
	t[*k1] = x[1];
/*<         nk1 = n-k1 >*/
	nk1 = *n - *k1;
/*<         nk2 = nk1+1 >*/
	nk2 = nk1 + 1;
/*<         t(nk2) = x(m) >*/
	t[nk2] = x[*m];
/*<         do 60 j=1,k >*/
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
/*<           i1 = nk2+j >*/
	    i1 = nk2 + j;
/*<           i2 = nk2-j >*/
	    i2 = nk2 - j;
/*<           j1 = k1+j >*/
	    j1 = *k1 + j;
/*<           j2 = k1-j >*/
	    j2 = *k1 - j;
/*<           t(i1) = t(j1)+per >*/
	    t[i1] = t[j1] + per;
/*<           t(j2) = t(i2)-per >*/
	    t[j2] = t[i2] - per;
/*<   60    continue >*/
/* L60: */
	}
/*  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares */
/*  periodic spline sinf(x). the observation matrix a is built up row */
/*  by row while taking into account condition (**) and is reduced to */
/*  triangular form by givens transformations . */
/*  at the same time fp=f(p=inf) is computed. */
/*  the n7 x n7 triangularised upper matrix a has the form */
/*            ! a1 '    ! */
/*        a = !    ' a2 ! */
/*            ! 0  '    ! */
/*  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular */
/*  matrix of bandwith k+1 ( n10 = n7-k). */
/*  initialization. */
/*<         do 70 i=1,nk1 >*/
	i__2 = nk1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           z(i) = 0. >*/
	    z__[i__] = 0.;
/*<           do 70 j=1,kk1 >*/
	    i__3 = kk1;
	    for (j = 1; j <= i__3; ++j) {
/*<             a1(i,j) = 0. >*/
		a1[i__ + j * a1_dim1] = 0.;
/*<   70    continue >*/
/* L70: */
	    }
	}
/*<         n7 = nk1-k >*/
	n7 = nk1 - *k;
/*<         n10 = n7-kk >*/
	n10 = n7 - kk;
/*<         jper = 0 >*/
	jper = 0;
/*<         fp = 0. >*/
	*fp = 0.;
/*<         l = k1 >*/
	l = *k1;
/*<         do 290 it=1,m1 >*/
	i__3 = m1;
	for (it = 1; it <= i__3; ++it) {
/*  fetch the current data point x(it),y(it) */
/*<           xi = x(it) >*/
	    xi = x[it];
/*<           wi = w(it) >*/
	    wi = w[it];
/*<           yi = y(it)*wi >*/
	    yi = y[it] * wi;
/*  search for knot interval t(l) <= xi < t(l+1). */
/*<   80      if(xi.lt.t(l+1)) go to 85 >*/
L80:
	    if (xi < t[l + 1]) {
		goto L85;
	    }
/*<           l = l+1 >*/
	    ++l;
/*<           go to 80 >*/
	    goto L80;
/*  evaluate the (k+1) non-zero b-splines at xi and store them in q. */
/*<   85      call fpbspl(t,n,k,xi,l,h) >*/
L85:
	    fpbspl_(&t[1], n, k, &xi, &l, h__);
/*<           do 90 i=1,k1 >*/
	    i__2 = *k1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             q(it,i) = h(i) >*/
		q[it + i__ * q_dim1] = h__[i__ - 1];
/*<             h(i) = h(i)*wi >*/
		h__[i__ - 1] *= wi;
/*<   90      continue >*/
/* L90: */
	    }
/*<           l5 = l-k1 >*/
	    l5 = l - *k1;
/*  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi */
/*<           if(l5.lt.n10) go to 285 >*/
	    if (l5 < n10) {
		goto L285;
	    }
/*<           if(jper.ne.0) go to 160 >*/
	    if (jper != 0) {
		goto L160;
	    }
/*  initialize the matrix a2. */
/*<           do 95 i=1,n7 >*/
	    i__2 = n7;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<           do 95 j=1,kk >*/
		i__4 = kk;
		for (j = 1; j <= i__4; ++j) {
/*<               a2(i,j) = 0. >*/
		    a2[i__ + j * a2_dim1] = 0.;
/*<   95      continue >*/
/* L95: */
		}
	    }
/*<           jk = n10+1 >*/
	    jk = n10 + 1;
/*<           do 110 i=1,kk >*/
	    i__4 = kk;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<             ik = jk >*/
		ik = jk;
/*<             do 100 j=1,kk1 >*/
		i__2 = kk1;
		for (j = 1; j <= i__2; ++j) {
/*<               if(ik.le.0) go to 105 >*/
		    if (ik <= 0) {
			goto L105;
		    }
/*<               a2(ik,i) = a1(ik,j) >*/
		    a2[ik + i__ * a2_dim1] = a1[ik + j * a1_dim1];
/*<               ik = ik-1 >*/
		    --ik;
/*<  100        continue >*/
/* L100: */
		}
/*<  105        jk = jk+1 >*/
L105:
		++jk;
/*<  110      continue >*/
/* L110: */
	    }
/*<           jper = 1 >*/
	    jper = 1;
/*  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi */
/*  we take account of condition (**) for setting up the new row */
/*  of the observation matrix a. this row is stored in the arrays h1 */
/*  (the part with respect to a1) and h2 (the part with */
/*  respect to a2). */
/*<  160      do 170 i=1,kk >*/
L160:
	    i__4 = kk;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<             h1(i) = 0. >*/
		h1[i__ - 1] = 0.;
/*<             h2(i) = 0. >*/
		h2[i__ - 1] = 0.;
/*<  170      continue >*/
/* L170: */
	    }
/*<           h1(kk1) = 0. >*/
	    h1[kk1 - 1] = 0.;
/*<           j = l5-n10 >*/
	    j = l5 - n10;
/*<           do 210 i=1,kk1 >*/
	    i__4 = kk1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<             j = j+1 >*/
		++j;
/*<             l0 = j >*/
		l0 = j;
/*<  180        l1 = l0-kk >*/
L180:
		l1 = l0 - kk;
/*<             if(l1.le.0) go to 200 >*/
		if (l1 <= 0) {
		    goto L200;
		}
/*<             if(l1.le.n10) go to 190 >*/
		if (l1 <= n10) {
		    goto L190;
		}
/*<             l0 = l1-n10 >*/
		l0 = l1 - n10;
/*<             go to 180 >*/
		goto L180;
/*<  190        h1(l1) = h(i) >*/
L190:
		h1[l1 - 1] = h__[i__ - 1];
/*<             go to 210 >*/
		goto L210;
/*<  200        h2(l0) = h2(l0)+h(i) >*/
L200:
		h2[l0 - 1] += h__[i__ - 1];
/*<  210      continue >*/
L210:
		;
	    }
/*  rotate the new row of the observation matrix into triangle */
/*  by givens transformations. */
/*<           if(n10.le.0) go to 250 >*/
	    if (n10 <= 0) {
		goto L250;
	    }
/*  rotation with the rows 1,2,...n10 of matrix a. */
/*<           do 240 j=1,n10 >*/
	    i__4 = n10;
	    for (j = 1; j <= i__4; ++j) {
/*<             piv = h1(1) >*/
		piv = h1[0];
/*<             if(piv.ne.0.) go to 214 >*/
		if (piv != 0.) {
		    goto L214;
		}
/*<             do 212 i=1,kk >*/
		i__2 = kk;
		for (i__ = 1; i__ <= i__2; ++i__) {
/*<               h1(i) = h1(i+1) >*/
		    h1[i__ - 1] = h1[i__];
/*<  212        continue >*/
/* L212: */
		}
/*<             h1(kk1) = 0. >*/
		h1[kk1 - 1] = 0.;
/*<             go to 240 >*/
		goto L240;
/*  calculate the parameters of the givens transformation. */
/*<  214        call fpgivs(piv,a1(j,1),cos,sin) >*/
L214:
		fpgivs_(&piv, &a1[j + a1_dim1], &cos__, &sin__);
/*  transformation to the right hand side. */
/*<             call fprota(cos,sin,yi,z(j)) >*/
		fprota_(&cos__, &sin__, &yi, &z__[j]);
/*  transformations to the left hand side with respect to a2. */
/*<             do 220 i=1,kk >*/
		i__2 = kk;
		for (i__ = 1; i__ <= i__2; ++i__) {
/*<               call fprota(cos,sin,h2(i),a2(j,i)) >*/
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &a2[j + i__ * 
			    a2_dim1]);
/*<  220        continue >*/
/* L220: */
		}
/*<             if(j.eq.n10) go to 250 >*/
		if (j == n10) {
		    goto L250;
		}
/*<             i2 = min0(n10-j,kk) >*/
/* Computing MIN */
		i__2 = n10 - j;
		i2 = min(i__2,kk);
/*  transformations to the left hand side with respect to a1. */
/*<             do 230 i=1,i2 >*/
		i__2 = i2;
		for (i__ = 1; i__ <= i__2; ++i__) {
/*<               i1 = i+1 >*/
		    i1 = i__ + 1;
/*<               call fprota(cos,sin,h1(i1),a1(j,i1)) >*/
		    fprota_(&cos__, &sin__, &h1[i1 - 1], &a1[j + i1 * a1_dim1]
			    );
/*<               h1(i) = h1(i1) >*/
		    h1[i__ - 1] = h1[i1 - 1];
/*<  230        continue >*/
/* L230: */
		}
/*<             h1(i1) = 0. >*/
		h1[i1 - 1] = 0.;
/*<  240      continue >*/
L240:
		;
	    }
/*  rotation with the rows n10+1,...n7 of matrix a. */
/*<  250      do 270 j=1,kk >*/
L250:
	    i__4 = kk;
	    for (j = 1; j <= i__4; ++j) {
/*<             ij = n10+j >*/
		ij = n10 + j;
/*<             if(ij.le.0) go to 270 >*/
		if (ij <= 0) {
		    goto L270;
		}
/*<             piv = h2(j) >*/
		piv = h2[j - 1];
/*<             if(piv.eq.0.) go to 270 >*/
		if (piv == 0.) {
		    goto L270;
		}
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a2(ij,j),cos,sin) >*/
		fpgivs_(&piv, &a2[ij + j * a2_dim1], &cos__, &sin__);
/*  transformations to right hand side. */
/*<             call fprota(cos,sin,yi,z(ij)) >*/
		fprota_(&cos__, &sin__, &yi, &z__[ij]);
/*<             if(j.eq.kk) go to 280 >*/
		if (j == kk) {
		    goto L280;
		}
/*<             j1 = j+1 >*/
		j1 = j + 1;
/*  transformations to left hand side. */
/*<             do 260 i=j1,kk >*/
		i__2 = kk;
		for (i__ = j1; i__ <= i__2; ++i__) {
/*<               call fprota(cos,sin,h2(i),a2(ij,i)) >*/
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &a2[ij + i__ * 
			    a2_dim1]);
/*<  260        continue >*/
/* L260: */
		}
/*<  270      continue >*/
L270:
		;
	    }
/*  add contribution of this row to the sum of squares of residual */
/*  right hand sides. */
/*<  280      fp = fp+yi**2 >*/
L280:
/* Computing 2nd power */
	    d__1 = yi;
	    *fp += d__1 * d__1;
/*<           go to 290 >*/
	    goto L290;
/*  rotation of the new row of the observation matrix into */
/*  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero */
/*  at xi. */
/*<  285      j = l5 >*/
L285:
	    j = l5;
/*<           do 140 i=1,kk1 >*/
	    i__4 = kk1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<             j = j+1 >*/
		++j;
/*<             piv = h(i) >*/
		piv = h__[i__ - 1];
/*<             if(piv.eq.0.) go to 140 >*/
		if (piv == 0.) {
		    goto L140;
		}
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a1(j,1),cos,sin) >*/
		fpgivs_(&piv, &a1[j + a1_dim1], &cos__, &sin__);
/*  transformations to right hand side. */
/*<             call fprota(cos,sin,yi,z(j)) >*/
		fprota_(&cos__, &sin__, &yi, &z__[j]);
/*<             if(i.eq.kk1) go to 150 >*/
		if (i__ == kk1) {
		    goto L150;
		}
/*<             i2 = 1 >*/
		i2 = 1;
/*<             i3 = i+1 >*/
		i3 = i__ + 1;
/*  transformations to left hand side. */
/*<             do 130 i1=i3,kk1 >*/
		i__2 = kk1;
		for (i1 = i3; i1 <= i__2; ++i1) {
/*<               i2 = i2+1 >*/
		    ++i2;
/*<               call fprota(cos,sin,h(i1),a1(j,i2)) >*/
		    fprota_(&cos__, &sin__, &h__[i1 - 1], &a1[j + i2 * 
			    a1_dim1]);
/*<  130        continue >*/
/* L130: */
		}
/*<  140      continue >*/
L140:
		;
	    }
/*  add contribution of this row to the sum of squares of residual */
/*  right hand sides. */
/*<  150      fp = fp+yi**2 >*/
L150:
/* Computing 2nd power */
	    d__1 = yi;
	    *fp += d__1 * d__1;
/*<  290    continue >*/
L290:
	    ;
	}
/*<         fpint(n) = fp0 >*/
	fpint[*n] = fp0;
/*<         fpint(n-1) = fpold >*/
	fpint[*n - 1] = fpold;
/*<         nrdata(n) = nplus >*/
	nrdata[*n] = nplus;
/*  backward substitution to obtain the b-spline coefficients c(j),j=1,.n */
/*<         call fpbacp(a1,a2,z,n7,kk,c,kk1,nest) >*/
	fpbacp_(&a1[a1_offset], &a2[a2_offset], &z__[1], &n7, &kk, &c__[1], &
		kk1, nest);
/*  calculate from condition (**) the coefficients c(j+n7),j=1,2,...k. */
/*<         do 295 i=1,k >*/
	i__3 = *k;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           j = i+n7 >*/
	    j = i__ + n7;
/*<           c(j) = c(i) >*/
	    c__[j] = c__[i__];
/*<  295    continue >*/
/* L295: */
	}
/*<         if(iopt.lt.0) go to 660 >*/
	if (*iopt < 0) {
	    goto L660;
	}
/*  test whether the approximation sinf(x) is an acceptable solution. */
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).lt.acc) go to 660 >*/
	if (abs(fpms) < acc) {
	    goto L660;
	}
/*  if f(p=inf) < s accept the choice of knots. */
/*<         if(fpms.lt.0.) go to 350 >*/
	if (fpms < 0.) {
	    goto L350;
	}
/*  if n=nmax, sinf(x) is an interpolating spline. */
/*<         if(n.eq.nmax) go to 630 >*/
	if (*n == nmax) {
	    goto L630;
	}
/*  increase the number of knots. */
/*  if n=nest we cannot increase the number of knots because of the */
/*  storage capacity limitation. */
/*<         if(n.eq.nest) go to 620 >*/
	if (*n == *nest) {
	    goto L620;
	}
/*  determine the number of knots nplus we are going to add. */
/*<         npl1 = nplus*2 >*/
	npl1 = nplus << 1;
/*<         rn = nplus >*/
	rn = (doublereal) nplus;
/*<         if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp) >*/
	if (fpold - *fp > acc) {
	    npl1 = (integer) (rn * fpms / (fpold - *fp));
	}
/*<         nplus = min0(nplus*2,max0(npl1,nplus/2,1)) >*/
/* Computing MIN */
/* Computing MAX */
	i__2 = npl1, i__5 = nplus / 2, i__2 = max(i__2,i__5);
	i__3 = nplus << 1, i__4 = max(i__2,1);
	nplus = min(i__3,i__4);
/*<         fpold = fp >*/
	fpold = *fp;
/*  compute the sum(wi*(yi-s(xi))**2) for each knot interval */
/*  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
/*<         fpart = 0. >*/
	fpart = 0.;
/*<         i = 1 >*/
	i__ = 1;
/*<         l = k1 >*/
	l = *k1;
/*<         do 320 it=1,m1 >*/
	i__3 = m1;
	for (it = 1; it <= i__3; ++it) {
/*<           if(x(it).lt.t(l)) go to 300 >*/
	    if (x[it] < t[l]) {
		goto L300;
	    }
/*<           new = 1 >*/
	    new__ = 1;
/*<           l = l+1 >*/
	    ++l;
/*<  300      term = 0. >*/
L300:
	    term = 0.;
/*<           l0 = l-k2 >*/
	    l0 = l - *k2;
/*<           do 310 j=1,k1 >*/
	    i__4 = *k1;
	    for (j = 1; j <= i__4; ++j) {
/*<             l0 = l0+1 >*/
		++l0;
/*<             term = term+c(l0)*q(it,j) >*/
		term += c__[l0] * q[it + j * q_dim1];
/*<  310      continue >*/
/* L310: */
	    }
/*<           term = (w(it)*(term-y(it)))**2 >*/
/* Computing 2nd power */
	    d__1 = w[it] * (term - y[it]);
	    term = d__1 * d__1;
/*<           fpart = fpart+term >*/
	    fpart += term;
/*<           if(new.eq.0) go to 320 >*/
	    if (new__ == 0) {
		goto L320;
	    }
/*<           if(l.gt.k2) go to 315 >*/
	    if (l > *k2) {
		goto L315;
	    }
/*<           fpint(nrint) = term >*/
	    fpint[nrint] = term;
/*<           new = 0 >*/
	    new__ = 0;
/*<           go to 320 >*/
	    goto L320;
/*<  315      store = term*half >*/
L315:
	    store = term * half;
/*<           fpint(i) = fpart-store >*/
	    fpint[i__] = fpart - store;
/*<           i = i+1 >*/
	    ++i__;
/*<           fpart = store >*/
	    fpart = store;
/*<           new = 0 >*/
	    new__ = 0;
/*<  320    continue >*/
L320:
	    ;
	}
/*<         fpint(nrint) = fpint(nrint)+fpart >*/
	fpint[nrint] += fpart;
/*<         do 330 l=1,nplus >*/
	i__3 = nplus;
	for (l = 1; l <= i__3; ++l) {
/*  add a new knot */
/*<           call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1) >*/
	    fpknot_(&x[1], m, &t[1], n, &fpint[1], &nrdata[1], &nrint, nest, &
		    c__1);
/*  if n=nmax we locate the knots as for interpolation. */
/*<           if(n.eq.nmax) go to 5 >*/
	    if (*n == nmax) {
		goto L5;
	    }
/*  test whether we cannot further increase the number of knots. */
/*<           if(n.eq.nest) go to 340 >*/
	    if (*n == *nest) {
		goto L340;
	    }
/*<  330    continue >*/
/* L330: */
	}
/*  restart the computations with the new set of knots. */
/*<  340  continue >*/
L340:
	;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 2: determination of the smoothing periodic spline sp(x).       c */
/*  *************************************************************       c */
/*  we have determined the number of knots and their position.          c */
/*  we now compute the b-spline coefficients of the smoothing spline    c */
/*  sp(x). the observation matrix a is extended by the rows of matrix   c */
/*  b expressing that the kth derivative discontinuities of sp(x) at    c */
/*  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c */
/*  ponding weights of these additional rows are set to 1/sqrt(p).      c */
/*  iteratively we then have to determine the value of p such that      c */
/*  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c */
/*  the least-squares constant function corresponds to p=0, and that    c */
/*  the least-squares periodic spline corresponds to p=infinity. the    c */
/*  iteration process which is proposed here, makes use of rational     c */
/*  interpolation. since f(p) is a convex and strictly decreasing       c */
/*  function of p, it can be approximated by a rational function        c */
/*  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c */
/*  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c */
/*  to calculate the new value of p such that r(p)=s. convergence is    c */
/*  guaranteed by taking f1>0 and f3<0.                                 c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  evaluate the discontinuity jump of the kth derivative of the */
/*  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b. */
/*<  350  call fpdisc(t,n,k2,b,nest) >*/
L350:
    fpdisc_(&t[1], n, k2, &b[b_offset], nest);
/*  initial value for p. */
/*<       p1 = 0. >*/
    p1 = 0.;
/*<       f1 = fp0-s >*/
    f1 = fp0 - *s;
/*<       p3 = -one >*/
    p3 = -one;
/*<       f3 = fpms >*/
    f3 = fpms;
/*<       n11 = n10-1 >*/
    n11 = n10 - 1;
/*<       n8 = n7-1 >*/
    n8 = n7 - 1;
/*<       p = 0. >*/
    p = 0.;
/*<       l = n7 >*/
    l = n7;
/*<       do 352 i=1,k >*/
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          j = k+1-i >*/
	j = *k + 1 - i__;
/*<          p = p+a2(l,j) >*/
	p += a2[l + j * a2_dim1];
/*<          l = l-1 >*/
	--l;
/*<          if(l.eq.0) go to 356 >*/
	if (l == 0) {
	    goto L356;
	}
/*<  352  continue >*/
/* L352: */
    }
/*<       do 354 i=1,n10 >*/
    i__1 = n10;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          p = p+a1(i,1) >*/
	p += a1[i__ + a1_dim1];
/*<  354  continue >*/
/* L354: */
    }
/*<  356  rn = n7 >*/
L356:
    rn = (doublereal) n7;
/*<       p = rn/p >*/
    p = rn / p;
/*<       ich1 = 0 >*/
    ich1 = 0;
/*<       ich3 = 0 >*/
    ich3 = 0;
/*  iteration process to find the root of f(p) = s. */
/*<       do 595 iter=1,maxit >*/
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*  form the matrix g  as the matrix a extended by the rows of matrix b. */
/*  the rows of matrix b with weight 1/p are rotated into */
/*  the triangularised observation matrix a. */
/*  after triangularisation our n7 x n7 matrix g takes the form */
/*            ! g1 '    ! */
/*        g = !    ' g2 ! */
/*            ! 0  '    ! */
/*  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular */
/*  matrix of bandwidth k+2. ( n11 = n7-k-1) */
/*<         pinv = one/p >*/
	pinv = one / p;
/*  store matrix a into g */
/*<         do 360 i=1,n7 >*/
	i__3 = n7;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           c(i) = z(i) >*/
	    c__[i__] = z__[i__];
/*<           g1(i,k1) = a1(i,k1) >*/
	    g1[i__ + *k1 * g1_dim1] = a1[i__ + *k1 * a1_dim1];
/*<           g1(i,k2) = 0. >*/
	    g1[i__ + *k2 * g1_dim1] = 0.;
/*<           g2(i,1) = 0. >*/
	    g2[i__ + g2_dim1] = 0.;
/*<           do 360 j=1,k >*/
	    i__4 = *k;
	    for (j = 1; j <= i__4; ++j) {
/*<             g1(i,j) = a1(i,j) >*/
		g1[i__ + j * g1_dim1] = a1[i__ + j * a1_dim1];
/*<             g2(i,j+1) = a2(i,j) >*/
		g2[i__ + (j + 1) * g2_dim1] = a2[i__ + j * a2_dim1];
/*<  360    continue >*/
/* L360: */
	    }
	}
/*<         l = n10 >*/
	l = n10;
/*<         do 370 j=1,k1 >*/
	i__4 = *k1;
	for (j = 1; j <= i__4; ++j) {
/*<           if(l.le.0) go to 375 >*/
	    if (l <= 0) {
		goto L375;
	    }
/*<           g2(l,1) = a1(l,j) >*/
	    g2[l + g2_dim1] = a1[l + j * a1_dim1];
/*<           l = l-1 >*/
	    --l;
/*<  370    continue >*/
/* L370: */
	}
/*<  375    do 540 it=1,n8 >*/
L375:
	i__4 = n8;
	for (it = 1; it <= i__4; ++it) {
/*  fetch a new row of matrix b and store it in the arrays h1 (the part */
/*  with respect to g1) and h2 (the part with respect to g2). */
/*<           yi = 0. >*/
	    yi = 0.;
/*<           do 380 i=1,k1 >*/
	    i__3 = *k1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             h1(i) = 0. >*/
		h1[i__ - 1] = 0.;
/*<             h2(i) = 0. >*/
		h2[i__ - 1] = 0.;
/*<  380      continue >*/
/* L380: */
	    }
/*<           h1(k2) = 0. >*/
	    h1[*k2 - 1] = 0.;
/*<           if(it.gt.n11) go to 420 >*/
	    if (it > n11) {
		goto L420;
	    }
/*<           l = it >*/
	    l = it;
/*<           l0 = it >*/
	    l0 = it;
/*<           do 390 j=1,k2 >*/
	    i__3 = *k2;
	    for (j = 1; j <= i__3; ++j) {
/*<             if(l0.eq.n10) go to 400 >*/
		if (l0 == n10) {
		    goto L400;
		}
/*<             h1(j) = b(it,j)*pinv >*/
		h1[j - 1] = b[it + j * b_dim1] * pinv;
/*<             l0 = l0+1 >*/
		++l0;
/*<  390      continue >*/
/* L390: */
	    }
/*<           go to 470 >*/
	    goto L470;
/*<  400      l0 = 1 >*/
L400:
	    l0 = 1;
/*<           do 410 l1=j,k2 >*/
	    i__3 = *k2;
	    for (l1 = j; l1 <= i__3; ++l1) {
/*<             h2(l0) = b(it,l1)*pinv >*/
		h2[l0 - 1] = b[it + l1 * b_dim1] * pinv;
/*<             l0 = l0+1 >*/
		++l0;
/*<  410      continue >*/
/* L410: */
	    }
/*<           go to 470 >*/
	    goto L470;
/*<  420      l = 1 >*/
L420:
	    l = 1;
/*<           i = it-n10 >*/
	    i__ = it - n10;
/*<           do 460 j=1,k2 >*/
	    i__3 = *k2;
	    for (j = 1; j <= i__3; ++j) {
/*<             i = i+1 >*/
		++i__;
/*<             l0 = i >*/
		l0 = i__;
/*<  430        l1 = l0-k1 >*/
L430:
		l1 = l0 - *k1;
/*<             if(l1.le.0) go to 450 >*/
		if (l1 <= 0) {
		    goto L450;
		}
/*<             if(l1.le.n11) go to 440 >*/
		if (l1 <= n11) {
		    goto L440;
		}
/*<             l0 = l1-n11 >*/
		l0 = l1 - n11;
/*<             go to 430 >*/
		goto L430;
/*<  440        h1(l1) = b(it,j)*pinv >*/
L440:
		h1[l1 - 1] = b[it + j * b_dim1] * pinv;
/*<             go to 460 >*/
		goto L460;
/*<  450        h2(l0) = h2(l0)+b(it,j)*pinv >*/
L450:
		h2[l0 - 1] += b[it + j * b_dim1] * pinv;
/*<  460      continue >*/
L460:
		;
	    }
/*<           if(n11.le.0) go to 510 >*/
	    if (n11 <= 0) {
		goto L510;
	    }
/*  rotate this row into triangle by givens transformations without */
/*  square roots. */
/*  rotation with the rows l,l+1,...n11. */
/*<  470      do 500 j=l,n11 >*/
L470:
	    i__3 = n11;
	    for (j = l; j <= i__3; ++j) {
/*<             piv = h1(1) >*/
		piv = h1[0];
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,g1(j,1),cos,sin) >*/
		fpgivs_(&piv, &g1[j + g1_dim1], &cos__, &sin__);
/*  transformation to right hand side. */
/*<             call fprota(cos,sin,yi,c(j)) >*/
		fprota_(&cos__, &sin__, &yi, &c__[j]);
/*  transformation to the left hand side with respect to g2. */
/*<             do 480 i=1,k1 >*/
		i__2 = *k1;
		for (i__ = 1; i__ <= i__2; ++i__) {
/*<               call fprota(cos,sin,h2(i),g2(j,i)) >*/
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &g2[j + i__ * 
			    g2_dim1]);
/*<  480        continue >*/
/* L480: */
		}
/*<             if(j.eq.n11) go to 510 >*/
		if (j == n11) {
		    goto L510;
		}
/*<             i2 = min0(n11-j,k1) >*/
/* Computing MIN */
		i__2 = n11 - j;
		i2 = min(i__2,*k1);
/*  transformation to the left hand side with respect to g1. */
/*<             do 490 i=1,i2 >*/
		i__2 = i2;
		for (i__ = 1; i__ <= i__2; ++i__) {
/*<               i1 = i+1 >*/
		    i1 = i__ + 1;
/*<               call fprota(cos,sin,h1(i1),g1(j,i1)) >*/
		    fprota_(&cos__, &sin__, &h1[i1 - 1], &g1[j + i1 * g1_dim1]
			    );
/*<               h1(i) = h1(i1) >*/
		    h1[i__ - 1] = h1[i1 - 1];
/*<  490        continue >*/
/* L490: */
		}
/*<             h1(i1) = 0. >*/
		h1[i1 - 1] = 0.;
/*<  500      continue >*/
/* L500: */
	    }
/*  rotation with the rows n11+1,...n7 */
/*<  510      do 530 j=1,k1 >*/
L510:
	    i__3 = *k1;
	    for (j = 1; j <= i__3; ++j) {
/*<             ij = n11+j >*/
		ij = n11 + j;
/*<             if(ij.le.0) go to 530 >*/
		if (ij <= 0) {
		    goto L530;
		}
/*<             piv = h2(j) >*/
		piv = h2[j - 1];
/*  calculate the parameters of the givens transformation */
/*<             call fpgivs(piv,g2(ij,j),cos,sin) >*/
		fpgivs_(&piv, &g2[ij + j * g2_dim1], &cos__, &sin__);
/*  transformation to the right hand side. */
/*<             call fprota(cos,sin,yi,c(ij)) >*/
		fprota_(&cos__, &sin__, &yi, &c__[ij]);
/*<             if(j.eq.k1) go to 540 >*/
		if (j == *k1) {
		    goto L540;
		}
/*<             j1 = j+1 >*/
		j1 = j + 1;
/*  transformation to the left hand side. */
/*<             do 520 i=j1,k1 >*/
		i__2 = *k1;
		for (i__ = j1; i__ <= i__2; ++i__) {
/*<               call fprota(cos,sin,h2(i),g2(ij,i)) >*/
		    fprota_(&cos__, &sin__, &h2[i__ - 1], &g2[ij + i__ * 
			    g2_dim1]);
/*<  520        continue >*/
/* L520: */
		}
/*<  530      continue >*/
L530:
		;
	    }
/*<  540    continue >*/
L540:
	    ;
	}
/*  backward substitution to obtain the b-spline coefficients */
/*  c(j),j=1,2,...n7 of sp(x). */
/*<         call fpbacp(g1,g2,c,n7,k1,c,k2,nest) >*/
	fpbacp_(&g1[g1_offset], &g2[g2_offset], &c__[1], &n7, k1, &c__[1], k2,
		 nest);
/*  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,. */
/*<         do 545 i=1,k >*/
	i__4 = *k;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<           j = i+n7 >*/
	    j = i__ + n7;
/*<           c(j) = c(i) >*/
	    c__[j] = c__[i__];
/*<  545    continue >*/
/* L545: */
	}
/*  computation of f(p). */
/*<         fp = 0. >*/
	*fp = 0.;
/*<         l = k1 >*/
	l = *k1;
/*<         do 570 it=1,m1 >*/
	i__4 = m1;
	for (it = 1; it <= i__4; ++it) {
/*<           if(x(it).lt.t(l)) go to 550 >*/
	    if (x[it] < t[l]) {
		goto L550;
	    }
/*<           l = l+1 >*/
	    ++l;
/*<  550      l0 = l-k2 >*/
L550:
	    l0 = l - *k2;
/*<           term = 0. >*/
	    term = 0.;
/*<           do 560 j=1,k1 >*/
	    i__3 = *k1;
	    for (j = 1; j <= i__3; ++j) {
/*<             l0 = l0+1 >*/
		++l0;
/*<             term = term+c(l0)*q(it,j) >*/
		term += c__[l0] * q[it + j * q_dim1];
/*<  560      continue >*/
/* L560: */
	    }
/*<           fp = fp+(w(it)*(term-y(it)))**2 >*/
/* Computing 2nd power */
	    d__1 = w[it] * (term - y[it]);
	    *fp += d__1 * d__1;
/*<  570    continue >*/
/* L570: */
	}
/*  test whether the approximation sp(x) is an acceptable solution. */
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).lt.acc) go to 660 >*/
	if (abs(fpms) < acc) {
	    goto L660;
	}
/*  test whether the maximal number of iterations is reached. */
/*<         if(iter.eq.maxit) go to 600 >*/
	if (iter == *maxit) {
	    goto L600;
	}
/*  carry out one more step of the iteration process. */
/*<         p2 = p >*/
	p2 = p;
/*<         f2 = fpms >*/
	f2 = fpms;
/*<         if(ich3.ne.0) go to 580 >*/
	if (ich3 != 0) {
	    goto L580;
	}
/*<         if((f2-f3) .gt. acc) go to 575 >*/
	if (f2 - f3 > acc) {
	    goto L575;
	}
/*  our initial choice of p is too large. */
/*<         p3 = p2 >*/
	p3 = p2;
/*<         f3 = f2 >*/
	f3 = f2;
/*<         p = p*con4 >*/
	p *= con4;
/*<         if(p.le.p1) p = p1*con9 +p2*con1 >*/
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
/*<         go to 595 >*/
	goto L595;
/*<  575    if(f2.lt.0.) ich3 = 1 >*/
L575:
	if (f2 < 0.) {
	    ich3 = 1;
	}
/*<  580    if(ich1.ne.0) go to 590 >*/
L580:
	if (ich1 != 0) {
	    goto L590;
	}
/*<         if((f1-f2) .gt. acc) go to 585 >*/
	if (f1 - f2 > acc) {
	    goto L585;
	}
/*  our initial choice of p is too small */
/*<         p1 = p2 >*/
	p1 = p2;
/*<         f1 = f2 >*/
	f1 = f2;
/*<         p = p/con4 >*/
	p /= con4;
/*<         if(p3.lt.0.) go to 595 >*/
	if (p3 < 0.) {
	    goto L595;
	}
/*<         if(p.ge.p3) p = p2*con1 +p3*con9 >*/
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
/*<         go to 595 >*/
	goto L595;
/*<  585    if(f2.gt.0.) ich1 = 1 >*/
L585:
	if (f2 > 0.) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
/*<  590    if(f2.ge.f1 .or. f2.le.f3) go to 610 >*/
L590:
	if (f2 >= f1 || f2 <= f3) {
	    goto L610;
	}
/*  find the new value for p. */
/*<         p = fprati(p1,f1,p2,f2,p3,f3) >*/
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
/*<  595  continue >*/
L595:
	;
    }
/*  error codes and messages. */
/*<  600  ier = 3 >*/
L600:
    *ier = 3;
/*<       go to 660 >*/
    goto L660;
/*<  610  ier = 2 >*/
L610:
    *ier = 2;
/*<       go to 660 >*/
    goto L660;
/*<  620  ier = 1 >*/
L620:
    *ier = 1;
/*<       go to 660 >*/
    goto L660;
/*<  630  ier = -1 >*/
L630:
    *ier = -1;
/*<       go to 660 >*/
    goto L660;
/*<  640  ier = -2 >*/
L640:
    *ier = -2;
/*  the least-squares constant function c1 is a solution of our problem. */
/*  a constant function is a spline of degree k with all b-spline */
/*  coefficients equal to that constant c1. */
/*<       do 650 i=1,k1 >*/
    i__1 = *k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         rn = k1-i >*/
	rn = (doublereal) (*k1 - i__);
/*<         t(i) = x(1)-rn*per >*/
	t[i__] = x[1] - rn * per;
/*<         c(i) = c1 >*/
	c__[i__] = c1;
/*<         j = i+k1 >*/
	j = i__ + *k1;
/*<         rn = i-1 >*/
	rn = (doublereal) (i__ - 1);
/*<         t(j) = x(m)+rn*per >*/
	t[j] = x[*m] + rn * per;
/*<  650  continue >*/
/* L650: */
    }
/*<       n = nmin >*/
    *n = nmin;
/*<       fp = fp0 >*/
    *fp = fp0;
/*<       fpint(n) = fp0 >*/
    fpint[*n] = fp0;
/*<       fpint(n-1) = 0. >*/
    fpint[*n - 1] = 0.;
/*<       nrdata(n) = 0 >*/
    nrdata[*n] = 0;
/*<  660  return >*/
L660:
    return 0;
/*<       end >*/
} /* fpperi_ */

