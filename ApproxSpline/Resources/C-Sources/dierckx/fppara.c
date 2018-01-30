/* fppara.f -- translated by f2c (version 20061008).
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
/* Subroutine */ int fppara_(integer *iopt, integer *idim, integer *m, 
	doublereal *u, integer *mx, doublereal *x, doublereal *w, doublereal *
	ub, doublereal *ue, integer *k, doublereal *s, integer *nest, 
	doublereal *tol, integer *maxit, integer *k1, integer *k2, integer *n,
	 doublereal *t, integer *nc, doublereal *c__, doublereal *fp, 
	doublereal *fpint, doublereal *z__, doublereal *a, doublereal *b, 
	doublereal *g, doublereal *q, integer *nrdata, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static doublereal h__[7];
    static integer i__, j, l;
    static doublereal p, f1, f2, f3;
    static integer i1, i2, i3, j1, j2;
    static doublereal p1, p2, p3;
    static integer k3, l0, n8, jj, it;
    static doublereal ui, rn, wi, xi[10], fp0;
    static integer mk1, nk1;
    static doublereal acc, fac, one, cos__, sin__;
    static integer new__;
    static doublereal piv;
    static integer ich1, ich3;
    static doublereal con1, con4, con9;
    static integer npl1;
    static doublereal half;
    static integer nmin, iter, nmax;
    static doublereal fpms, term, pinv, fpold, fpart;
    static integer nrint;
    static doublereal store;
    static integer nplus;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), fpdisc_(doublereal *, 
	    integer *, integer *, doublereal *, integer *);
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
/*<       real ub,ue,s,tol,fp >*/
/*<       integer iopt,idim,m,mx,k,nest,maxit,k1,k2,n,nc,ier >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer nrdata(nest) >*/
/*  ..local scalars.. */
/*<    >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real h(7),xi(10) >*/
/*  ..function references */
/*<       real abs,fprati >*/
/*<       integer max0,min0 >*/
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota */
/*  .. */
/*  set constants */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --w;
    --u;
    --x;
    --nrdata;
    --fpint;
    --t;
    q_dim1 = *m;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *nest;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    g_dim1 = *nest;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    b_dim1 = *nest;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --z__;
    --c__;

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
/*  given a set of knots we compute the least-squares curve sinf(u),    c */
/*  and the corresponding sum of squared residuals fp=f(p=inf).         c */
/*  if iopt=-1 sinf(u) is the requested curve.                          c */
/*  if iopt=0 or iopt=1 we check whether we can accept the knots:       c */
/*    if fp <=s we will continue with the current set of knots.         c */
/*    if fp > s we will increase the number of knots and compute the    c */
/*       corresponding least-squares curve until finally fp<=s.         c */
/*    the initial choice of knots depends on the value of s and iopt.   c */
/*    if s=0 we have spline interpolation; in that case the number of   c */
/*    knots equals nmax = m+k+1.                                        c */
/*    if s > 0 and                                                      c */
/*      iopt=0 we first compute the least-squares polynomial curve of   c */
/*      degree k; n = nmin = 2*k+2                                      c */
/*      iopt=1 we start with the set of knots found at the last         c */
/*      call of the routine, except for the case that s > fp0; then     c */
/*      we compute directly the polynomial curve of degree k.           c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  determine nmin, the number of knots for polynomial approximation. */
/*<       nmin = 2*k1 >*/
    nmin = *k1 << 1;
/*<       if(iopt.lt.0) go to 60 >*/
    if (*iopt < 0) {
	goto L60;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
/*<       acc = tol*s >*/
    acc = *tol * *s;
/*  determine nmax, the number of knots for spline interpolation. */
/*<       nmax = m+k1 >*/
    nmax = *m + *k1;
/*<       if(s.gt.0.) go to 45 >*/
    if (*s > 0.) {
	goto L45;
    }
/*  if s=0, s(u) is an interpolating curve. */
/*  test whether the required storage space exceeds the available one. */
/*<       n = nmax >*/
    *n = nmax;
/*<       if(nmax.gt.nest) go to 420 >*/
    if (nmax > *nest) {
	goto L420;
    }
/*  find the position of the interior knots in case of interpolation. */
/*<   10  mk1 = m-k1 >*/
L10:
    mk1 = *m - *k1;
/*<       if(mk1.eq.0) go to 60 >*/
    if (mk1 == 0) {
	goto L60;
    }
/*<       k3 = k/2 >*/
    k3 = *k / 2;
/*<       i = k2 >*/
    i__ = *k2;
/*<       j = k3+2 >*/
    j = k3 + 2;
/*<       if(k3*2.eq.k) go to 30 >*/
    if (k3 << 1 == *k) {
	goto L30;
    }
/*<       do 20 l=1,mk1 >*/
    i__1 = mk1;
    for (l = 1; l <= i__1; ++l) {
/*<         t(i) = u(j) >*/
	t[i__] = u[j];
/*<         i = i+1 >*/
	++i__;
/*<         j = j+1 >*/
	++j;
/*<   20  continue >*/
/* L20: */
    }
/*<       go to 60 >*/
    goto L60;
/*<   30  do 40 l=1,mk1 >*/
L30:
    i__1 = mk1;
    for (l = 1; l <= i__1; ++l) {
/*<         t(i) = (u(j)+u(j-1))*half >*/
	t[i__] = (u[j] + u[j - 1]) * half;
/*<         i = i+1 >*/
	++i__;
/*<         j = j+1 >*/
	++j;
/*<   40  continue >*/
/* L40: */
    }
/*<       go to 60 >*/
    goto L60;
/*  if s>0 our initial choice of knots depends on the value of iopt. */
/*  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares */
/*  polynomial curve which is a spline curve without interior knots. */
/*  if iopt=1 and fp0>s we start computing the least squares spline curve */
/*  according to the set of knots found at the last call of the routine. */
/*<   45  if(iopt.eq.0) go to 50 >*/
L45:
    if (*iopt == 0) {
	goto L50;
    }
/*<       if(n.eq.nmin) go to 50 >*/
    if (*n == nmin) {
	goto L50;
    }
/*<       fp0 = fpint(n) >*/
    fp0 = fpint[*n];
/*<       fpold = fpint(n-1) >*/
    fpold = fpint[*n - 1];
/*<       nplus = nrdata(n) >*/
    nplus = nrdata[*n];
/*<       if(fp0.gt.s) go to 60 >*/
    if (fp0 > *s) {
	goto L60;
    }
/*<   50  n = nmin >*/
L50:
    *n = nmin;
/*<       fpold = 0. >*/
    fpold = 0.;
/*<       nplus = 0 >*/
    nplus = 0;
/*<       nrdata(1) = m-2 >*/
    nrdata[1] = *m - 2;
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
/*<   60  do 200 iter = 1,m >*/
L60:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*<         if(n.eq.nmin) ier = -2 >*/
	if (*n == nmin) {
	    *ier = -2;
	}
/*  find nrint, tne number of knot intervals. */
/*<         nrint = n-nmin+1 >*/
	nrint = *n - nmin + 1;
/*  find the position of the additional knots which are needed for */
/*  the b-spline representation of s(u). */
/*<         nk1 = n-k1 >*/
	nk1 = *n - *k1;
/*<         i = n >*/
	i__ = *n;
/*<         do 70 j=1,k1 >*/
	i__2 = *k1;
	for (j = 1; j <= i__2; ++j) {
/*<           t(j) = ub >*/
	    t[j] = *ub;
/*<           t(i) = ue >*/
	    t[i__] = *ue;
/*<           i = i-1 >*/
	    --i__;
/*<   70    continue >*/
/* L70: */
	}
/*  compute the b-spline coefficients of the least-squares spline curve */
/*  sinf(u). the observation matrix a is built up row by row and */
/*  reduced to upper triangular form by givens transformations. */
/*  at the same time fp=f(p=inf) is computed. */
/*<         fp = 0. >*/
	*fp = 0.;
/*  initialize the b-spline coefficients and the observation matrix a. */
/*<         do 75 i=1,nc >*/
	i__2 = *nc;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           z(i) = 0. >*/
	    z__[i__] = 0.;
/*<   75    continue >*/
/* L75: */
	}
/*<         do 80 i=1,nk1 >*/
	i__2 = nk1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           do 80 j=1,k1 >*/
	    i__3 = *k1;
	    for (j = 1; j <= i__3; ++j) {
/*<             a(i,j) = 0. >*/
		a[i__ + j * a_dim1] = 0.;
/*<   80    continue >*/
/* L80: */
	    }
	}
/*<         l = k1 >*/
	l = *k1;
/*<         jj = 0 >*/
	jj = 0;
/*<         do 130 it=1,m >*/
	i__3 = *m;
	for (it = 1; it <= i__3; ++it) {
/*  fetch the current data point u(it),x(it). */
/*<           ui = u(it) >*/
	    ui = u[it];
/*<           wi = w(it) >*/
	    wi = w[it];
/*<           do 83 j=1,idim >*/
	    i__2 = *idim;
	    for (j = 1; j <= i__2; ++j) {
/*<              jj = jj+1 >*/
		++jj;
/*<              xi(j) = x(jj)*wi >*/
		xi[j - 1] = x[jj] * wi;
/*<   83      continue >*/
/* L83: */
	    }
/*  search for knot interval t(l) <= ui < t(l+1). */
/*<   85      if(ui.lt.t(l+1) .or. l.eq.nk1) go to 90 >*/
L85:
	    if (ui < t[l + 1] || l == nk1) {
		goto L90;
	    }
/*<           l = l+1 >*/
	    ++l;
/*<           go to 85 >*/
	    goto L85;
/*  evaluate the (k+1) non-zero b-splines at ui and store them in q. */
/*<   90      call fpbspl(t,n,k,ui,l,h) >*/
L90:
	    fpbspl_(&t[1], n, k, &ui, &l, h__);
/*<           do 95 i=1,k1 >*/
	    i__2 = *k1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             q(it,i) = h(i) >*/
		q[it + i__ * q_dim1] = h__[i__ - 1];
/*<             h(i) = h(i)*wi >*/
		h__[i__ - 1] *= wi;
/*<   95      continue >*/
/* L95: */
	    }
/*  rotate the new row of the observation matrix into triangle. */
/*<           j = l-k1 >*/
	    j = l - *k1;
/*<           do 110 i=1,k1 >*/
	    i__2 = *k1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<             j = j+1 >*/
		++j;
/*<             piv = h(i) >*/
		piv = h__[i__ - 1];
/*<             if(piv.eq.0.) go to 110 >*/
		if (piv == 0.) {
		    goto L110;
		}
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a(j,1),cos,sin) >*/
		fpgivs_(&piv, &a[j + a_dim1], &cos__, &sin__);
/*  transformations to right hand side. */
/*<             j1 = j >*/
		j1 = j;
/*<             do 97 j2 =1,idim >*/
		i__4 = *idim;
		for (j2 = 1; j2 <= i__4; ++j2) {
/*<                call fprota(cos,sin,xi(j2),z(j1)) >*/
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &z__[j1]);
/*<                j1 = j1+n >*/
		    j1 += *n;
/*<   97        continue >*/
/* L97: */
		}
/*<             if(i.eq.k1) go to 120 >*/
		if (i__ == *k1) {
		    goto L120;
		}
/*<             i2 = 1 >*/
		i2 = 1;
/*<             i3 = i+1 >*/
		i3 = i__ + 1;
/*<             do 100 i1 = i3,k1 >*/
		i__4 = *k1;
		for (i1 = i3; i1 <= i__4; ++i1) {
/*<               i2 = i2+1 >*/
		    ++i2;
/*  transformations to left hand side. */
/*<               call fprota(cos,sin,h(i1),a(j,i2)) >*/
		    fprota_(&cos__, &sin__, &h__[i1 - 1], &a[j + i2 * a_dim1])
			    ;
/*<  100        continue >*/
/* L100: */
		}
/*<  110      continue >*/
L110:
		;
	    }
/*  add contribution of this row to the sum of squares of residual */
/*  right hand sides. */
/*<  120      do 125 j2=1,idim >*/
L120:
	    i__2 = *idim;
	    for (j2 = 1; j2 <= i__2; ++j2) {
/*<              fp  = fp+xi(j2)**2 >*/
/* Computing 2nd power */
		d__1 = xi[j2 - 1];
		*fp += d__1 * d__1;
/*<  125      continue >*/
/* L125: */
	    }
/*<  130    continue >*/
/* L130: */
	}
/*<         if(ier.eq.(-2)) fp0 = fp >*/
	if (*ier == -2) {
	    fp0 = *fp;
	}
/*<         fpint(n) = fp0 >*/
	fpint[*n] = fp0;
/*<         fpint(n-1) = fpold >*/
	fpint[*n - 1] = fpold;
/*<         nrdata(n) = nplus >*/
	nrdata[*n] = nplus;
/*  backward substitution to obtain the b-spline coefficients. */
/*<         j1 = 1 >*/
	j1 = 1;
/*<         do 135 j2=1,idim >*/
	i__3 = *idim;
	for (j2 = 1; j2 <= i__3; ++j2) {
/*<            call fpback(a,z(j1),nk1,k1,c(j1),nest) >*/
	    fpback_(&a[a_offset], &z__[j1], &nk1, k1, &c__[j1], nest);
/*<            j1 = j1+n >*/
	    j1 += *n;
/*<  135    continue >*/
/* L135: */
	}
/*  test whether the approximation sinf(u) is an acceptable solution. */
/*<         if(iopt.lt.0) go to 440 >*/
	if (*iopt < 0) {
	    goto L440;
	}
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).lt.acc) go to 440 >*/
	if (abs(fpms) < acc) {
	    goto L440;
	}
/*  if f(p=inf) < s accept the choice of knots. */
/*<         if(fpms.lt.0.) go to 250 >*/
	if (fpms < 0.) {
	    goto L250;
	}
/*  if n = nmax, sinf(u) is an interpolating spline curve. */
/*<         if(n.eq.nmax) go to 430 >*/
	if (*n == nmax) {
	    goto L430;
	}
/*  increase the number of knots. */
/*  if n=nest we cannot increase the number of knots because of */
/*  the storage capacity limitation. */
/*<         if(n.eq.nest) go to 420 >*/
	if (*n == *nest) {
	    goto L420;
	}
/*  determine the number of knots nplus we are going to add. */
/*<         if(ier.eq.0) go to 140 >*/
	if (*ier == 0) {
	    goto L140;
	}
/*<         nplus = 1 >*/
	nplus = 1;
/*<         ier = 0 >*/
	*ier = 0;
/*<         go to 150 >*/
	goto L150;
/*<  140    npl1 = nplus*2 >*/
L140:
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
	i__4 = npl1, i__5 = nplus / 2, i__4 = max(i__4,i__5);
	i__3 = nplus << 1, i__2 = max(i__4,1);
	nplus = min(i__3,i__2);
/*<  150    fpold = fp >*/
L150:
	fpold = *fp;
/*  compute the sum of squared residuals for each knot interval */
/*  t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
/*<         fpart = 0. >*/
	fpart = 0.;
/*<         i = 1 >*/
	i__ = 1;
/*<         l = k2 >*/
	l = *k2;
/*<         new = 0 >*/
	new__ = 0;
/*<         jj = 0 >*/
	jj = 0;
/*<         do 180 it=1,m >*/
	i__3 = *m;
	for (it = 1; it <= i__3; ++it) {
/*<           if(u(it).lt.t(l) .or. l.gt.nk1) go to 160 >*/
	    if (u[it] < t[l] || l > nk1) {
		goto L160;
	    }
/*<           new = 1 >*/
	    new__ = 1;
/*<           l = l+1 >*/
	    ++l;
/*<  160      term = 0. >*/
L160:
	    term = 0.;
/*<           l0 = l-k2 >*/
	    l0 = l - *k2;
/*<           do 175 j2=1,idim >*/
	    i__2 = *idim;
	    for (j2 = 1; j2 <= i__2; ++j2) {
/*<             fac = 0. >*/
		fac = 0.;
/*<             j1 = l0 >*/
		j1 = l0;
/*<             do 170 j=1,k1 >*/
		i__4 = *k1;
		for (j = 1; j <= i__4; ++j) {
/*<               j1 = j1+1 >*/
		    ++j1;
/*<               fac = fac+c(j1)*q(it,j) >*/
		    fac += c__[j1] * q[it + j * q_dim1];
/*<  170        continue >*/
/* L170: */
		}
/*<             jj = jj+1 >*/
		++jj;
/*<             term = term+(w(it)*(fac-x(jj)))**2 >*/
/* Computing 2nd power */
		d__1 = w[it] * (fac - x[jj]);
		term += d__1 * d__1;
/*<             l0 = l0+n >*/
		l0 += *n;
/*<  175      continue >*/
/* L175: */
	    }
/*<           fpart = fpart+term >*/
	    fpart += term;
/*<           if(new.eq.0) go to 180 >*/
	    if (new__ == 0) {
		goto L180;
	    }
/*<           store = term*half >*/
	    store = term * half;
/*<           fpint(i) = fpart-store >*/
	    fpint[i__] = fpart - store;
/*<           i = i+1 >*/
	    ++i__;
/*<           fpart = store >*/
	    fpart = store;
/*<           new = 0 >*/
	    new__ = 0;
/*<  180    continue >*/
L180:
	    ;
	}
/*<         fpint(nrint) = fpart >*/
	fpint[nrint] = fpart;
/*<         do 190 l=1,nplus >*/
	i__3 = nplus;
	for (l = 1; l <= i__3; ++l) {
/*  add a new knot. */
/*<           call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1) >*/
	    fpknot_(&u[1], m, &t[1], n, &fpint[1], &nrdata[1], &nrint, nest, &
		    c__1);
/*  if n=nmax we locate the knots as for interpolation */
/*<           if(n.eq.nmax) go to 10 >*/
	    if (*n == nmax) {
		goto L10;
	    }
/*  test whether we cannot further increase the number of knots. */
/*<           if(n.eq.nest) go to 200 >*/
	    if (*n == *nest) {
		goto L200;
	    }
/*<  190    continue >*/
/* L190: */
	}
/*  restart the computations with the new set of knots. */
/*<  200  continue >*/
L200:
	;
    }
/*  test whether the least-squares kth degree polynomial curve is a */
/*  solution of our approximation problem. */
/*<  250  if(ier.eq.(-2)) go to 440 >*/
L250:
    if (*ier == -2) {
	goto L440;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 2: determination of the smoothing spline curve sp(u).          c */
/*  **********************************************************          c */
/*  we have determined the number of knots and their position.          c */
/*  we now compute the b-spline coefficients of the smoothing curve     c */
/*  sp(u). the observation matrix a is extended by the rows of matrix   c */
/*  b expressing that the kth derivative discontinuities of sp(u) at    c */
/*  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c */
/*  ponding weights of these additional rows are set to 1/p.            c */
/*  iteratively we then have to determine the value of p such that f(p),c */
/*  the sum of squared residuals be = s. we already know that the least c */
/*  squares kth degree polynomial curve corresponds to p=0, and that    c */
/*  the least-squares spline curve corresponds to p=infinity. the       c */
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
/*<       call fpdisc(t,n,k2,b,nest) >*/
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
/*<       p = 0. >*/
    p = 0.;
/*<       do 252 i=1,nk1 >*/
    i__1 = nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          p = p+a(i,1) >*/
	p += a[i__ + a_dim1];
/*<  252  continue >*/
/* L252: */
    }
/*<       rn = nk1 >*/
    rn = (doublereal) nk1;
/*<       p = rn/p >*/
    p = rn / p;
/*<       ich1 = 0 >*/
    ich1 = 0;
/*<       ich3 = 0 >*/
    ich3 = 0;
/*<       n8 = n-nmin >*/
    n8 = *n - nmin;
/*  iteration process to find the root of f(p) = s. */
/*<       do 360 iter=1,maxit >*/
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*  the rows of matrix b with weight 1/p are rotated into the */
/*  triangularised observation matrix a which is stored in g. */
/*<         pinv = one/p >*/
	pinv = one / p;
/*<         do 255 i=1,nc >*/
	i__3 = *nc;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           c(i) = z(i) >*/
	    c__[i__] = z__[i__];
/*<  255    continue >*/
/* L255: */
	}
/*<         do 260 i=1,nk1 >*/
	i__3 = nk1;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           g(i,k2) = 0. >*/
	    g[i__ + *k2 * g_dim1] = 0.;
/*<           do 260 j=1,k1 >*/
	    i__2 = *k1;
	    for (j = 1; j <= i__2; ++j) {
/*<             g(i,j) = a(i,j) >*/
		g[i__ + j * g_dim1] = a[i__ + j * a_dim1];
/*<  260    continue >*/
/* L260: */
	    }
	}
/*<         do 300 it=1,n8 >*/
	i__2 = n8;
	for (it = 1; it <= i__2; ++it) {
/*  the row of matrix b is rotated into triangle by givens transformation */
/*<           do 270 i=1,k2 >*/
	    i__3 = *k2;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             h(i) = b(it,i)*pinv >*/
		h__[i__ - 1] = b[it + i__ * b_dim1] * pinv;
/*<  270      continue >*/
/* L270: */
	    }
/*<           do 275 j=1,idim >*/
	    i__3 = *idim;
	    for (j = 1; j <= i__3; ++j) {
/*<             xi(j) = 0. >*/
		xi[j - 1] = 0.;
/*<  275      continue >*/
/* L275: */
	    }
/*<           do 290 j=it,nk1 >*/
	    i__3 = nk1;
	    for (j = it; j <= i__3; ++j) {
/*<             piv = h(1) >*/
		piv = h__[0];
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,g(j,1),cos,sin) >*/
		fpgivs_(&piv, &g[j + g_dim1], &cos__, &sin__);
/*  transformations to right hand side. */
/*<             j1 = j >*/
		j1 = j;
/*<             do 277 j2=1,idim >*/
		i__4 = *idim;
		for (j2 = 1; j2 <= i__4; ++j2) {
/*<               call fprota(cos,sin,xi(j2),c(j1)) >*/
		    fprota_(&cos__, &sin__, &xi[j2 - 1], &c__[j1]);
/*<               j1 = j1+n >*/
		    j1 += *n;
/*<  277        continue >*/
/* L277: */
		}
/*<             if(j.eq.nk1) go to 300 >*/
		if (j == nk1) {
		    goto L300;
		}
/*<             i2 = k1 >*/
		i2 = *k1;
/*<             if(j.gt.n8) i2 = nk1-j >*/
		if (j > n8) {
		    i2 = nk1 - j;
		}
/*<             do 280 i=1,i2 >*/
		i__4 = i2;
		for (i__ = 1; i__ <= i__4; ++i__) {
/*  transformations to left hand side. */
/*<               i1 = i+1 >*/
		    i1 = i__ + 1;
/*<               call fprota(cos,sin,h(i1),g(j,i1)) >*/
		    fprota_(&cos__, &sin__, &h__[i1 - 1], &g[j + i1 * g_dim1])
			    ;
/*<               h(i) = h(i1) >*/
		    h__[i__ - 1] = h__[i1 - 1];
/*<  280        continue >*/
/* L280: */
		}
/*<             h(i2+1) = 0. >*/
		h__[i2] = 0.;
/*<  290      continue >*/
/* L290: */
	    }
/*<  300    continue >*/
L300:
	    ;
	}
/*  backward substitution to obtain the b-spline coefficients. */
/*<         j1 = 1 >*/
	j1 = 1;
/*<         do 305 j2=1,idim >*/
	i__2 = *idim;
	for (j2 = 1; j2 <= i__2; ++j2) {
/*<           call fpback(g,c(j1),nk1,k2,c(j1),nest) >*/
	    fpback_(&g[g_offset], &c__[j1], &nk1, k2, &c__[j1], nest);
/*<           j1 =j1+n >*/
	    j1 += *n;
/*<  305    continue >*/
/* L305: */
	}
/*  computation of f(p). */
/*<         fp = 0. >*/
	*fp = 0.;
/*<         l = k2 >*/
	l = *k2;
/*<         jj = 0 >*/
	jj = 0;
/*<         do 330 it=1,m >*/
	i__2 = *m;
	for (it = 1; it <= i__2; ++it) {
/*<           if(u(it).lt.t(l) .or. l.gt.nk1) go to 310 >*/
	    if (u[it] < t[l] || l > nk1) {
		goto L310;
	    }
/*<           l = l+1 >*/
	    ++l;
/*<  310      l0 = l-k2 >*/
L310:
	    l0 = l - *k2;
/*<           term = 0. >*/
	    term = 0.;
/*<           do 325 j2=1,idim >*/
	    i__3 = *idim;
	    for (j2 = 1; j2 <= i__3; ++j2) {
/*<             fac = 0. >*/
		fac = 0.;
/*<             j1 = l0 >*/
		j1 = l0;
/*<             do 320 j=1,k1 >*/
		i__4 = *k1;
		for (j = 1; j <= i__4; ++j) {
/*<               j1 = j1+1 >*/
		    ++j1;
/*<               fac = fac+c(j1)*q(it,j) >*/
		    fac += c__[j1] * q[it + j * q_dim1];
/*<  320        continue >*/
/* L320: */
		}
/*<             jj = jj+1 >*/
		++jj;
/*<             term = term+(fac-x(jj))**2 >*/
/* Computing 2nd power */
		d__1 = fac - x[jj];
		term += d__1 * d__1;
/*<             l0 = l0+n >*/
		l0 += *n;
/*<  325      continue >*/
/* L325: */
	    }
/*<           fp = fp+term*w(it)**2 >*/
/* Computing 2nd power */
	    d__1 = w[it];
	    *fp += term * (d__1 * d__1);
/*<  330    continue >*/
/* L330: */
	}
/*  test whether the approximation sp(u) is an acceptable solution. */
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).lt.acc) go to 440 >*/
	if (abs(fpms) < acc) {
	    goto L440;
	}
/*  test whether the maximal number of iterations is reached. */
/*<         if(iter.eq.maxit) go to 400 >*/
	if (iter == *maxit) {
	    goto L400;
	}
/*  carry out one more step of the iteration process. */
/*<         p2 = p >*/
	p2 = p;
/*<         f2 = fpms >*/
	f2 = fpms;
/*<         if(ich3.ne.0) go to 340 >*/
	if (ich3 != 0) {
	    goto L340;
	}
/*<         if((f2-f3).gt.acc) go to 335 >*/
	if (f2 - f3 > acc) {
	    goto L335;
	}
/*  our initial choice of p is too large. */
/*<         p3 = p2 >*/
	p3 = p2;
/*<         f3 = f2 >*/
	f3 = f2;
/*<         p = p*con4 >*/
	p *= con4;
/*<         if(p.le.p1) p=p1*con9 + p2*con1 >*/
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
/*<         go to 360 >*/
	goto L360;
/*<  335    if(f2.lt.0.) ich3=1 >*/
L335:
	if (f2 < 0.) {
	    ich3 = 1;
	}
/*<  340    if(ich1.ne.0) go to 350 >*/
L340:
	if (ich1 != 0) {
	    goto L350;
	}
/*<         if((f1-f2).gt.acc) go to 345 >*/
	if (f1 - f2 > acc) {
	    goto L345;
	}
/*  our initial choice of p is too small */
/*<         p1 = p2 >*/
	p1 = p2;
/*<         f1 = f2 >*/
	f1 = f2;
/*<         p = p/con4 >*/
	p /= con4;
/*<         if(p3.lt.0.) go to 360 >*/
	if (p3 < 0.) {
	    goto L360;
	}
/*<         if(p.ge.p3) p = p2*con1 + p3*con9 >*/
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
/*<         go to 360 >*/
	goto L360;
/*<  345    if(f2.gt.0.) ich1=1 >*/
L345:
	if (f2 > 0.) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
/*<  350    if(f2.ge.f1 .or. f2.le.f3) go to 410 >*/
L350:
	if (f2 >= f1 || f2 <= f3) {
	    goto L410;
	}
/*  find the new value for p. */
/*<         p = fprati(p1,f1,p2,f2,p3,f3) >*/
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
/*<  360  continue >*/
L360:
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
/*<  440  return >*/
L440:
    return 0;
/*<       end >*/
} /* fppara_ */

