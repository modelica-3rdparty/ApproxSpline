/* fpcosp.f -- translated by f2c (version 20061008).
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

/*<    >*/
/* Subroutine */ int fpcosp_(integer *m, doublereal *x, doublereal *y, 
	doublereal *w, integer *n, doublereal *t, doublereal *e, integer *
	maxtr, integer *maxbin, doublereal *c__, doublereal *sq, doublereal *
	sx, logical *bind, integer *nm, integer *mb, doublereal *a, 
	doublereal *b, doublereal *const__, doublereal *z__, doublereal *zz, 
	doublereal *u, doublereal *q, integer *info, integer *up, integer *
	left, integer *right, integer *jbind, integer *ibind, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal f, h__[4];
    static integer i__, j, k, l, i1, j1, j2, j3, k1, k2, k3, k4, k5, k6, l1, 
	    l2, l3, n1, n4, n6;
    static doublereal wi, xi;
    static integer lp1, kdim, merk, nbind, count;
    extern /* Subroutine */ int fpadno_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), fpdeno_(integer *, integer *, integer *, integer *, 
	    integer *, integer *), fpbspl_(doublereal *, integer *, integer *,
	     doublereal *, integer *, doublereal *);
    static integer number;
    extern /* Subroutine */ int fpseno_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);

/*  .. */
/*  ..scalar arguments.. */
/*<       real sq >*/
/*<       integer m,n,maxtr,maxbin,nm,mb,ier >*/
/*  ..array arguments.. */
/*<    >*/
/*<    >*/
/*<       logical bind(n) >*/
/*  ..local scalars.. */
/*<    >*/
/*<       real f,wi,xi >*/
/*  ..local array.. */
/*<       real h(4) >*/
/*  ..subroutine references.. */
/*    fpbspl,fpadno,fpdeno,fpfrno,fpseno */
/*  .. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  if we use the b-spline representation of s(x) our approximation     c */
/*  problem results in a quadratic programming problem:                 c */
/*    find the b-spline coefficients c(j),j=1,2,...n-4 such that        c */
/*        (1) sumi((wi*(yi-sumj(cj*nj(xi))))**2),i=1,2,...m is minimal  c */
/*        (2) sumj(cj*n''j(t(l+3)))*e(l) <= 0, l=1,2,...n-6.            c */
/*  to solve this problem we use the theil-van de panne procedure.      c */
/*  if the inequality constraints (2) are numbered from 1 to n-6,       c */
/*  this algorithm finds a subset of constraints ibind(1)..ibind(nbind) c */
/*  such that the solution of the minimization problem (1) with these   c */
/*  constraints in equality form, satisfies all constraints. such a     c */
/*  feasible solution is optimal if the lagrange parameters associated  c */
/*  with that problem with equality constraints, are all positive.      c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  determine n6, the number of inequality constraints. */
/*<       n6 = n-6 >*/
    /* Parameter adjustments */
    q_dim1 = *m;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --sx;
    --w;
    --y;
    --x;
    --zz;
    --z__;
    --const__;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --bind;
    --c__;
    --e;
    --t;
    --right;
    --left;
    --up;
    --info;
    --u;
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --ibind;
    --jbind;

    /* Function Body */
    n6 = *n - 6;
/*  fix the parameters which determine these constraints. */
/*<       do 10 i=1,n6 >*/
    i__1 = n6;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         const(i) = e(i)*(t(i+4)-t(i+1))/(t(i+5)-t(i+2)) >*/
	const__[i__] = e[i__] * (t[i__ + 4] - t[i__ + 1]) / (t[i__ + 5] - t[
		i__ + 2]);
/*<   10  continue >*/
/* L10: */
    }
/*  initialize the triply linked tree which is used to find the subset */
/*  of constraints ibind(1),...ibind(nbind). */
/*<       count = 1 >*/
    count = 1;
/*<       info(1) = 0 >*/
    info[1] = 0;
/*<       left(1) = 0 >*/
    left[1] = 0;
/*<       right(1) = 0 >*/
    right[1] = 0;
/*<       up(1) = 1 >*/
    up[1] = 1;
/*<       merk = 1 >*/
    merk = 1;
/*  set up the normal equations  n'nc=n'y  where n denotes the m x (n-4) */
/*  observation matrix with elements ni,j = wi*nj(xi)  and y is the */
/*  column vector with elements yi*wi. */
/*  from the properties of the b-splines nj(x),j=1,2,...n-4, it follows */
/*  that  n'n  is a (n-4) x (n-4)  positive definit bandmatrix of */
/*  bandwidth 7. the matrices n'n and n'y are built up in a and z. */
/*<       n4 = n-4 >*/
    n4 = *n - 4;
/*  initialization */
/*<       do 20 i=1,n4 >*/
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         z(i) = 0. >*/
	z__[i__] = 0.;
/*<         do 20 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           a(i,j) = 0. >*/
	    a[i__ + j * a_dim1] = 0.;
/*<   20  continue >*/
/* L20: */
	}
    }
/*<       l = 4 >*/
    l = 4;
/*<       lp1 = l+1 >*/
    lp1 = l + 1;
/*<       do 70 i=1,m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch the current row of the observation matrix. */
/*<         xi = x(i) >*/
	xi = x[i__];
/*<         wi = w(i)**2 >*/
/* Computing 2nd power */
	d__1 = w[i__];
	wi = d__1 * d__1;
/*  search for knot interval  t(l) <= xi < t(l+1) */
/*<   30    if(xi.lt.t(lp1) .or. l.eq.n4) go to 40 >*/
L30:
	if (xi < t[lp1] || l == n4) {
	    goto L40;
	}
/*<         l = lp1 >*/
	l = lp1;
/*<         lp1 = l+1 >*/
	lp1 = l + 1;
/*<         go to 30 >*/
	goto L30;
/*  evaluate the four non-zero cubic b-splines nj(xi),j=l-3,...l. */
/*<   40    call fpbspl(t,n,3,xi,l,h) >*/
L40:
	fpbspl_(&t[1], n, &c__3, &xi, &l, h__);
/*  store in q these values h(1),h(2),...h(4). */
/*<         do 50 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           q(i,j) = h(j) >*/
	    q[i__ + j * q_dim1] = h__[j - 1];
/*<   50    continue >*/
/* L50: */
	}
/*  add the contribution of the current row of the observation matrix */
/*  n to the normal equations. */
/*<         l3 = l-3 >*/
	l3 = l - 3;
/*<         k1 = 0 >*/
	k1 = 0;
/*<         do 60 j1 = l3,l >*/
	i__2 = l;
	for (j1 = l3; j1 <= i__2; ++j1) {
/*<           k1 = k1+1 >*/
	    ++k1;
/*<           f = h(k1) >*/
	    f = h__[k1 - 1];
/*<           z(j1) = z(j1)+f*wi*y(i) >*/
	    z__[j1] += f * wi * y[i__];
/*<           k2 = k1 >*/
	    k2 = k1;
/*<           j2 = 4 >*/
	    j2 = 4;
/*<           do 60 j3 = j1,l >*/
	    i__3 = l;
	    for (j3 = j1; j3 <= i__3; ++j3) {
/*<             a(j3,j2) = a(j3,j2)+f*wi*h(k2) >*/
		a[j3 + j2 * a_dim1] += f * wi * h__[k2 - 1];
/*<             k2 = k2+1 >*/
		++k2;
/*<             j2 = j2-1 >*/
		--j2;
/*<   60    continue >*/
/* L60: */
	    }
	}
/*<   70  continue >*/
/* L70: */
    }
/*  since n'n is a symmetric matrix it can be factorized as */
/*        (3)  n'n = (r1)'(d1)(r1) */
/*  with d1 a diagonal matrix and r1 an (n-4) x (n-4)  unit upper */
/*  triangular matrix of bandwidth 4. the matrices r1 and d1 are built */
/*  up in a. at the same time we solve the systems of equations */
/*        (4)  (r1)'(z2) = n'y */
/*        (5)  (d1) (z1) = (z2) */
/*  the vectors z2 and z1 are kept in zz and z. */
/*<       do 140 i=1,n4 >*/
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         k1 = 1 >*/
	k1 = 1;
/*<         if(i.lt.4) k1 = 5-i >*/
	if (i__ < 4) {
	    k1 = 5 - i__;
	}
/*<         k2 = i-4+k1 >*/
	k2 = i__ - 4 + k1;
/*<         k3 = k2 >*/
	k3 = k2;
/*<         do 100 j=k1,4 >*/
	for (j = k1; j <= 4; ++j) {
/*<           k4 = j-1 >*/
	    k4 = j - 1;
/*<           k5 = 4-j+k1 >*/
	    k5 = 4 - j + k1;
/*<           f = a(i,j) >*/
	    f = a[i__ + j * a_dim1];
/*<           if(k1.gt.k4) go to 90 >*/
	    if (k1 > k4) {
		goto L90;
	    }
/*<           k6 = k2 >*/
	    k6 = k2;
/*<           do 80 k=k1,k4 >*/
	    i__3 = k4;
	    for (k = k1; k <= i__3; ++k) {
/*<             f = f-a(i,k)*a(k3,k5)*a(k6,4) >*/
		f -= a[i__ + k * a_dim1] * a[k3 + k5 * a_dim1] * a[k6 + (
			a_dim1 << 2)];
/*<             k5 = k5+1 >*/
		++k5;
/*<             k6 = k6+1 >*/
		++k6;
/*<   80      continue >*/
/* L80: */
	    }
/*<   90      if(j.eq.4) go to 110 >*/
L90:
	    if (j == 4) {
		goto L110;
	    }
/*<           a(i,j) = f/a(k3,4) >*/
	    a[i__ + j * a_dim1] = f / a[k3 + (a_dim1 << 2)];
/*<           k3 = k3+1 >*/
	    ++k3;
/*<  100    continue >*/
/* L100: */
	}
/*<  110    a(i,4) = f >*/
L110:
	a[i__ + (a_dim1 << 2)] = f;
/*<         f = z(i) >*/
	f = z__[i__];
/*<         if(i.eq.1) go to 130 >*/
	if (i__ == 1) {
	    goto L130;
	}
/*<         k4 = i >*/
	k4 = i__;
/*<         do 120 j=k1,3 >*/
	for (j = k1; j <= 3; ++j) {
/*<           k = k1+3-j >*/
	    k = k1 + 3 - j;
/*<           k4 = k4-1 >*/
	    --k4;
/*<           f = f-a(i,k)*z(k4)*a(k4,4) >*/
	    f -= a[i__ + k * a_dim1] * z__[k4] * a[k4 + (a_dim1 << 2)];
/*<  120    continue >*/
/* L120: */
	}
/*<  130    z(i) = f/a(i,4) >*/
L130:
	z__[i__] = f / a[i__ + (a_dim1 << 2)];
/*<         zz(i) = f >*/
	zz[i__] = f;
/*<  140  continue >*/
/* L140: */
    }
/*  start computing the least-squares cubic spline without taking account */
/*  of any constraint. */
/*<       nbind = 0 >*/
    nbind = 0;
/*<       n1 = 1 >*/
    n1 = 1;
/*<       ibind(1) = 0 >*/
    ibind[1] = 0;
/*  main loop for the least-squares problems with different subsets of */
/*  the constraints (2) in equality form. the resulting b-spline coeff. */
/*  c and lagrange parameters u are the solution of the system */
/*            ! n'n  b' ! ! c !   ! n'y ! */
/*        (6) !         ! !   ! = !     ! */
/*            !  b   0  ! ! u !   !  0  ! */
/*  z1 is stored into array c. */
/*<  150  do 160 i=1,n4 >*/
L150:
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         c(i) = z(i) >*/
	c__[i__] = z__[i__];
/*<  160  continue >*/
/* L160: */
    }
/*  if there are no equality constraints, compute the coeff. c directly. */
/*<       if(nbind.eq.0) go to 370 >*/
    if (nbind == 0) {
	goto L370;
    }
/*  initialization */
/*<       kdim = n4+nbind >*/
    kdim = n4 + nbind;
/*<       do 170 i=1,nbind >*/
    i__1 = nbind;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         do 170 j=1,kdim >*/
	i__3 = kdim;
	for (j = 1; j <= i__3; ++j) {
/*<           b(j,i) = 0. >*/
	    b[j + i__ * b_dim1] = 0.;
/*<  170  continue >*/
/* L170: */
	}
    }
/*  matrix b is built up,expressing that the constraints nrs ibind(1),... */
/*  ibind(nbind) must be satisfied in equality form. */
/*<       do 180 i=1,nbind >*/
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         l = ibind(i) >*/
	l = ibind[i__];
/*<         b(l,i) = e(l) >*/
	b[l + i__ * b_dim1] = e[l];
/*<         b(l+1,i) = -(e(l)+const(l)) >*/
	b[l + 1 + i__ * b_dim1] = -(e[l] + const__[l]);
/*<         b(l+2,i) = const(l) >*/
	b[l + 2 + i__ * b_dim1] = const__[l];
/*<  180  continue >*/
/* L180: */
    }
/*  find the matrix (b1) as the solution of the system of equations */
/*        (7)  (r1)'(d1)(b1) = b' */
/*  (b1) is built up in the upper part of the array b(rows 1,...n-4). */
/*<       do 220 k1=1,nbind >*/
    i__3 = nbind;
    for (k1 = 1; k1 <= i__3; ++k1) {
/*<         l = ibind(k1) >*/
	l = ibind[k1];
/*<         do 210 i=l,n4 >*/
	i__1 = n4;
	for (i__ = l; i__ <= i__1; ++i__) {
/*<           f = b(i,k1) >*/
	    f = b[i__ + k1 * b_dim1];
/*<           if(i.eq.1) go to 200 >*/
	    if (i__ == 1) {
		goto L200;
	    }
/*<           k2 = 3 >*/
	    k2 = 3;
/*<           if(i.lt.4) k2 = i-1 >*/
	    if (i__ < 4) {
		k2 = i__ - 1;
	    }
/*<           do 190 k3=1,k2 >*/
	    i__2 = k2;
	    for (k3 = 1; k3 <= i__2; ++k3) {
/*<             l1 = i-k3 >*/
		l1 = i__ - k3;
/*<             l2 = 4-k3 >*/
		l2 = 4 - k3;
/*<             f = f-b(l1,k1)*a(i,l2)*a(l1,4) >*/
		f -= b[l1 + k1 * b_dim1] * a[i__ + l2 * a_dim1] * a[l1 + (
			a_dim1 << 2)];
/*<  190      continue >*/
/* L190: */
	    }
/*<  200      b(i,k1) = f/a(i,4) >*/
L200:
	    b[i__ + k1 * b_dim1] = f / a[i__ + (a_dim1 << 2)];
/*<  210    continue >*/
/* L210: */
	}
/*<  220  continue >*/
/* L220: */
    }
/*  factorization of the symmetric matrix  -(b1)'(d1)(b1) */
/*        (8)  -(b1)'(d1)(b1) = (r2)'(d2)(r2) */
/*  with (d2) a diagonal matrix and (r2) an nbind x nbind unit upper */
/*  triangular matrix. the matrices r2 and d2 are built up in the lower */
/*  part of the array b (rows n-3,n-2,...n-4+nbind). */
/*<       do 270 i=1,nbind >*/
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         i1 = i-1 >*/
	i1 = i__ - 1;
/*<         do 260 j=i,nbind >*/
	i__1 = nbind;
	for (j = i__; j <= i__1; ++j) {
/*<           f = 0. >*/
	    f = 0.;
/*<           do 230 k=1,n4 >*/
	    i__2 = n4;
	    for (k = 1; k <= i__2; ++k) {
/*<             f = f+b(k,i)*b(k,j)*a(k,4) >*/
		f += b[k + i__ * b_dim1] * b[k + j * b_dim1] * a[k + (a_dim1 
			<< 2)];
/*<  230      continue >*/
/* L230: */
	    }
/*<           k1 = n4+1 >*/
	    k1 = n4 + 1;
/*<           if(i1.eq.0) go to 250 >*/
	    if (i1 == 0) {
		goto L250;
	    }
/*<           do 240 k=1,i1 >*/
	    i__2 = i1;
	    for (k = 1; k <= i__2; ++k) {
/*<             f = f+b(k1,i)*b(k1,j)*b(k1,k) >*/
		f += b[k1 + i__ * b_dim1] * b[k1 + j * b_dim1] * b[k1 + k * 
			b_dim1];
/*<             k1 = k1+1 >*/
		++k1;
/*<  240      continue >*/
/* L240: */
	    }
/*<  250      b(k1,j) = -f >*/
L250:
	    b[k1 + j * b_dim1] = -f;
/*<           if(j.eq.i) go to 260 >*/
	    if (j == i__) {
		goto L260;
	    }
/*<           b(k1,j) = b(k1,j)/b(k1,i) >*/
	    b[k1 + j * b_dim1] /= b[k1 + i__ * b_dim1];
/*<  260    continue >*/
L260:
	    ;
	}
/*<  270  continue >*/
/* L270: */
    }
/*  according to (3),(7) and (8) the system of equations (6) becomes */
/*         ! (r1)'    0  ! ! (d1)    0  ! ! (r1)  (b1) ! ! c !   ! n'y ! */
/*    (9)  !             ! !            ! !            ! !   ! = !     ! */
/*         ! (b1)'  (r2)'! !   0   (d2) ! !   0   (r2) ! ! u !   !  0  ! */
/*  backward substitution to obtain the b-spline coefficients c(j),j=1,.. */
/*  n-4 and the lagrange parameters u(j),j=1,2,...nbind. */
/*  first step of the backward substitution: solve the system */
/*             ! (r1)'(d1)      0     ! ! (c1) !   ! n'y ! */
/*        (10) !                      ! !      ! = !     ! */
/*             ! (b1)'(d1)  (r2)'(d2) ! ! (u1) !   !  0  ! */
/*  from (4) and (5) we know that this is equivalent to */
/*        (11)  (c1) = (z1) */
/*        (12)  (r2)'(d2)(u1) = -(b1)'(z2) */
/*<       do 310 i=1,nbind >*/
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         f = 0. >*/
	f = 0.;
/*<         do 280 j=1,n4 >*/
	i__1 = n4;
	for (j = 1; j <= i__1; ++j) {
/*<           f = f+b(j,i)*zz(j) >*/
	    f += b[j + i__ * b_dim1] * zz[j];
/*<  280    continue >*/
/* L280: */
	}
/*<         i1 = i-1 >*/
	i1 = i__ - 1;
/*<         k1 = n4+1 >*/
	k1 = n4 + 1;
/*<         if(i1.eq.0) go to 300 >*/
	if (i1 == 0) {
	    goto L300;
	}
/*<         do 290 j=1,i1 >*/
	i__1 = i1;
	for (j = 1; j <= i__1; ++j) {
/*<           f = f+u(j)*b(k1,i)*b(k1,j) >*/
	    f += u[j] * b[k1 + i__ * b_dim1] * b[k1 + j * b_dim1];
/*<           k1 = k1+1 >*/
	    ++k1;
/*<  290    continue >*/
/* L290: */
	}
/*<  300    u(i) = -f/b(k1,i) >*/
L300:
	u[i__] = -f / b[k1 + i__ * b_dim1];
/*<  310  continue >*/
/* L310: */
    }
/*  second step of the backward substitution: solve the system */
/*             ! (r1)  (b1) ! ! c !   ! c1 ! */
/*        (13) !            ! !   ! = !    ! */
/*             !   0   (r2) ! ! u !   ! u1 ! */
/*<       k1 = nbind >*/
    k1 = nbind;
/*<       k2 = kdim >*/
    k2 = kdim;
/*  find the lagrange parameters u. */
/*<       do 340 i=1,nbind >*/
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         f = u(k1) >*/
	f = u[k1];
/*<         if(i.eq.1) go to 330 >*/
	if (i__ == 1) {
	    goto L330;
	}
/*<         k3 = k1+1 >*/
	k3 = k1 + 1;
/*<         do 320 j=k3,nbind >*/
	i__1 = nbind;
	for (j = k3; j <= i__1; ++j) {
/*<           f = f-u(j)*b(k2,j) >*/
	    f -= u[j] * b[k2 + j * b_dim1];
/*<  320    continue >*/
/* L320: */
	}
/*<  330    u(k1) = f >*/
L330:
	u[k1] = f;
/*<         k1 = k1-1 >*/
	--k1;
/*<         k2 = k2-1 >*/
	--k2;
/*<  340  continue >*/
/* L340: */
    }
/*  find the b-spline coefficients c. */
/*<       do 360 i=1,n4 >*/
    i__3 = n4;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         f = c(i) >*/
	f = c__[i__];
/*<         do 350 j=1,nbind >*/
	i__1 = nbind;
	for (j = 1; j <= i__1; ++j) {
/*<           f = f-u(j)*b(i,j) >*/
	    f -= u[j] * b[i__ + j * b_dim1];
/*<  350    continue >*/
/* L350: */
	}
/*<         c(i) = f >*/
	c__[i__] = f;
/*<  360  continue >*/
/* L360: */
    }
/*<  370  k1 = n4 >*/
L370:
    k1 = n4;
/*<       do 390 i=2,n4 >*/
    i__3 = n4;
    for (i__ = 2; i__ <= i__3; ++i__) {
/*<         k1 = k1-1 >*/
	--k1;
/*<         f = c(k1) >*/
	f = c__[k1];
/*<         k2 = 1 >*/
	k2 = 1;
/*<         if(i.lt.5) k2 = 5-i >*/
	if (i__ < 5) {
	    k2 = 5 - i__;
	}
/*<         k3 = k1 >*/
	k3 = k1;
/*<         l = 3 >*/
	l = 3;
/*<         do 380 j=k2,3 >*/
	for (j = k2; j <= 3; ++j) {
/*<           k3 = k3+1 >*/
	    ++k3;
/*<           f = f-a(k3,l)*c(k3) >*/
	    f -= a[k3 + l * a_dim1] * c__[k3];
/*<           l = l-1 >*/
	    --l;
/*<  380    continue >*/
/* L380: */
	}
/*<         c(k1) = f >*/
	c__[k1] = f;
/*<  390  continue >*/
/* L390: */
    }
/*  test whether the solution of the least-squares problem with the */
/*  constraints ibind(1),...ibind(nbind) in equality form, satisfies */
/*  all of the constraints (2). */
/*<       k = 1 >*/
    k = 1;
/*  number counts the number of violated inequality constraints. */
/*<       number = 0 >*/
    number = 0;
/*<       do 440 j=1,n6 >*/
    i__3 = n6;
    for (j = 1; j <= i__3; ++j) {
/*<         l = ibind(k) >*/
	l = ibind[k];
/*<         k = k+1 >*/
	++k;
/*<         if(j.eq.l) go to 440 >*/
	if (j == l) {
	    goto L440;
	}
/*<         k = k-1 >*/
	--k;
/*  test whether constraint j is satisfied */
/*<         f = e(j)*(c(j)-c(j+1))+const(j)*(c(j+2)-c(j+1)) >*/
	f = e[j] * (c__[j] - c__[j + 1]) + const__[j] * (c__[j + 2] - c__[j + 
		1]);
/*<         if(f.le.0.) go to 440 >*/
	if (f <= 0.) {
	    goto L440;
	}
/*  if constraint j is not satisfied, add a branch of length nbind+1 */
/*  to the tree. the nodes of this branch contain in their information */
/*  field the number of the constraints ibind(1),...ibind(nbind) and j, */
/*  arranged in increasing order. */
/*<         number = number+1 >*/
	++number;
/*<         k1 = k-1 >*/
	k1 = k - 1;
/*<         if(k1.eq.0) go to 410 >*/
	if (k1 == 0) {
	    goto L410;
	}
/*<         do 400 i=1,k1 >*/
	i__1 = k1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<           jbind(i) = ibind(i) >*/
	    jbind[i__] = ibind[i__];
/*<  400    continue >*/
/* L400: */
	}
/*<  410    jbind(k) = j >*/
L410:
	jbind[k] = j;
/*<         if(l.eq.0) go to 430 >*/
	if (l == 0) {
	    goto L430;
	}
/*<         do 420 i=k,nbind >*/
	i__1 = nbind;
	for (i__ = k; i__ <= i__1; ++i__) {
/*<           jbind(i+1) = ibind(i) >*/
	    jbind[i__ + 1] = ibind[i__];
/*<  420    continue >*/
/* L420: */
	}
/*<  430    call fpadno(maxtr,up,left,right,info,count,merk,jbind,n1,ier) >*/
L430:
	fpadno_(maxtr, &up[1], &left[1], &right[1], &info[1], &count, &merk, &
		jbind[1], &n1, ier);
/*  test whether the storage space which is required for the tree,exceeds */
/*  the available storage space. */
/*<         if(ier.ne.0) go to 560 >*/
	if (*ier != 0) {
	    goto L560;
	}
/*<  440  continue >*/
L440:
	;
    }
/*  test whether the solution of the least-squares problem with equality */
/*  constraints is a feasible solution. */
/*<       if(number.eq.0) go to 470 >*/
    if (number == 0) {
	goto L470;
    }
/*  test whether there are still cases with nbind constraints in */
/*  equality form to be considered. */
/*<  450  if(merk.gt.1) go to 460 >*/
L450:
    if (merk > 1) {
	goto L460;
    }
/*<       nbind = n1 >*/
    nbind = n1;
/*  test whether the number of knots where s''(x)=0 exceeds maxbin. */
/*<       if(nbind.gt.maxbin) go to 550 >*/
    if (nbind > *maxbin) {
	goto L550;
    }
/*<       n1 = n1+1 >*/
    ++n1;
/*<       ibind(n1) = 0 >*/
    ibind[n1] = 0;
/*  search which cases with nbind constraints in equality form */
/*  are going to be considered. */
/*<       call fpdeno(maxtr,up,left,right,nbind,merk) >*/
    fpdeno_(maxtr, &up[1], &left[1], &right[1], &nbind, &merk);
/*  test whether the quadratic programming problem has a solution. */
/*<       if(merk.eq.1) go to 570 >*/
    if (merk == 1) {
	goto L570;
    }
/*  find a new case with nbind constraints in equality form. */
/*<  460  call fpseno(maxtr,up,left,right,info,merk,ibind,nbind) >*/
L460:
    fpseno_(maxtr, &up[1], &left[1], &right[1], &info[1], &merk, &ibind[1], &
	    nbind);
/*<       go to 150 >*/
    goto L150;
/*  test whether the feasible solution is optimal. */
/*<  470  ier = 0 >*/
L470:
    *ier = 0;
/*<       do 480 i=1,n6 >*/
    i__3 = n6;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         bind(i) = .false. >*/
	bind[i__] = FALSE_;
/*<  480  continue >*/
/* L480: */
    }
/*<       if(nbind.eq.0) go to 500 >*/
    if (nbind == 0) {
	goto L500;
    }
/*<       do 490 i=1,nbind >*/
    i__3 = nbind;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<         if(u(i).le.0.) go to 450 >*/
	if (u[i__] <= 0.) {
	    goto L450;
	}
/*<         j = ibind(i) >*/
	j = ibind[i__];
/*<         bind(j) = .true. >*/
	bind[j] = TRUE_;
/*<  490  continue >*/
/* L490: */
    }
/*  evaluate s(x) at the data points x(i) and calculate the weighted */
/*  sum of squared residual right hand sides sq. */
/*<  500  sq = 0. >*/
L500:
    *sq = 0.;
/*<       l = 4 >*/
    l = 4;
/*<       lp1 = 5 >*/
    lp1 = 5;
/*<       do 530 i=1,m >*/
    i__3 = *m;
    for (i__ = 1; i__ <= i__3; ++i__) {
/*<  510    if(x(i).lt.t(lp1) .or. l.eq.n4) go to 520 >*/
L510:
	if (x[i__] < t[lp1] || l == n4) {
	    goto L520;
	}
/*<         l = lp1 >*/
	l = lp1;
/*<         lp1 = l+1 >*/
	lp1 = l + 1;
/*<         go to 510 >*/
	goto L510;
/*<  520    sx(i) = c(l-3)*q(i,1)+c(l-2)*q(i,2)+c(l-1)*q(i,3)+c(l)*q(i,4) >*/
L520:
	sx[i__] = c__[l - 3] * q[i__ + q_dim1] + c__[l - 2] * q[i__ + (q_dim1 
		<< 1)] + c__[l - 1] * q[i__ + q_dim1 * 3] + c__[l] * q[i__ + (
		q_dim1 << 2)];
/*<         sq = sq+(w(i)*(y(i)-sx(i)))**2 >*/
/* Computing 2nd power */
	d__1 = w[i__] * (y[i__] - sx[i__]);
	*sq += d__1 * d__1;
/*<  530  continue >*/
/* L530: */
    }
/*<       go to 600 >*/
    goto L600;
/*  error codes and messages. */
/*<  550  ier = 1 >*/
L550:
    *ier = 1;
/*<       go to 600 >*/
    goto L600;
/*<  560  ier = 2 >*/
L560:
    *ier = 2;
/*<       go to 600 >*/
    goto L600;
/*<  570  ier = 3 >*/
L570:
    *ier = 3;
/*<  600  return >*/
L600:
    return 0;
/*<       end >*/
} /* fpcosp_ */

