/* fptrpe.f -- translated by f2c (version 20061008).
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

/*<       subroutine fptrpe(m,mm,idim,n,nr,sp,p,b,z,a,aa,q,right) >*/
/* Subroutine */ int fptrpe_(integer *m, integer *mm, integer *idim, integer *
	n, integer *nr, doublereal *sp, doublereal *p, doublereal *b, 
	doublereal *z__, doublereal *a, doublereal *aa, doublereal *q, 
	doublereal *right)
{
    /* System generated locals */
    integer sp_dim1, sp_offset, b_dim1, b_offset, a_dim1, a_offset, aa_dim1, 
	    aa_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal h__[5];
    static integer i__, j, l;
    static doublereal h1[5], h2[4];
    static integer i2, i3, l0, m1, m2, m3, n1, n4, l1, i1, n7, j1, n11;
    static doublereal co;
    static integer ii, jj, jk, ik, ij;
    static doublereal si;
    static integer it, mid, nmd;
    static doublereal one, piv;
    static integer jper;
    static doublereal pinv;
    static integer irot, nrold, number;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*  subroutine fptrpe reduces the (m+n-7) x (n-7) cyclic bandmatrix a */
/*  to upper triangular form and applies the same givens transformations */
/*  to the (m) x (mm) x (idim) matrix z to obtain the (n-7) x (mm) x */
/*  (idim) matrix q. */
/*  .. */
/*  ..scalar arguments.. */
/*<       real p >*/
/*<       integer m,mm,idim,n >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer nr(m) >*/
/*  ..local scalars.. */
/*<       real co,pinv,piv,si,one >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real h(5),h1(5),h2(4) >*/
/*  ..subroutine references.. */
/*    fpgivs,fprota */
/*  .. */
/*<       one = 1 >*/
    /* Parameter adjustments */
    sp_dim1 = *m;
    sp_offset = 1 + sp_dim1;
    sp -= sp_offset;
    --nr;
    --right;
    --z__;
    --q;
    aa_dim1 = *n;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    one = 1.;
/*<       if(p.gt.0.) pinv = one/p >*/
    if (*p > 0.) {
	pinv = one / *p;
    }
/*<       n4 = n-4 >*/
    n4 = *n - 4;
/*<       n7 = n-7 >*/
    n7 = *n - 7;
/*<       n11 = n-11 >*/
    n11 = *n - 11;
/*<       mid = mm*idim >*/
    mid = *mm * *idim;
/*<       m2 = m*mm >*/
    m2 = *m * *mm;
/*<       m3 = n7*mm >*/
    m3 = n7 * *mm;
/*<       m1 = m-1 >*/
    m1 = *m - 1;
/*  we determine the matrix (a) and then we reduce her to */
/*  upper triangular form (r) using givens rotations. */
/*  we apply the same transformations to the rows of matrix */
/*  z to obtain the (mm) x (n-7) matrix g. */
/*  we store matrix (r) into a and aa, g into q. */
/*  the n7 x n7 upper triangular matrix (r) has the form */
/*             | a1 '     | */
/*       (r) = |    ' a2  | */
/*             |  0 '     | */
/*  with (a2) a n7 x 4 matrix and (a1) a n11 x n11 upper */
/*  triangular matrix of bandwidth 5. */
/*  initialization. */
/*<       nmd = n7*mid >*/
    nmd = n7 * mid;
/*<       do 50 i=1,nmd >*/
    i__1 = nmd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         q(i) = 0. >*/
	q[i__] = 0.;
/*<   50  continue >*/
/* L50: */
    }
/*<       do 100 i=1,n4 >*/
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         a(i,5) = 0. >*/
	a[i__ + a_dim1 * 5] = 0.;
/*<         do 100 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           a(i,j) = 0. >*/
	    a[i__ + j * a_dim1] = 0.;
/*<           aa(i,j) = 0. >*/
	    aa[i__ + j * aa_dim1] = 0.;
/*<  100  continue >*/
/* L100: */
	}
    }
/*<       jper = 0 >*/
    jper = 0;
/*<       nrold = 0 >*/
    nrold = 0;
/*<       do 760 it=1,m1 >*/
    i__1 = m1;
    for (it = 1; it <= i__1; ++it) {
/*<         number = nr(it) >*/
	number = nr[it];
/*<  120    if(nrold.eq.number) go to 180 >*/
L120:
	if (nrold == number) {
	    goto L180;
	}
/*<         if(p.le.0.) go to 740 >*/
	if (*p <= 0.) {
	    goto L740;
	}
/*  fetch a new row of matrix (b). */
/*<         n1 = nrold+1 >*/
	n1 = nrold + 1;
/*<         do 140 j=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<           h(j) = b(n1,j)*pinv >*/
	    h__[j - 1] = b[n1 + j * b_dim1] * pinv;
/*<  140    continue >*/
/* L140: */
	}
/*  find the appropiate row of q. */
/*<         do 160 j=1,mid >*/
	i__2 = mid;
	for (j = 1; j <= i__2; ++j) {
/*<           right(j) = 0. >*/
	    right[j] = 0.;
/*<  160    continue >*/
/* L160: */
	}
/*<         go to 240 >*/
	goto L240;
/*  fetch a new row of matrix (sp) */
/*<  180    h(5) = 0. >*/
L180:
	h__[4] = 0.;
/*<         do 200 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           h(j) = sp(it,j) >*/
	    h__[j - 1] = sp[it + j * sp_dim1];
/*<  200    continue >*/
/* L200: */
	}
/*  find the appropiate row of q. */
/*<         j = 0 >*/
	j = 0;
/*<         do 220 ii=1,idim >*/
	i__2 = *idim;
	for (ii = 1; ii <= i__2; ++ii) {
/*<           l = (ii-1)*m2+(it-1)*mm >*/
	    l = (ii - 1) * m2 + (it - 1) * *mm;
/*<           do 220 jj=1,mm >*/
	    i__3 = *mm;
	    for (jj = 1; jj <= i__3; ++jj) {
/*<             j = j+1 >*/
		++j;
/*<             l = l+1 >*/
		++l;
/*<             right(j) = z(l) >*/
		right[j] = z__[l];
/*<  220    continue >*/
/* L220: */
	    }
	}
/*  test whether there are non-zero values in the new row of (a) */
/*  corresponding to the b-splines n(j,*),j=n7+1,...,n4. */
/*<  240     if(nrold.lt.n11) go to 640 >*/
L240:
	if (nrold < n11) {
	    goto L640;
	}
/*<          if(jper.ne.0) go to 320 >*/
	if (jper != 0) {
	    goto L320;
	}
/*  initialize the matrix (aa). */
/*<          jk = n11+1 >*/
	jk = n11 + 1;
/*<          do 300 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<             ik = jk >*/
	    ik = jk;
/*<             do 260 j=1,5 >*/
	    for (j = 1; j <= 5; ++j) {
/*<                if(ik.le.0) go to 280 >*/
		if (ik <= 0) {
		    goto L280;
		}
/*<                aa(ik,i) = a(ik,j) >*/
		aa[ik + i__ * aa_dim1] = a[ik + j * a_dim1];
/*<                ik = ik-1 >*/
		--ik;
/*<  260        continue >*/
/* L260: */
	    }
/*<  280        jk = jk+1 >*/
L280:
	    ++jk;
/*<  300     continue >*/
/* L300: */
	}
/*<          jper = 1 >*/
	jper = 1;
/*  if one of the non-zero elements of the new row corresponds to one of */
/*  the b-splines n(j;*),j=n7+1,...,n4,we take account of the periodicity */
/*  conditions for setting up this row of (a). */
/*<  320     do 340 i=1,4 >*/
L320:
	for (i__ = 1; i__ <= 4; ++i__) {
/*<             h1(i) = 0. >*/
	    h1[i__ - 1] = 0.;
/*<             h2(i) = 0. >*/
	    h2[i__ - 1] = 0.;
/*<  340     continue >*/
/* L340: */
	}
/*<          h1(5) = 0. >*/
	h1[4] = 0.;
/*<          j = nrold-n11 >*/
	j = nrold - n11;
/*<          do 420 i=1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<             j = j+1 >*/
	    ++j;
/*<             l0 = j >*/
	    l0 = j;
/*<  360        l1 = l0-4 >*/
L360:
	    l1 = l0 - 4;
/*<             if(l1.le.0) go to 400 >*/
	    if (l1 <= 0) {
		goto L400;
	    }
/*<             if(l1.le.n11) go to 380 >*/
	    if (l1 <= n11) {
		goto L380;
	    }
/*<             l0 = l1-n11 >*/
	    l0 = l1 - n11;
/*<             go to 360 >*/
	    goto L360;
/*<  380        h1(l1) = h(i) >*/
L380:
	    h1[l1 - 1] = h__[i__ - 1];
/*<             go to 420 >*/
	    goto L420;
/*<  400        h2(l0) = h2(l0) + h(i) >*/
L400:
	    h2[l0 - 1] += h__[i__ - 1];
/*<  420     continue >*/
L420:
	    ;
	}
/*  rotate the new row of (a) into triangle. */
/*<          if(n11.le.0) go to 560 >*/
	if (n11 <= 0) {
	    goto L560;
	}
/*  rotations with the rows 1,2,...,n11 of (a). */
/*<          do 540 irot=1,n11 >*/
	i__3 = n11;
	for (irot = 1; irot <= i__3; ++irot) {
/*<             piv = h1(1) >*/
	    piv = h1[0];
/*<             i2 = min0(n11-irot,4) >*/
/* Computing MIN */
	    i__2 = n11 - irot;
	    i2 = min(i__2,4);
/*<             if(piv.eq.0.) go to 500 >*/
	    if (piv == 0.) {
		goto L500;
	    }
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a(irot,1),co,si) >*/
	    fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the columns of matrix q. */
/*<             j = 0 >*/
	    j = 0;
/*<             do 440 ii=1,idim >*/
	    i__2 = *idim;
	    for (ii = 1; ii <= i__2; ++ii) {
/*<                l = (ii-1)*m3+irot >*/
		l = (ii - 1) * m3 + irot;
/*<                do 440 jj=1,mm >*/
		i__4 = *mm;
		for (jj = 1; jj <= i__4; ++jj) {
/*<                  j = j+1 >*/
		    ++j;
/*<                  call fprota(co,si,right(j),q(l)) >*/
		    fprota_(&co, &si, &right[j], &q[l]);
/*<                  l = l+n7 >*/
		    l += n7;
/*<  440        continue >*/
/* L440: */
		}
	    }
/*  apply that transformation to the rows of (a) with respect to aa. */
/*<             do 460 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<                call fprota(co,si,h2(i),aa(irot,i)) >*/
		fprota_(&co, &si, &h2[i__ - 1], &aa[irot + i__ * aa_dim1]);
/*<  460        continue >*/
/* L460: */
	    }
/*  apply that transformation to the rows of (a) with respect to a. */
/*<             if(i2.eq.0) go to 560 >*/
	    if (i2 == 0) {
		goto L560;
	    }
/*<             do 480 i=1,i2 >*/
	    i__4 = i2;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<                i1 = i+1 >*/
		i1 = i__ + 1;
/*<                call fprota(co,si,h1(i1),a(irot,i1)) >*/
		fprota_(&co, &si, &h1[i1 - 1], &a[irot + i1 * a_dim1]);
/*<  480        continue >*/
/* L480: */
	    }
/*<  500        do 520 i=1,i2 >*/
L500:
	    i__4 = i2;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<                h1(i) = h1(i+1) >*/
		h1[i__ - 1] = h1[i__];
/*<  520        continue >*/
/* L520: */
	    }
/*<             h1(i2+1) = 0. >*/
	    h1[i2] = 0.;
/*<  540     continue >*/
/* L540: */
	}
/*  rotations with the rows n11+1,...,n7 of a. */
/*<  560     do 620 irot=1,4 >*/
L560:
	for (irot = 1; irot <= 4; ++irot) {
/*<             ij = n11+irot >*/
	    ij = n11 + irot;
/*<             if(ij.le.0) go to 620 >*/
	    if (ij <= 0) {
		goto L620;
	    }
/*<             piv = h2(irot) >*/
	    piv = h2[irot - 1];
/*<             if(piv.eq.0.) go to 620 >*/
	    if (piv == 0.) {
		goto L620;
	    }
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,aa(ij,irot),co,si) >*/
	    fpgivs_(&piv, &aa[ij + irot * aa_dim1], &co, &si);
/*  apply that transformation to the columns of matrix q. */
/*<             j = 0 >*/
	    j = 0;
/*<             do 580 ii=1,idim >*/
	    i__3 = *idim;
	    for (ii = 1; ii <= i__3; ++ii) {
/*<                l = (ii-1)*m3+ij >*/
		l = (ii - 1) * m3 + ij;
/*<                do 580 jj=1,mm >*/
		i__4 = *mm;
		for (jj = 1; jj <= i__4; ++jj) {
/*<                  j = j+1 >*/
		    ++j;
/*<                  call fprota(co,si,right(j),q(l)) >*/
		    fprota_(&co, &si, &right[j], &q[l]);
/*<                  l = l+n7 >*/
		    l += n7;
/*<  580        continue >*/
/* L580: */
		}
	    }
/*<             if(irot.eq.4) go to 620 >*/
	    if (irot == 4) {
		goto L620;
	    }
/*  apply that transformation to the rows of (a) with respect to aa. */
/*<             j1 = irot+1 >*/
	    j1 = irot + 1;
/*<             do 600 i=j1,4 >*/
	    for (i__ = j1; i__ <= 4; ++i__) {
/*<                call fprota(co,si,h2(i),aa(ij,i)) >*/
		fprota_(&co, &si, &h2[i__ - 1], &aa[ij + i__ * aa_dim1]);
/*<  600        continue >*/
/* L600: */
	    }
/*<  620     continue >*/
L620:
	    ;
	}
/*<          go to 720 >*/
	goto L720;
/*  rotation into triangle of the new row of (a), in case the elements */
/*  corresponding to the b-splines n(j;*),j=n7+1,...,n4 are all zero. */
/*<  640     irot =nrold >*/
L640:
	irot = nrold;
/*<          do 700 i=1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<             irot = irot+1 >*/
	    ++irot;
/*<             piv = h(i) >*/
	    piv = h__[i__ - 1];
/*<             if(piv.eq.0.) go to 700 >*/
	    if (piv == 0.) {
		goto L700;
	    }
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a(irot,1),co,si) >*/
	    fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
/*<             j = 0 >*/
	    j = 0;
/*<             do 660 ii=1,idim >*/
	    i__4 = *idim;
	    for (ii = 1; ii <= i__4; ++ii) {
/*<                l = (ii-1)*m3+irot >*/
		l = (ii - 1) * m3 + irot;
/*<                do 660 jj=1,mm >*/
		i__3 = *mm;
		for (jj = 1; jj <= i__3; ++jj) {
/*<                  j = j+1 >*/
		    ++j;
/*<                  call fprota(co,si,right(j),q(l)) >*/
		    fprota_(&co, &si, &right[j], &q[l]);
/*<                  l = l+n7 >*/
		    l += n7;
/*<  660        continue >*/
/* L660: */
		}
	    }
/*  apply that transformation to the rows of (a). */
/*<             if(i.eq.5) go to 700 >*/
	    if (i__ == 5) {
		goto L700;
	    }
/*<             i2 = 1 >*/
	    i2 = 1;
/*<             i3 = i+1 >*/
	    i3 = i__ + 1;
/*<             do 680 j=i3,5 >*/
	    for (j = i3; j <= 5; ++j) {
/*<                i2 = i2+1 >*/
		++i2;
/*<                call fprota(co,si,h(j),a(irot,i2)) >*/
		fprota_(&co, &si, &h__[j - 1], &a[irot + i2 * a_dim1]);
/*<  680        continue >*/
/* L680: */
	    }
/*<  700     continue >*/
L700:
	    ;
	}
/*<  720     if(nrold.eq.number) go to 760 >*/
L720:
	if (nrold == number) {
	    goto L760;
	}
/*<  740     nrold = nrold+1 >*/
L740:
	++nrold;
/*<          go to 120 >*/
	goto L120;
/*<  760  continue >*/
L760:
	;
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fptrpe_ */

