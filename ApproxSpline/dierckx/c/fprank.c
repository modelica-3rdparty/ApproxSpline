/* fprank.f -- translated by f2c (version 20061008).
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

/*<       subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h) >*/
/* Subroutine */ int fprank_(doublereal *a, doublereal *f, integer *n, 
	integer *m, integer *na, doublereal *tol, doublereal *c__, doublereal 
	*sq, integer *rank, doublereal *aa, doublereal *ff, doublereal *h__)
{
    /* System generated locals */
    integer a_dim1, a_offset, aa_dim1, aa_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i1, i2, j1, j2, j3, m1, ii, ij, jj, kk, nl;
    static doublereal yi, fac, cos__, sin__, piv, stor1, stor2, stor3, store;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*  subroutine fprank finds the minimum norm solution of a least- */
/*  squares problem in case of rank deficiency. */

/*  input parameters: */
/*    a : array, which contains the non-zero elements of the observation */
/*        matrix after triangularization by givens transformations. */
/*    f : array, which contains the transformed right hand side. */
/*    n : integer,wich contains the dimension of a. */
/*    m : integer, which denotes the bandwidth of a. */
/*  tol : real value, giving a threshold to determine the rank of a. */

/*  output parameters: */
/*    c : array, which contains the minimum norm solution. */
/*   sq : real value, giving the contribution of reducing the rank */
/*        to the sum of squared residuals. */
/* rank : integer, which contains the rank of matrix a. */

/*  ..scalar arguments.. */
/*<       integer n,m,na,rank >*/
/*<       real tol,sq >*/
/*  ..array arguments.. */
/*<       real a(na,m),f(n),c(n),aa(n,m),ff(n),h(m) >*/
/*  ..local scalars.. */
/*<       integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl >*/
/*<       real cos,fac,piv,sin,yi >*/
/*<       double precision store,stor1,stor2,stor3 >*/
/*  ..function references.. */
/*<       integer min0 >*/
/*  ..subroutine references.. */
/*    fpgivs,fprota */
/*  .. */
/*<       m1 = m-1 >*/
    /* Parameter adjustments */
    --ff;
    --c__;
    --f;
    --h__;
    aa_dim1 = *n;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    a_dim1 = *na;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    m1 = *m - 1;
/*  the rank deficiency nl is considered to be the number of sufficient */
/*  small diagonal elements of a. */
/*<       nl = 0 >*/
    nl = 0;
/*<       sq = 0. >*/
    *sq = 0.;
/*<       do 90 i=1,n >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         if(a(i,1).gt.tol) go to 90 >*/
	if (a[i__ + a_dim1] > *tol) {
	    goto L90;
	}
/*  if a sufficient small diagonal element is found, we put it to */
/*  zero. the remainder of the row corresponding to that zero diagonal */
/*  element is then rotated into triangle by givens rotations . */
/*  the rank deficiency is increased by one. */
/*<         nl = nl+1 >*/
	++nl;
/*<         if(i.eq.n) go to 90 >*/
	if (i__ == *n) {
	    goto L90;
	}
/*<         yi = f(i) >*/
	yi = f[i__];
/*<         do 10 j=1,m1 >*/
	i__2 = m1;
	for (j = 1; j <= i__2; ++j) {
/*<           h(j) = a(i,j+1) >*/
	    h__[j] = a[i__ + (j + 1) * a_dim1];
/*<   10    continue >*/
/* L10: */
	}
/*<         h(m) = 0. >*/
	h__[*m] = 0.;
/*<         i1 = i+1 >*/
	i1 = i__ + 1;
/*<         do 60 ii=i1,n >*/
	i__2 = *n;
	for (ii = i1; ii <= i__2; ++ii) {
/*<           i2 = min0(n-ii,m1) >*/
/* Computing MIN */
	    i__3 = *n - ii;
	    i2 = min(i__3,m1);
/*<           piv = h(1) >*/
	    piv = h__[1];
/*<           if(piv.eq.0.) go to 30 >*/
	    if (piv == 0.) {
		goto L30;
	    }
/*<           call fpgivs(piv,a(ii,1),cos,sin) >*/
	    fpgivs_(&piv, &a[ii + a_dim1], &cos__, &sin__);
/*<           call fprota(cos,sin,yi,f(ii)) >*/
	    fprota_(&cos__, &sin__, &yi, &f[ii]);
/*<           if(i2.eq.0) go to 70 >*/
	    if (i2 == 0) {
		goto L70;
	    }
/*<           do 20 j=1,i2 >*/
	    i__3 = i2;
	    for (j = 1; j <= i__3; ++j) {
/*<             j1 = j+1 >*/
		j1 = j + 1;
/*<             call fprota(cos,sin,h(j1),a(ii,j1)) >*/
		fprota_(&cos__, &sin__, &h__[j1], &a[ii + j1 * a_dim1]);
/*<             h(j) = h(j1) >*/
		h__[j] = h__[j1];
/*<   20      continue >*/
/* L20: */
	    }
/*<           go to 50 >*/
	    goto L50;
/*<   30      if(i2.eq.0) go to 70 >*/
L30:
	    if (i2 == 0) {
		goto L70;
	    }
/*<           do 40 j=1,i2 >*/
	    i__3 = i2;
	    for (j = 1; j <= i__3; ++j) {
/*<             h(j) = h(j+1) >*/
		h__[j] = h__[j + 1];
/*<   40      continue >*/
/* L40: */
	    }
/*<   50      h(i2+1) = 0. >*/
L50:
	    h__[i2 + 1] = 0.;
/*<   60    continue >*/
/* L60: */
	}
/*  add to the sum of squared residuals the contribution of deleting */
/*  the row with small diagonal element. */
/*<   70    sq = sq+yi**2 >*/
L70:
/* Computing 2nd power */
	d__1 = yi;
	*sq += d__1 * d__1;
/*<   90  continue >*/
L90:
	;
    }
/*  rank denotes the rank of a. */
/*<       rank = n-nl >*/
    *rank = *n - nl;
/*  let b denote the (rank*n) upper trapezoidal matrix which can be */
/*  obtained from the (n*n) upper triangular matrix a by deleting */
/*  the rows and interchanging the columns corresponding to a zero */
/*  diagonal element. if this matrix is factorized using givens */
/*  transformations as  b = (r) (u)  where */
/*    r is a (rank*rank) upper triangular matrix, */
/*    u is a (rank*n) orthonormal matrix */
/*  then the minimal least-squares solution c is given by c = b' v, */
/*  where v is the solution of the system  (r) (r)' v = g  and */
/*  g denotes the vector obtained from the old right hand side f, by */
/*  removing the elements corresponding to a zero diagonal element of a. */
/*  initialization. */
/*<       do 100 i=1,rank >*/
    i__1 = *rank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         do 100 j=1,m >*/
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
/*<           aa(i,j) = 0. >*/
	    aa[i__ + j * aa_dim1] = 0.;
/*<  100  continue >*/
/* L100: */
	}
    }
/*  form in aa the upper triangular matrix obtained from a by */
/*  removing rows and columns with zero diagonal elements. form in ff */
/*  the new right hand side by removing the elements of the old right */
/*  hand side corresponding to a deleted row. */
/*<       ii = 0 >*/
    ii = 0;
/*<       do 120 i=1,n >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         if(a(i,1).le.tol) go to 120 >*/
	if (a[i__ + a_dim1] <= *tol) {
	    goto L120;
	}
/*<         ii = ii+1 >*/
	++ii;
/*<         ff(ii) = f(i) >*/
	ff[ii] = f[i__];
/*<         aa(ii,1) = a(i,1) >*/
	aa[ii + aa_dim1] = a[i__ + a_dim1];
/*<         jj = ii >*/
	jj = ii;
/*<         kk = 1 >*/
	kk = 1;
/*<         j = i >*/
	j = i__;
/*<         j1 = min0(j-1,m1) >*/
/* Computing MIN */
	i__1 = j - 1;
	j1 = min(i__1,m1);
/*<         if(j1.eq.0) go to 120 >*/
	if (j1 == 0) {
	    goto L120;
	}
/*<         do 110 k=1,j1 >*/
	i__1 = j1;
	for (k = 1; k <= i__1; ++k) {
/*<           j = j-1 >*/
	    --j;
/*<           if(a(j,1).le.tol) go to 110 >*/
	    if (a[j + a_dim1] <= *tol) {
		goto L110;
	    }
/*<           kk = kk+1 >*/
	    ++kk;
/*<           jj = jj-1 >*/
	    --jj;
/*<           aa(jj,kk) = a(j,k+1) >*/
	    aa[jj + kk * aa_dim1] = a[j + (k + 1) * a_dim1];
/*<  110    continue >*/
L110:
	    ;
	}
/*<  120  continue >*/
L120:
	;
    }
/*  form successively in h the columns of a with a zero diagonal element. */
/*<       ii = 0 >*/
    ii = 0;
/*<       do 200 i=1,n >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         ii = ii+1 >*/
	++ii;
/*<         if(a(i,1).gt.tol) go to 200 >*/
	if (a[i__ + a_dim1] > *tol) {
	    goto L200;
	}
/*<         ii = ii-1 >*/
	--ii;
/*<         if(ii.eq.0) go to 200 >*/
	if (ii == 0) {
	    goto L200;
	}
/*<         jj = 1 >*/
	jj = 1;
/*<         j = i >*/
	j = i__;
/*<         j1 = min0(j-1,m1) >*/
/* Computing MIN */
	i__1 = j - 1;
	j1 = min(i__1,m1);
/*<         do 130 k=1,j1 >*/
	i__1 = j1;
	for (k = 1; k <= i__1; ++k) {
/*<           j = j-1 >*/
	    --j;
/*<           if(a(j,1).le.tol) go to 130 >*/
	    if (a[j + a_dim1] <= *tol) {
		goto L130;
	    }
/*<           h(jj) = a(j,k+1) >*/
	    h__[jj] = a[j + (k + 1) * a_dim1];
/*<           jj = jj+1 >*/
	    ++jj;
/*<  130    continue >*/
L130:
	    ;
	}
/*<         do 140 kk=jj,m >*/
	i__1 = *m;
	for (kk = jj; kk <= i__1; ++kk) {
/*<           h(kk) = 0. >*/
	    h__[kk] = 0.;
/*<  140    continue >*/
/* L140: */
	}
/*  rotate this column into aa by givens transformations. */
/*<         jj = ii >*/
	jj = ii;
/*<         do 190 i1=1,ii >*/
	i__1 = ii;
	for (i1 = 1; i1 <= i__1; ++i1) {
/*<           j1 = min0(jj-1,m1) >*/
/* Computing MIN */
	    i__3 = jj - 1;
	    j1 = min(i__3,m1);
/*<           piv = h(1) >*/
	    piv = h__[1];
/*<           if(piv.ne.0.) go to 160 >*/
	    if (piv != 0.) {
		goto L160;
	    }
/*<           if(j1.eq.0) go to 200 >*/
	    if (j1 == 0) {
		goto L200;
	    }
/*<           do 150 j2=1,j1 >*/
	    i__3 = j1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
/*<             j3 = j2+1 >*/
		j3 = j2 + 1;
/*<             h(j2) = h(j3) >*/
		h__[j2] = h__[j3];
/*<  150      continue >*/
/* L150: */
	    }
/*<           go to 180 >*/
	    goto L180;
/*<  160      call fpgivs(piv,aa(jj,1),cos,sin) >*/
L160:
	    fpgivs_(&piv, &aa[jj + aa_dim1], &cos__, &sin__);
/*<           if(j1.eq.0) go to 200 >*/
	    if (j1 == 0) {
		goto L200;
	    }
/*<           kk = jj >*/
	    kk = jj;
/*<           do 170 j2=1,j1 >*/
	    i__3 = j1;
	    for (j2 = 1; j2 <= i__3; ++j2) {
/*<             j3 = j2+1 >*/
		j3 = j2 + 1;
/*<             kk = kk-1 >*/
		--kk;
/*<             call fprota(cos,sin,h(j3),aa(kk,j3)) >*/
		fprota_(&cos__, &sin__, &h__[j3], &aa[kk + j3 * aa_dim1]);
/*<             h(j2) = h(j3) >*/
		h__[j2] = h__[j3];
/*<  170      continue >*/
/* L170: */
	    }
/*<  180      jj = jj-1 >*/
L180:
	    --jj;
/*<           h(j3) = 0. >*/
	    h__[j3] = 0.;
/*<  190    continue >*/
/* L190: */
	}
/*<  200  continue >*/
L200:
	;
    }
/*  solve the system (aa) (f1) = ff */
/*<       ff(rank) = ff(rank)/aa(rank,1) >*/
    ff[*rank] /= aa[*rank + aa_dim1];
/*<       i = rank-1 >*/
    i__ = *rank - 1;
/*<       if(i.eq.0) go to 230 >*/
    if (i__ == 0) {
	goto L230;
    }
/*<       do 220 j=2,rank >*/
    i__2 = *rank;
    for (j = 2; j <= i__2; ++j) {
/*<         store = ff(i) >*/
	store = ff[i__];
/*<         i1 = min0(j-1,m1) >*/
/* Computing MIN */
	i__1 = j - 1;
	i1 = min(i__1,m1);
/*<         k = i >*/
	k = i__;
/*<         do 210 ii=1,i1 >*/
	i__1 = i1;
	for (ii = 1; ii <= i__1; ++ii) {
/*<           k = k+1 >*/
	    ++k;
/*<           stor1 = ff(k) >*/
	    stor1 = ff[k];
/*<           stor2 = aa(i,ii+1) >*/
	    stor2 = aa[i__ + (ii + 1) * aa_dim1];
/*<           store = store-stor1*stor2 >*/
	    store -= stor1 * stor2;
/*<  210    continue >*/
/* L210: */
	}
/*<         stor1 = aa(i,1) >*/
	stor1 = aa[i__ + aa_dim1];
/*<         ff(i) = store/stor1 >*/
	ff[i__] = store / stor1;
/*<         i = i-1 >*/
	--i__;
/*<  220  continue >*/
/* L220: */
    }
/*  solve the system  (aa)' (f2) = f1 */
/*<  230  ff(1) = ff(1)/aa(1,1) >*/
L230:
    ff[1] /= aa[aa_dim1 + 1];
/*<       if(rank.eq.1) go to 260 >*/
    if (*rank == 1) {
	goto L260;
    }
/*<       do 250 j=2,rank >*/
    i__2 = *rank;
    for (j = 2; j <= i__2; ++j) {
/*<         store = ff(j) >*/
	store = ff[j];
/*<         i1 = min0(j-1,m1) >*/
/* Computing MIN */
	i__1 = j - 1;
	i1 = min(i__1,m1);
/*<         k = j >*/
	k = j;
/*<         do 240 ii=1,i1 >*/
	i__1 = i1;
	for (ii = 1; ii <= i__1; ++ii) {
/*<           k = k-1 >*/
	    --k;
/*<           stor1 = ff(k) >*/
	    stor1 = ff[k];
/*<           stor2 = aa(k,ii+1) >*/
	    stor2 = aa[k + (ii + 1) * aa_dim1];
/*<           store = store-stor1*stor2 >*/
	    store -= stor1 * stor2;
/*<  240    continue >*/
/* L240: */
	}
/*<         stor1 = aa(j,1) >*/
	stor1 = aa[j + aa_dim1];
/*<         ff(j) = store/stor1 >*/
	ff[j] = store / stor1;
/*<  250  continue >*/
/* L250: */
    }
/*  premultiply f2 by the transpoze of a. */
/*<  260  k = 0 >*/
L260:
    k = 0;
/*<       do 280 i=1,n >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         store = 0. >*/
	store = 0.;
/*<         if(a(i,1).gt.tol) k = k+1 >*/
	if (a[i__ + a_dim1] > *tol) {
	    ++k;
	}
/*<         j1 = min0(i,m) >*/
	j1 = min(i__,*m);
/*<         kk = k >*/
	kk = k;
/*<         ij = i+1 >*/
	ij = i__ + 1;
/*<         do 270 j=1,j1 >*/
	i__1 = j1;
	for (j = 1; j <= i__1; ++j) {
/*<           ij = ij-1 >*/
	    --ij;
/*<           if(a(ij,1).le.tol) go to 270 >*/
	    if (a[ij + a_dim1] <= *tol) {
		goto L270;
	    }
/*<           stor1 = a(ij,j) >*/
	    stor1 = a[ij + j * a_dim1];
/*<           stor2 = ff(kk) >*/
	    stor2 = ff[kk];
/*<           store = store+stor1*stor2 >*/
	    store += stor1 * stor2;
/*<           kk = kk-1 >*/
	    --kk;
/*<  270    continue >*/
L270:
	    ;
	}
/*<         c(i) = store >*/
	c__[i__] = store;
/*<  280  continue >*/
/* L280: */
    }
/*  add to the sum of squared residuals the contribution of putting */
/*  to zero the small diagonal elements of matrix (a). */
/*<       stor3 = 0. >*/
    stor3 = 0.;
/*<       do 310 i=1,n >*/
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         if(a(i,1).gt.tol) go to 310 >*/
	if (a[i__ + a_dim1] > *tol) {
	    goto L310;
	}
/*<         store = f(i) >*/
	store = f[i__];
/*<         i1 = min0(n-i,m1) >*/
/* Computing MIN */
	i__1 = *n - i__;
	i1 = min(i__1,m1);
/*<         if(i1.eq.0) go to 300 >*/
	if (i1 == 0) {
	    goto L300;
	}
/*<         do 290 j=1,i1 >*/
	i__1 = i1;
	for (j = 1; j <= i__1; ++j) {
/*<           ij = i+j >*/
	    ij = i__ + j;
/*<           stor1 = c(ij) >*/
	    stor1 = c__[ij];
/*<           stor2 = a(i,j+1) >*/
	    stor2 = a[i__ + (j + 1) * a_dim1];
/*<           store = store-stor1*stor2 >*/
	    store -= stor1 * stor2;
/*<  290    continue >*/
/* L290: */
	}
/*<  300    fac = a(i,1)*c(i) >*/
L300:
	fac = a[i__ + a_dim1] * c__[i__];
/*<         stor1 = a(i,1) >*/
	stor1 = a[i__ + a_dim1];
/*<         stor2 = c(i) >*/
	stor2 = c__[i__];
/*<         stor1 = stor1*stor2 >*/
	stor1 *= stor2;
/*<         stor3 = stor3+stor1*(stor1-store-store) >*/
	stor3 += stor1 * (stor1 - store - store);
/*<  310  continue >*/
L310:
	;
    }
/*<       fac = stor3 >*/
    fac = stor3;
/*<       sq = sq+fac >*/
    *sq += fac;
/*<       return >*/
    return 0;
/*<       end >*/
} /* fprank_ */

