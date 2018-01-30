/* fpgrsp.f -- translated by f2c (version 20061008).
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
static integer c__5 = 5;
static integer c__4 = 4;

/*<    >*/
/* Subroutine */ int fpgrsp_(integer *ifsu, integer *ifsv, integer *ifbu, 
	integer *ifbv, integer *iback, doublereal *u, integer *mu, doublereal 
	*v, integer *mv, doublereal *r__, integer *mr, doublereal *dr, 
	integer *iop0, integer *iop1, doublereal *tu, integer *nu, doublereal 
	*tv, integer *nv, doublereal *p, doublereal *c__, integer *nc, 
	doublereal *sq, doublereal *fp, doublereal *fpu, doublereal *fpv, 
	integer *mm, integer *mvnu, doublereal *spu, doublereal *spv, 
	doublereal *right, doublereal *q, doublereal *au, doublereal *av1, 
	doublereal *av2, doublereal *bu, doublereal *bv, doublereal *a0, 
	doublereal *a1, doublereal *b0, doublereal *b1, doublereal *c0, 
	doublereal *c1, doublereal *cosi, integer *nru, integer *nrv)
{
    /* System generated locals */
    integer spu_dim1, spu_offset, spv_dim1, spv_offset, au_dim1, au_offset, 
	    av1_dim1, av1_offset, av2_dim1, av2_offset, bu_dim1, bu_offset, 
	    bv_dim1, bv_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal h__[5];
    static integer i__, j, k, l, i0, i1, i2, i3, j0, j1, k1, k2, l0, l1, l2, 
	    n1;
    static doublereal h1[5], h2[4];
    static integer ic;
    static doublereal co;
    static integer ii, ij, ik;
    static doublereal si;
    static integer iq, it, ir, jj, jk, nu4, nv4, nu7, nu8, nu9, nv7, nv8;
    static doublereal fac, dr01, dr02, dr03, arg, dr11, dr12, dr13, one;
    static integer nv11;
    static doublereal piv, fac0, fac1;
    static integer mvv, nuu;
    static doublereal half;
    static integer ncof, jper;
    static doublereal term, pinv;
    static integer irot, numu, numv, numu1, numv1;
    static doublereal three;
    static integer nrold;
    extern /* Subroutine */ int fpcyt1_(doublereal *, integer *, integer *), 
	    fpcyt2_(doublereal *, integer *, doublereal *, doublereal *, 
	    integer *), fpback_(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *), fpbacp_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), fpdisc_(doublereal *, integer *, integer *,
	     doublereal *, integer *), fpbspl_(doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer number;
    extern /* Subroutine */ int fprota_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), fpgivs_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer nroldu, nroldv;

/*  .. */
/*  ..scalar arguments.. */
/*<       real p,sq,fp >*/
/*<    >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer nru(mu),nrv(mv) >*/
/*  ..local scalars.. */
/*<    >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real h(5),h1(5),h2(4) >*/
/*  ..function references.. */
/*<       integer min0 >*/
/*<       real cos,sin >*/
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota */
/*  .. */
/*  let */
/*               |     (spu)      |            |     (spv)      | */
/*        (au) = | -------------- |     (av) = | -------------- | */
/*               | sqrt(1/p) (bu) |            | sqrt(1/p) (bv) | */

/*                                | r  ' 0 | */
/*                            q = | ------ | */
/*                                | 0  ' 0 | */

/*  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline */
/*                coefficients. */
/*       r      : the mu x mv matrix which contains the function values. */
/*       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices */
/*                according to the least-squares problems in the u-,resp. */
/*                v-direction. */
/*       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices */
/*                containing the discontinuity jumps of the derivatives */
/*                of the b-splines in the u-,resp.v-variable at the knots */
/*  the b-spline coefficients of the smoothing spline are then calculated */
/*  as the least-squares solution of the following over-determined linear */
/*  system of equations */

/*  (1)  (av) c (au)' = q */

/*  subject to the constraints */

/*  (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4 */

/*  (3)  if iop0 = 0  c(1,j) = dr(1) */
/*          iop0 = 1  c(1,j) = dr(1) */
/*                    c(2,j) = dr(1)+(dr(2)*cosi(1,j)+dr(3)*cosi(2,j))* */
/*                            tu(5)/3. = c0(j) , j=1,2,...nv-4 */

/*  (4)  if iop1 = 0  c(nu-4,j) = dr(4) */
/*          iop1 = 1  c(nu-4,j) = dr(4) */
/*                    c(nu-5,j) = dr(4)+(dr(5)*cosi(1,j)+dr(6)*cosi(2,j)) */
/*                                *(tu(nu-4)-tu(nu-3))/3. = c1(j) */

/*  set constants */
/*<       one = 1 >*/
    /* Parameter adjustments */
    --nru;
    spu_dim1 = *mu;
    spu_offset = 1 + spu_dim1;
    spu -= spu_offset;
    --u;
    --nrv;
    a1 -= 3;
    a0 -= 3;
    spv_dim1 = *mv;
    spv_offset = 1 + spv_dim1;
    spv -= spv_offset;
    --v;
    --r__;
    --dr;
    bu_dim1 = *nu;
    bu_offset = 1 + bu_dim1;
    bu -= bu_offset;
    au_dim1 = *nu;
    au_offset = 1 + au_dim1;
    au -= au_offset;
    --fpu;
    --tu;
    cosi -= 3;
    --c1;
    --c0;
    b1 -= 3;
    b0 -= 3;
    bv_dim1 = *nv;
    bv_offset = 1 + bv_dim1;
    bv -= bv_offset;
    av2_dim1 = *nv;
    av2_offset = 1 + av2_dim1;
    av2 -= av2_offset;
    av1_dim1 = *nv;
    av1_offset = 1 + av1_dim1;
    av1 -= av1_offset;
    --fpv;
    --tv;
    --c__;
    --right;
    --q;

    /* Function Body */
    one = 1.;
/*<       three = 3 >*/
    three = 3.;
/*<       half = 0.5 >*/
    half = .5;
/*  initialization */
/*<       nu4 = nu-4 >*/
    nu4 = *nu - 4;
/*<       nu7 = nu-7 >*/
    nu7 = *nu - 7;
/*<       nu8 = nu-8 >*/
    nu8 = *nu - 8;
/*<       nu9 = nu-9 >*/
    nu9 = *nu - 9;
/*<       nv4 = nv-4 >*/
    nv4 = *nv - 4;
/*<       nv7 = nv-7 >*/
    nv7 = *nv - 7;
/*<       nv8 = nv-8 >*/
    nv8 = *nv - 8;
/*<       nv11 = nv-11 >*/
    nv11 = *nv - 11;
/*<       nuu = nu4-iop0-iop1-2 >*/
    nuu = nu4 - *iop0 - *iop1 - 2;
/*<       if(p.gt.0.) pinv = one/p >*/
    if (*p > 0.) {
	pinv = one / *p;
    }
/*  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv,iop0,iop1 */
/*  and on the value of p whether the matrices (spu), (spv), (bu), (bv), */
/*  (cosi) still must be determined. */
/*<       if(ifsu.ne.0) go to 30 >*/
    if (*ifsu != 0) {
	goto L30;
    }
/*  calculate the non-zero elements of the matrix (spu) which is the ob- */
/*  servation matrix according to the least-squares spline approximation */
/*  problem in the u-direction. */
/*<       l = 4 >*/
    l = 4;
/*<       l1 = 5 >*/
    l1 = 5;
/*<       number = 0 >*/
    number = 0;
/*<       do 25 it=1,mu >*/
    i__1 = *mu;
    for (it = 1; it <= i__1; ++it) {
/*<         arg = u(it) >*/
	arg = u[it];
/*<   10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 15 >*/
L10:
	if (arg < tu[l1] || l == nu4) {
	    goto L15;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         number = number+1 >*/
	++number;
/*<         go to 10 >*/
	goto L10;
/*<   15    call fpbspl(tu,nu,3,arg,l,h) >*/
L15:
	fpbspl_(&tu[1], nu, &c__3, &arg, &l, h__);
/*<         do 20 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           spu(it,i) = h(i) >*/
	    spu[it + i__ * spu_dim1] = h__[i__ - 1];
/*<   20    continue >*/
/* L20: */
	}
/*<         nru(it) = number >*/
	nru[it] = number;
/*<   25  continue >*/
/* L25: */
    }
/*<       ifsu = 1 >*/
    *ifsu = 1;
/*  calculate the non-zero elements of the matrix (spv) which is the ob- */
/*  servation matrix according to the least-squares spline approximation */
/*  problem in the v-direction. */
/*<   30  if(ifsv.ne.0) go to 85 >*/
L30:
    if (*ifsv != 0) {
	goto L85;
    }
/*<       l = 4 >*/
    l = 4;
/*<       l1 = 5 >*/
    l1 = 5;
/*<       number = 0 >*/
    number = 0;
/*<       do 50 it=1,mv >*/
    i__1 = *mv;
    for (it = 1; it <= i__1; ++it) {
/*<         arg = v(it) >*/
	arg = v[it];
/*<   35    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 40 >*/
L35:
	if (arg < tv[l1] || l == nv4) {
	    goto L40;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         number = number+1 >*/
	++number;
/*<         go to 35 >*/
	goto L35;
/*<   40    call fpbspl(tv,nv,3,arg,l,h) >*/
L40:
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
/*<         do 45 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           spv(it,i) = h(i) >*/
	    spv[it + i__ * spv_dim1] = h__[i__ - 1];
/*<   45    continue >*/
/* L45: */
	}
/*<         nrv(it) = number >*/
	nrv[it] = number;
/*<   50  continue >*/
/* L50: */
    }
/*<       ifsv = 1 >*/
    *ifsv = 1;
/*<       if(iop0.eq.0 .and. iop1.eq.0) go to 85 >*/
    if (*iop0 == 0 && *iop1 == 0) {
	goto L85;
    }
/*  calculate the coefficients of the interpolating splines for cos(v) */
/*  and sin(v). */
/*<       do 55 i=1,nv4 >*/
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          cosi(1,i) = 0. >*/
	cosi[(i__ << 1) + 1] = 0.;
/*<          cosi(2,i) = 0. >*/
	cosi[(i__ << 1) + 2] = 0.;
/*<   55  continue >*/
/* L55: */
    }
/*<       if(nv7.lt.4) go to 85 >*/
    if (nv7 < 4) {
	goto L85;
    }
/*<       do 65 i=1,nv7 >*/
    i__1 = nv7;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          l = i+3 >*/
	l = i__ + 3;
/*<          arg = tv(l) >*/
	arg = tv[l];
/*<          call fpbspl(tv,nv,3,arg,l,h) >*/
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
/*<          do 60 j=1,3 >*/
	for (j = 1; j <= 3; ++j) {
/*<             av1(i,j) = h(j) >*/
	    av1[i__ + j * av1_dim1] = h__[j - 1];
/*<   60     continue >*/
/* L60: */
	}
/*<          cosi(1,i) = cos(arg) >*/
	cosi[(i__ << 1) + 1] = cos(arg);
/*<          cosi(2,i) = sin(arg) >*/
	cosi[(i__ << 1) + 2] = sin(arg);
/*<   65  continue >*/
/* L65: */
    }
/*<       call fpcyt1(av1,nv7,nv) >*/
    fpcyt1_(&av1[av1_offset], &nv7, nv);
/*<       do 80 j=1,2 >*/
    for (j = 1; j <= 2; ++j) {
/*<          do 70 i=1,nv7 >*/
	i__1 = nv7;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             right(i) = cosi(j,i) >*/
	    right[i__] = cosi[j + (i__ << 1)];
/*<   70     continue >*/
/* L70: */
	}
/*<          call fpcyt2(av1,nv7,right,right,nv) >*/
	fpcyt2_(&av1[av1_offset], &nv7, &right[1], &right[1], nv);
/*<          do 75 i=1,nv7 >*/
	i__1 = nv7;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             cosi(j,i+1) = right(i) >*/
	    cosi[j + (i__ + 1 << 1)] = right[i__];
/*<   75     continue >*/
/* L75: */
	}
/*<          cosi(j,1) = cosi(j,nv7+1) >*/
	cosi[j + 2] = cosi[j + (nv7 + 1 << 1)];
/*<          cosi(j,nv7+2) = cosi(j,2) >*/
	cosi[j + (nv7 + 2 << 1)] = cosi[j + 4];
/*<          cosi(j,nv4) = cosi(j,3) >*/
	cosi[j + (nv4 << 1)] = cosi[j + 6];
/*<   80  continue >*/
/* L80: */
    }
/*<   85  if(p.le.0.) go to  150 >*/
L85:
    if (*p <= 0.) {
	goto L150;
    }
/*  calculate the non-zero elements of the matrix (bu). */
/*<       if(ifbu.ne.0 .or. nu8.eq.0) go to 90 >*/
    if (*ifbu != 0 || nu8 == 0) {
	goto L90;
    }
/*<       call fpdisc(tu,nu,5,bu,nu) >*/
    fpdisc_(&tu[1], nu, &c__5, &bu[bu_offset], nu);
/*<       ifbu = 1 >*/
    *ifbu = 1;
/*  calculate the non-zero elements of the matrix (bv). */
/*<   90  if(ifbv.ne.0 .or. nv8.eq.0) go to 150 >*/
L90:
    if (*ifbv != 0 || nv8 == 0) {
	goto L150;
    }
/*<       call fpdisc(tv,nv,5,bv,nv) >*/
    fpdisc_(&tv[1], nv, &c__5, &bv[bv_offset], nv);
/*<       ifbv = 1 >*/
    *ifbv = 1;
/*  substituting (2),(3) and (4) into (1), we obtain the overdetermined */
/*  system */
/*         (5)  (avv) (cc) (auu)' = (qq) */
/*  from which the nuu*nv7 remaining coefficients */
/*         c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7. */
/*  the elements of (cc), are then determined in the least-squares sense. */
/*  simultaneously, we compute the resulting sum of squared residuals sq. */
/*<  150  dr01 = dr(1) >*/
L150:
    dr01 = dr[1];
/*<       dr11 = dr(4) >*/
    dr11 = dr[4];
/*<       do 155 i=1,mv >*/
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          a0(1,i) = dr01 >*/
	a0[(i__ << 1) + 1] = dr01;
/*<          a1(1,i) = dr11 >*/
	a1[(i__ << 1) + 1] = dr11;
/*<  155  continue >*/
/* L155: */
    }
/*<       if(nv8.eq.0 .or. p.le.0.) go to 165 >*/
    if (nv8 == 0 || *p <= 0.) {
	goto L165;
    }
/*<       do 160 i=1,nv8 >*/
    i__1 = nv8;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          b0(1,i) = 0. >*/
	b0[(i__ << 1) + 1] = 0.;
/*<          b1(1,i) = 0. >*/
	b1[(i__ << 1) + 1] = 0.;
/*<  160  continue >*/
/* L160: */
    }
/*<  165  mvv = mv >*/
L165:
    mvv = *mv;
/*<       if(iop0.eq.0) go to 195 >*/
    if (*iop0 == 0) {
	goto L195;
    }
/*<       fac = (tu(5)-tu(4))/three >*/
    fac = (tu[5] - tu[4]) / three;
/*<       dr02 = dr(2)*fac >*/
    dr02 = dr[2] * fac;
/*<       dr03 = dr(3)*fac >*/
    dr03 = dr[3] * fac;
/*<       do 170 i=1,nv4 >*/
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          c0(i) = dr01+dr02*cosi(1,i)+dr03*cosi(2,i) >*/
	c0[i__] = dr01 + dr02 * cosi[(i__ << 1) + 1] + dr03 * cosi[(i__ << 1) 
		+ 2];
/*<  170  continue >*/
/* L170: */
    }
/*<       do 180 i=1,mv >*/
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          number = nrv(i) >*/
	number = nrv[i__];
/*<          fac = 0. >*/
	fac = 0.;
/*<          do 175 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             number = number+1 >*/
	    ++number;
/*<             fac = fac+c0(number)*spv(i,j) >*/
	    fac += c0[number] * spv[i__ + j * spv_dim1];
/*<  175     continue >*/
/* L175: */
	}
/*<          a0(2,i) = fac >*/
	a0[(i__ << 1) + 2] = fac;
/*<  180  continue >*/
/* L180: */
    }
/*<       if(nv8.eq.0 .or. p.le.0.) go to 195 >*/
    if (nv8 == 0 || *p <= 0.) {
	goto L195;
    }
/*<       do 190 i=1,nv8 >*/
    i__1 = nv8;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          number = i >*/
	number = i__;
/*<          fac = 0. >*/
	fac = 0.;
/*<          do 185 j=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             fac = fac+c0(number)*bv(i,j) >*/
	    fac += c0[number] * bv[i__ + j * bv_dim1];
/*<             number = number+1 >*/
	    ++number;
/*<  185     continue >*/
/* L185: */
	}
/*<          b0(2,i) = fac*pinv >*/
	b0[(i__ << 1) + 2] = fac * pinv;
/*<  190  continue >*/
/* L190: */
    }
/*<       mvv = mv+nv8 >*/
    mvv = *mv + nv8;
/*<  195  if(iop1.eq.0) go to 225 >*/
L195:
    if (*iop1 == 0) {
	goto L225;
    }
/*<       fac = (tu(nu4)-tu(nu4+1))/three >*/
    fac = (tu[nu4] - tu[nu4 + 1]) / three;
/*<       dr12 = dr(5)*fac >*/
    dr12 = dr[5] * fac;
/*<       dr13 = dr(6)*fac >*/
    dr13 = dr[6] * fac;
/*<       do 200 i=1,nv4 >*/
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          c1(i) = dr11+dr12*cosi(1,i)+dr13*cosi(2,i) >*/
	c1[i__] = dr11 + dr12 * cosi[(i__ << 1) + 1] + dr13 * cosi[(i__ << 1) 
		+ 2];
/*<  200  continue >*/
/* L200: */
    }
/*<       do 210 i=1,mv >*/
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          number = nrv(i) >*/
	number = nrv[i__];
/*<          fac = 0. >*/
	fac = 0.;
/*<          do 205 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             number = number+1 >*/
	    ++number;
/*<             fac = fac+c1(number)*spv(i,j) >*/
	    fac += c1[number] * spv[i__ + j * spv_dim1];
/*<  205     continue >*/
/* L205: */
	}
/*<          a1(2,i) = fac >*/
	a1[(i__ << 1) + 2] = fac;
/*<  210  continue >*/
/* L210: */
    }
/*<       if(nv8.eq.0 .or. p.le.0.) go to 225 >*/
    if (nv8 == 0 || *p <= 0.) {
	goto L225;
    }
/*<       do 220 i=1,nv8 >*/
    i__1 = nv8;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          number = i >*/
	number = i__;
/*<          fac = 0. >*/
	fac = 0.;
/*<          do 215 j=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<             fac = fac+c1(number)*bv(i,j) >*/
	    fac += c1[number] * bv[i__ + j * bv_dim1];
/*<             number = number+1 >*/
	    ++number;
/*<  215     continue >*/
/* L215: */
	}
/*<          b1(2,i) = fac*pinv >*/
	b1[(i__ << 1) + 2] = fac * pinv;
/*<  220  continue >*/
/* L220: */
    }
/*<       mvv = mv+nv8 >*/
    mvv = *mv + nv8;
/*  we first determine the matrices (auu) and (qq). then we reduce the */
/*  matrix (auu) to an unit upper triangular form (ru) using givens */
/*  rotations without square roots. we apply the same transformations to */
/*  the rows of matrix qq to obtain the mv x nuu matrix g. */
/*  we store matrix (ru) into au and g into q. */
/*<  225  l = mvv*nuu >*/
L225:
    l = mvv * nuu;
/*  initialization. */
/*<       sq = 0. >*/
    *sq = 0.;
/*<       if(l.eq.0) go to 245 >*/
    if (l == 0) {
	goto L245;
    }
/*<       do 230 i=1,l >*/
    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         q(i) = 0. >*/
	q[i__] = 0.;
/*<  230  continue >*/
/* L230: */
    }
/*<       do 240 i=1,nuu >*/
    i__1 = nuu;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         do 240 j=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<           au(i,j) = 0. >*/
	    au[i__ + j * au_dim1] = 0.;
/*<  240  continue >*/
/* L240: */
	}
    }
/*<       l = 0 >*/
    l = 0;
/*<  245  nrold = 0 >*/
L245:
    nrold = 0;
/*<       n1 = nrold+1 >*/
    n1 = nrold + 1;
/*<       do 420 it=1,mu >*/
    i__1 = *mu;
    for (it = 1; it <= i__1; ++it) {
/*<         number = nru(it) >*/
	number = nru[it];
/*  find the appropriate column of q. */
/*<  250    do 260 j=1,mvv >*/
L250:
	i__2 = mvv;
	for (j = 1; j <= i__2; ++j) {
/*<            right(j) = 0. >*/
	    right[j] = 0.;
/*<  260    continue >*/
/* L260: */
	}
/*<         if(nrold.eq.number) go to 280 >*/
	if (nrold == number) {
	    goto L280;
	}
/*<         if(p.le.0.) go to 410 >*/
	if (*p <= 0.) {
	    goto L410;
	}
/*  fetch a new row of matrix (bu). */
/*<         do 270 j=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<           h(j) = bu(n1,j)*pinv >*/
	    h__[j - 1] = bu[n1 + j * bu_dim1] * pinv;
/*<  270    continue >*/
/* L270: */
	}
/*<         i0 = 1 >*/
	i0 = 1;
/*<         i1 = 5 >*/
	i1 = 5;
/*<         go to 310 >*/
	goto L310;
/*  fetch a new row of matrix (spu). */
/*<  280    do 290 j=1,4 >*/
L280:
	for (j = 1; j <= 4; ++j) {
/*<           h(j) = spu(it,j) >*/
	    h__[j - 1] = spu[it + j * spu_dim1];
/*<  290    continue >*/
/* L290: */
	}
/*  find the appropriate column of q. */
/*<         do 300 j=1,mv >*/
	i__2 = *mv;
	for (j = 1; j <= i__2; ++j) {
/*<           l = l+1 >*/
	    ++l;
/*<           right(j) = r(l) >*/
	    right[j] = r__[l];
/*<  300    continue >*/
/* L300: */
	}
/*<         i0 = 1 >*/
	i0 = 1;
/*<         i1 = 4 >*/
	i1 = 4;
/*<  310    j0 = n1 >*/
L310:
	j0 = n1;
/*<         j1 = nu7-number >*/
	j1 = nu7 - number;
/*  take into account that we eliminate the constraints (3) */
/*<  315     if(j0-1.gt.iop0) go to 335 >*/
L315:
	if (j0 - 1 > *iop0) {
	    goto L335;
	}
/*<          fac0 = h(i0) >*/
	fac0 = h__[i0 - 1];
/*<          do 320 j=1,mv >*/
	i__2 = *mv;
	for (j = 1; j <= i__2; ++j) {
/*<             right(j) = right(j)-fac0*a0(j0,j) >*/
	    right[j] -= fac0 * a0[j0 + (j << 1)];
/*<  320     continue >*/
/* L320: */
	}
/*<          if(mv.eq.mvv) go to 330 >*/
	if (*mv == mvv) {
	    goto L330;
	}
/*<          j = mv >*/
	j = *mv;
/*<          do 325 jj=1,nv8 >*/
	i__2 = nv8;
	for (jj = 1; jj <= i__2; ++jj) {
/*<             j = j+1 >*/
	    ++j;
/*<             right(j) = right(j)-fac0*b0(j0,jj) >*/
	    right[j] -= fac0 * b0[j0 + (jj << 1)];
/*<  325     continue >*/
/* L325: */
	}
/*<  330     j0 = j0+1 >*/
L330:
	++j0;
/*<          i0 = i0+1 >*/
	++i0;
/*<          go to 315 >*/
	goto L315;
/*  take into account that we eliminate the constraints (4) */
/*<  335     if(j1-1.gt.iop1) go to 360 >*/
L335:
	if (j1 - 1 > *iop1) {
	    goto L360;
	}
/*<          fac1 = h(i1) >*/
	fac1 = h__[i1 - 1];
/*<          do 340 j=1,mv >*/
	i__2 = *mv;
	for (j = 1; j <= i__2; ++j) {
/*<             right(j) = right(j)-fac1*a1(j1,j) >*/
	    right[j] -= fac1 * a1[j1 + (j << 1)];
/*<  340     continue >*/
/* L340: */
	}
/*<          if(mv.eq.mvv) go to 350 >*/
	if (*mv == mvv) {
	    goto L350;
	}
/*<          j = mv >*/
	j = *mv;
/*<          do 345 jj=1,nv8 >*/
	i__2 = nv8;
	for (jj = 1; jj <= i__2; ++jj) {
/*<             j = j+1 >*/
	    ++j;
/*<             right(j) = right(j)-fac1*b1(j1,jj) >*/
	    right[j] -= fac1 * b1[j1 + (jj << 1)];
/*<  345     continue >*/
/* L345: */
	}
/*<  350     j1 = j1+1 >*/
L350:
	++j1;
/*<          i1 = i1-1 >*/
	--i1;
/*<          go to 335 >*/
	goto L335;
/*<  360     irot = nrold-iop0-1 >*/
L360:
	irot = nrold - *iop0 - 1;
/*<          if(irot.lt.0) irot = 0 >*/
	if (irot < 0) {
	    irot = 0;
	}
/*  rotate the new row of matrix (auu) into triangle. */
/*<         if(i0.gt.i1) go to 390 >*/
	if (i0 > i1) {
	    goto L390;
	}
/*<         do 385 i=i0,i1 >*/
	i__2 = i1;
	for (i__ = i0; i__ <= i__2; ++i__) {
/*<           irot = irot+1 >*/
	    ++irot;
/*<           piv = h(i) >*/
	    piv = h__[i__ - 1];
/*<           if(piv.eq.0.) go to 385 >*/
	    if (piv == 0.) {
		goto L385;
	    }
/*  calculate the parameters of the givens transformation. */
/*<           call fpgivs(piv,au(irot,1),co,si) >*/
	    fpgivs_(&piv, &au[irot + au_dim1], &co, &si);
/*  apply that transformation to the rows of matrix (qq). */
/*<           iq = (irot-1)*mvv >*/
	    iq = (irot - 1) * mvv;
/*<           do 370 j=1,mvv >*/
	    i__3 = mvv;
	    for (j = 1; j <= i__3; ++j) {
/*<             iq = iq+1 >*/
		++iq;
/*<             call fprota(co,si,right(j),q(iq)) >*/
		fprota_(&co, &si, &right[j], &q[iq]);
/*<  370      continue >*/
/* L370: */
	    }
/*  apply that transformation to the columns of (auu). */
/*<           if(i.eq.i1) go to 385 >*/
	    if (i__ == i1) {
		goto L385;
	    }
/*<           i2 = 1 >*/
	    i2 = 1;
/*<           i3 = i+1 >*/
	    i3 = i__ + 1;
/*<           do 380 j=i3,i1 >*/
	    i__3 = i1;
	    for (j = i3; j <= i__3; ++j) {
/*<             i2 = i2+1 >*/
		++i2;
/*<             call fprota(co,si,h(j),au(irot,i2)) >*/
		fprota_(&co, &si, &h__[j - 1], &au[irot + i2 * au_dim1]);
/*<  380      continue >*/
/* L380: */
	    }
/*<  385    continue >*/
L385:
	    ;
	}
/*  we update the sum of squared residuals. */
/*<  390    do 395 j=1,mvv >*/
L390:
	i__2 = mvv;
	for (j = 1; j <= i__2; ++j) {
/*<           sq = sq+right(j)**2 >*/
/* Computing 2nd power */
	    d__1 = right[j];
	    *sq += d__1 * d__1;
/*<  395    continue >*/
/* L395: */
	}
/*<  400    if(nrold.eq.number) go to 420 >*/
/* L400: */
	if (nrold == number) {
	    goto L420;
	}
/*<  410    nrold = n1 >*/
L410:
	nrold = n1;
/*<         n1 = n1+1 >*/
	++n1;
/*<         go to 250 >*/
	goto L250;
/*<  420  continue >*/
L420:
	;
    }
/*<       if(nuu.eq.0) go to 800 >*/
    if (nuu == 0) {
	goto L800;
    }
/*  we determine the matrix (avv) and then we reduce her to an unit */
/*  upper triangular form (rv) using givens rotations without square */
/*  roots. we apply the same transformations to the columns of matrix */
/*  g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h. */
/*  we store matrix (rv) into av1 and av2, h into c. */
/*  the nv7 x nv7 triangular unit upper matrix (rv) has the form */
/*              | av1 '     | */
/*       (rv) = |     ' av2 | */
/*              |  0  '     | */
/*  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper */
/*  triangular matrix of bandwidth 5. */
/*<       ncof = nuu*nv7 >*/
    ncof = nuu * nv7;
/*  initialization. */
/*<       do 430 i=1,ncof >*/
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         c(i) = 0. >*/
	c__[i__] = 0.;
/*<  430  continue >*/
/* L430: */
    }
/*<       do 440 i=1,nv4 >*/
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         av1(i,5) = 0. >*/
	av1[i__ + av1_dim1 * 5] = 0.;
/*<         do 440 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           av1(i,j) = 0. >*/
	    av1[i__ + j * av1_dim1] = 0.;
/*<           av2(i,j) = 0. >*/
	    av2[i__ + j * av2_dim1] = 0.;
/*<  440  continue >*/
/* L440: */
	}
    }
/*<       jper = 0 >*/
    jper = 0;
/*<       nrold = 0 >*/
    nrold = 0;
/*<       do 770 it=1,mv >*/
    i__1 = *mv;
    for (it = 1; it <= i__1; ++it) {
/*<         number = nrv(it) >*/
	number = nrv[it];
/*<  450    if(nrold.eq.number) go to 480 >*/
L450:
	if (nrold == number) {
	    goto L480;
	}
/*<         if(p.le.0.) go to 760 >*/
	if (*p <= 0.) {
	    goto L760;
	}
/*  fetch a new row of matrix (bv). */
/*<         n1 = nrold+1 >*/
	n1 = nrold + 1;
/*<         do 460 j=1,5 >*/
	for (j = 1; j <= 5; ++j) {
/*<           h(j) = bv(n1,j)*pinv >*/
	    h__[j - 1] = bv[n1 + j * bv_dim1] * pinv;
/*<  460    continue >*/
/* L460: */
	}
/*  find the appropiate row of g. */
/*<         do 465 j=1,nuu >*/
	i__2 = nuu;
	for (j = 1; j <= i__2; ++j) {
/*<           right(j) = 0. >*/
	    right[j] = 0.;
/*<  465    continue >*/
/* L465: */
	}
/*<         if(mv.eq.mvv) go to 510 >*/
	if (*mv == mvv) {
	    goto L510;
	}
/*<         l = mv+n1 >*/
	l = *mv + n1;
/*<         do 470 j=1,nuu >*/
	i__2 = nuu;
	for (j = 1; j <= i__2; ++j) {
/*<           right(j) = q(l) >*/
	    right[j] = q[l];
/*<           l = l+mvv >*/
	    l += mvv;
/*<  470    continue >*/
/* L470: */
	}
/*<         go to 510 >*/
	goto L510;
/*  fetch a new row of matrix (spv) */
/*<  480    h(5) = 0. >*/
L480:
	h__[4] = 0.;
/*<         do 490 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<           h(j) = spv(it,j) >*/
	    h__[j - 1] = spv[it + j * spv_dim1];
/*<  490    continue >*/
/* L490: */
	}
/*  find the appropiate row of g. */
/*<         l = it >*/
	l = it;
/*<         do 500 j=1,nuu >*/
	i__2 = nuu;
	for (j = 1; j <= i__2; ++j) {
/*<           right(j) = q(l) >*/
	    right[j] = q[l];
/*<           l = l+mvv >*/
	    l += mvv;
/*<  500    continue >*/
/* L500: */
	}
/*  test whether there are non-zero values in the new row of (avv) */
/*  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4. */
/*<  510     if(nrold.lt.nv11) go to 710 >*/
L510:
	if (nrold < nv11) {
	    goto L710;
	}
/*<          if(jper.ne.0) go to 550 >*/
	if (jper != 0) {
	    goto L550;
	}
/*  initialize the matrix (av2). */
/*<          jk = nv11+1 >*/
	jk = nv11 + 1;
/*<          do 540 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<             ik = jk >*/
	    ik = jk;
/*<             do 520 j=1,5 >*/
	    for (j = 1; j <= 5; ++j) {
/*<                if(ik.le.0) go to 530 >*/
		if (ik <= 0) {
		    goto L530;
		}
/*<                av2(ik,i) = av1(ik,j) >*/
		av2[ik + i__ * av2_dim1] = av1[ik + j * av1_dim1];
/*<                ik = ik-1 >*/
		--ik;
/*<  520        continue >*/
/* L520: */
	    }
/*<  530        jk = jk+1 >*/
L530:
	    ++jk;
/*<  540     continue >*/
/* L540: */
	}
/*<          jper = 1 >*/
	jper = 1;
/*  if one of the non-zero elements of the new row corresponds to one of */
/*  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition */
/*  (2) for setting up this row of (avv). the row is stored in h1( the */
/*  part with respect to av1) and h2 (the part with respect to av2). */
/*<  550     do 560 i=1,4 >*/
L550:
	for (i__ = 1; i__ <= 4; ++i__) {
/*<             h1(i) = 0. >*/
	    h1[i__ - 1] = 0.;
/*<             h2(i) = 0. >*/
	    h2[i__ - 1] = 0.;
/*<  560     continue >*/
/* L560: */
	}
/*<          h1(5) = 0. >*/
	h1[4] = 0.;
/*<          j = nrold-nv11 >*/
	j = nrold - nv11;
/*<          do 600 i=1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<             j = j+1 >*/
	    ++j;
/*<             l0 = j >*/
	    l0 = j;
/*<  570        l1 = l0-4 >*/
L570:
	    l1 = l0 - 4;
/*<             if(l1.le.0) go to 590 >*/
	    if (l1 <= 0) {
		goto L590;
	    }
/*<             if(l1.le.nv11) go to 580 >*/
	    if (l1 <= nv11) {
		goto L580;
	    }
/*<             l0 = l1-nv11 >*/
	    l0 = l1 - nv11;
/*<             go to 570 >*/
	    goto L570;
/*<  580        h1(l1) = h(i) >*/
L580:
	    h1[l1 - 1] = h__[i__ - 1];
/*<             go to 600 >*/
	    goto L600;
/*<  590        h2(l0) = h2(l0) + h(i) >*/
L590:
	    h2[l0 - 1] += h__[i__ - 1];
/*<  600     continue >*/
L600:
	    ;
	}
/*  rotate the new row of (avv) into triangle. */
/*<          if(nv11.le.0) go to 670 >*/
	if (nv11 <= 0) {
	    goto L670;
	}
/*  rotations with the rows 1,2,...,nv11 of (avv). */
/*<          do 660 j=1,nv11 >*/
	i__2 = nv11;
	for (j = 1; j <= i__2; ++j) {
/*<             piv = h1(1) >*/
	    piv = h1[0];
/*<             i2 = min0(nv11-j,4) >*/
/* Computing MIN */
	    i__3 = nv11 - j;
	    i2 = min(i__3,4);
/*<             if(piv.eq.0.) go to 640 >*/
	    if (piv == 0.) {
		goto L640;
	    }
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,av1(j,1),co,si) >*/
	    fpgivs_(&piv, &av1[j + av1_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
/*<             ic = j >*/
	    ic = j;
/*<             do 610 i=1,nuu >*/
	    i__3 = nuu;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<                call fprota(co,si,right(i),c(ic)) >*/
		fprota_(&co, &si, &right[i__], &c__[ic]);
/*<                ic = ic+nv7 >*/
		ic += nv7;
/*<  610        continue >*/
/* L610: */
	    }
/*  apply that transformation to the rows of (avv) with respect to av2. */
/*<             do 620 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<                call fprota(co,si,h2(i),av2(j,i)) >*/
		fprota_(&co, &si, &h2[i__ - 1], &av2[j + i__ * av2_dim1]);
/*<  620        continue >*/
/* L620: */
	    }
/*  apply that transformation to the rows of (avv) with respect to av1. */
/*<             if(i2.eq.0) go to 670 >*/
	    if (i2 == 0) {
		goto L670;
	    }
/*<             do 630 i=1,i2 >*/
	    i__3 = i2;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<                i1 = i+1 >*/
		i1 = i__ + 1;
/*<                call fprota(co,si,h1(i1),av1(j,i1)) >*/
		fprota_(&co, &si, &h1[i1 - 1], &av1[j + i1 * av1_dim1]);
/*<  630        continue >*/
/* L630: */
	    }
/*<  640        do 650 i=1,i2 >*/
L640:
	    i__3 = i2;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<                h1(i) = h1(i+1) >*/
		h1[i__ - 1] = h1[i__];
/*<  650        continue >*/
/* L650: */
	    }
/*<             h1(i2+1) = 0. >*/
	    h1[i2] = 0.;
/*<  660     continue >*/
/* L660: */
	}
/*  rotations with the rows nv11+1,...,nv7 of avv. */
/*<  670     do 700 j=1,4 >*/
L670:
	for (j = 1; j <= 4; ++j) {
/*<             ij = nv11+j >*/
	    ij = nv11 + j;
/*<             if(ij.le.0) go to 700 >*/
	    if (ij <= 0) {
		goto L700;
	    }
/*<             piv = h2(j) >*/
	    piv = h2[j - 1];
/*<             if(piv.eq.0.) go to 700 >*/
	    if (piv == 0.) {
		goto L700;
	    }
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,av2(ij,j),co,si) >*/
	    fpgivs_(&piv, &av2[ij + j * av2_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
/*<             ic = ij >*/
	    ic = ij;
/*<             do 680 i=1,nuu >*/
	    i__2 = nuu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/*<                call fprota(co,si,right(i),c(ic)) >*/
		fprota_(&co, &si, &right[i__], &c__[ic]);
/*<                ic = ic+nv7 >*/
		ic += nv7;
/*<  680        continue >*/
/* L680: */
	    }
/*<             if(j.eq.4) go to 700 >*/
	    if (j == 4) {
		goto L700;
	    }
/*  apply that transformation to the rows of (avv) with respect to av2. */
/*<             j1 = j+1 >*/
	    j1 = j + 1;
/*<             do 690 i=j1,4 >*/
	    for (i__ = j1; i__ <= 4; ++i__) {
/*<                call fprota(co,si,h2(i),av2(ij,i)) >*/
		fprota_(&co, &si, &h2[i__ - 1], &av2[ij + i__ * av2_dim1]);
/*<  690        continue >*/
/* L690: */
	    }
/*<  700     continue >*/
L700:
	    ;
	}
/*  we update the sum of squared residuals. */
/*<          do 705 i=1,nuu >*/
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<            sq = sq+right(i)**2 >*/
/* Computing 2nd power */
	    d__1 = right[i__];
	    *sq += d__1 * d__1;
/*<  705     continue >*/
/* L705: */
	}
/*<          go to 750 >*/
	goto L750;
/*  rotation into triangle of the new row of (avv), in case the elements */
/*  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero. */
/*<  710     irot =nrold >*/
L710:
	irot = nrold;
/*<          do 740 i=1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<             irot = irot+1 >*/
	    ++irot;
/*<             piv = h(i) >*/
	    piv = h__[i__ - 1];
/*<             if(piv.eq.0.) go to 740 >*/
	    if (piv == 0.) {
		goto L740;
	    }
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,av1(irot,1),co,si) >*/
	    fpgivs_(&piv, &av1[irot + av1_dim1], &co, &si);
/*  apply that transformation to the columns of matrix g. */
/*<             ic = irot >*/
	    ic = irot;
/*<             do 720 j=1,nuu >*/
	    i__2 = nuu;
	    for (j = 1; j <= i__2; ++j) {
/*<                call fprota(co,si,right(j),c(ic)) >*/
		fprota_(&co, &si, &right[j], &c__[ic]);
/*<                ic = ic+nv7 >*/
		ic += nv7;
/*<  720        continue >*/
/* L720: */
	    }
/*  apply that transformation to the rows of (avv). */
/*<             if(i.eq.5) go to 740 >*/
	    if (i__ == 5) {
		goto L740;
	    }
/*<             i2 = 1 >*/
	    i2 = 1;
/*<             i3 = i+1 >*/
	    i3 = i__ + 1;
/*<             do 730 j=i3,5 >*/
	    for (j = i3; j <= 5; ++j) {
/*<                i2 = i2+1 >*/
		++i2;
/*<                call fprota(co,si,h(j),av1(irot,i2)) >*/
		fprota_(&co, &si, &h__[j - 1], &av1[irot + i2 * av1_dim1]);
/*<  730        continue >*/
/* L730: */
	    }
/*<  740     continue >*/
L740:
	    ;
	}
/*  we update the sum of squared residuals. */
/*<          do 745 i=1,nuu >*/
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<            sq = sq+right(i)**2 >*/
/* Computing 2nd power */
	    d__1 = right[i__];
	    *sq += d__1 * d__1;
/*<  745     continue >*/
/* L745: */
	}
/*<  750     if(nrold.eq.number) go to 770 >*/
L750:
	if (nrold == number) {
	    goto L770;
	}
/*<  760     nrold = nrold+1 >*/
L760:
	++nrold;
/*<          go to 450 >*/
	goto L450;
/*<  770  continue >*/
L770:
	;
    }
/*  test whether the b-spline coefficients must be determined. */
/*<       if(iback.ne.0) return >*/
    if (*iback != 0) {
	return 0;
    }
/*  backward substitution to obtain the b-spline coefficients as the */
/*  solution of the linear system    (rv) (cr) (ru)' = h. */
/*  first step: solve the system  (rv) (c1) = h. */
/*<       k = 1 >*/
    k = 1;
/*<       do 780 i=1,nuu >*/
    i__1 = nuu;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          call fpbacp(av1,av2,c(k),nv7,4,c(k),5,nv) >*/
	fpbacp_(&av1[av1_offset], &av2[av2_offset], &c__[k], &nv7, &c__4, &
		c__[k], &c__5, nv);
/*<          k = k+nv7 >*/
	k += nv7;
/*<  780  continue >*/
/* L780: */
    }
/*  second step: solve the system  (cr) (ru)' = (c1). */
/*<       k = 0 >*/
    k = 0;
/*<       do 795 j=1,nv7 >*/
    i__1 = nv7;
    for (j = 1; j <= i__1; ++j) {
/*<         k = k+1 >*/
	++k;
/*<         l = k >*/
	l = k;
/*<         do 785 i=1,nuu >*/
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           right(i) = c(l) >*/
	    right[i__] = c__[l];
/*<           l = l+nv7 >*/
	    l += nv7;
/*<  785    continue >*/
/* L785: */
	}
/*<         call fpback(au,right,nuu,5,right,nu) >*/
	fpback_(&au[au_offset], &right[1], &nuu, &c__5, &right[1], nu);
/*<         l = k >*/
	l = k;
/*<         do 790 i=1,nuu >*/
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<            c(l) = right(i) >*/
	    c__[l] = right[i__];
/*<            l = l+nv7 >*/
	    l += nv7;
/*<  790    continue >*/
/* L790: */
	}
/*<  795  continue >*/
/* L795: */
    }
/*  calculate from the conditions (2)-(3)-(4), the remaining b-spline */
/*  coefficients. */
/*<  800  ncof = nu4*nv4 >*/
L800:
    ncof = nu4 * nv4;
/*<       j = ncof >*/
    j = ncof;
/*<       do 805 l=1,nv4 >*/
    i__1 = nv4;
    for (l = 1; l <= i__1; ++l) {
/*<          q(l) = dr01 >*/
	q[l] = dr01;
/*<          q(j) = dr11 >*/
	q[j] = dr11;
/*<          j = j-1 >*/
	--j;
/*<  805  continue >*/
/* L805: */
    }
/*<       i = nv4 >*/
    i__ = nv4;
/*<       j = 0 >*/
    j = 0;
/*<       if(iop0.eq.0) go to 815 >*/
    if (*iop0 == 0) {
	goto L815;
    }
/*<       do 810 l=1,nv4 >*/
    i__1 = nv4;
    for (l = 1; l <= i__1; ++l) {
/*<          i = i+1 >*/
	++i__;
/*<          q(i) = c0(l) >*/
	q[i__] = c0[l];
/*<  810  continue >*/
/* L810: */
    }
/*<  815  if(nuu.eq.0) go to 835 >*/
L815:
    if (nuu == 0) {
	goto L835;
    }
/*<       do 830 l=1,nuu >*/
    i__1 = nuu;
    for (l = 1; l <= i__1; ++l) {
/*<          ii = i >*/
	ii = i__;
/*<          do 820 k=1,nv7 >*/
	i__2 = nv7;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             j = j+1 >*/
	    ++j;
/*<             q(i) = c(j) >*/
	    q[i__] = c__[j];
/*<  820     continue >*/
/* L820: */
	}
/*<          do 825 k=1,3 >*/
	for (k = 1; k <= 3; ++k) {
/*<             ii = ii+1 >*/
	    ++ii;
/*<             i = i+1 >*/
	    ++i__;
/*<             q(i) = q(ii) >*/
	    q[i__] = q[ii];
/*<  825     continue >*/
/* L825: */
	}
/*<  830  continue >*/
/* L830: */
    }
/*<  835  if(iop1.eq.0) go to 845 >*/
L835:
    if (*iop1 == 0) {
	goto L845;
    }
/*<       do 840 l=1,nv4 >*/
    i__1 = nv4;
    for (l = 1; l <= i__1; ++l) {
/*<          i = i+1 >*/
	++i__;
/*<          q(i) = c1(l) >*/
	q[i__] = c1[l];
/*<  840  continue >*/
/* L840: */
    }
/*<  845  do 850 i=1,ncof >*/
L845:
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          c(i) = q(i) >*/
	c__[i__] = q[i__];
/*<  850  continue >*/
/* L850: */
    }
/*  calculate the quantities */
/*    res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv */
/*    fp = sumi=1,mu(sumj=1,mv(res(i,j))) */
/*    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7 */
/*                  tu(r+3) <= u(i) <= tu(r+4) */
/*    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7 */
/*                  tv(r+3) <= v(j) <= tv(r+4) */
/*<       fp = 0. >*/
    *fp = 0.;
/*<       do 890 i=1,nu >*/
    i__1 = *nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         fpu(i) = 0. >*/
	fpu[i__] = 0.;
/*<  890  continue >*/
/* L890: */
    }
/*<       do 900 i=1,nv >*/
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         fpv(i) = 0. >*/
	fpv[i__] = 0.;
/*<  900  continue >*/
/* L900: */
    }
/*<       ir = 0 >*/
    ir = 0;
/*<       nroldu = 0 >*/
    nroldu = 0;
/*  main loop for the different grid points. */
/*<       do 950 i1=1,mu >*/
    i__1 = *mu;
    for (i1 = 1; i1 <= i__1; ++i1) {
/*<         numu = nru(i1) >*/
	numu = nru[i1];
/*<         numu1 = numu+1 >*/
	numu1 = numu + 1;
/*<         nroldv = 0 >*/
	nroldv = 0;
/*<         do 940 i2=1,mv >*/
	i__2 = *mv;
	for (i2 = 1; i2 <= i__2; ++i2) {
/*<           numv = nrv(i2) >*/
	    numv = nrv[i2];
/*<           numv1 = numv+1 >*/
	    numv1 = numv + 1;
/*<           ir = ir+1 >*/
	    ++ir;
/*  evaluate s(u,v) at the current grid point by making the sum of the */
/*  cross products of the non-zero b-splines at (u,v), multiplied with */
/*  the appropiate b-spline coefficients. */
/*<           term = 0. >*/
	    term = 0.;
/*<           k1 = numu*nv4+numv >*/
	    k1 = numu * nv4 + numv;
/*<           do 920 l1=1,4 >*/
	    for (l1 = 1; l1 <= 4; ++l1) {
/*<             k2 = k1 >*/
		k2 = k1;
/*<             fac = spu(i1,l1) >*/
		fac = spu[i1 + l1 * spu_dim1];
/*<             do 910 l2=1,4 >*/
		for (l2 = 1; l2 <= 4; ++l2) {
/*<               k2 = k2+1 >*/
		    ++k2;
/*<               term = term+fac*spv(i2,l2)*c(k2) >*/
		    term += fac * spv[i2 + l2 * spv_dim1] * c__[k2];
/*<  910        continue >*/
/* L910: */
		}
/*<             k1 = k1+nv4 >*/
		k1 += nv4;
/*<  920      continue >*/
/* L920: */
	    }
/*  calculate the squared residual at the current grid point. */
/*<           term = (r(ir)-term)**2 >*/
/* Computing 2nd power */
	    d__1 = r__[ir] - term;
	    term = d__1 * d__1;
/*  adjust the different parameters. */
/*<           fp = fp+term >*/
	    *fp += term;
/*<           fpu(numu1) = fpu(numu1)+term >*/
	    fpu[numu1] += term;
/*<           fpv(numv1) = fpv(numv1)+term >*/
	    fpv[numv1] += term;
/*<           fac = term*half >*/
	    fac = term * half;
/*<           if(numv.eq.nroldv) go to 930 >*/
	    if (numv == nroldv) {
		goto L930;
	    }
/*<           fpv(numv1) = fpv(numv1)-fac >*/
	    fpv[numv1] -= fac;
/*<           fpv(numv) = fpv(numv)+fac >*/
	    fpv[numv] += fac;
/*<  930      nroldv = numv >*/
L930:
	    nroldv = numv;
/*<           if(numu.eq.nroldu) go to 940 >*/
	    if (numu == nroldu) {
		goto L940;
	    }
/*<           fpu(numu1) = fpu(numu1)-fac >*/
	    fpu[numu1] -= fac;
/*<           fpu(numu) = fpu(numu)+fac >*/
	    fpu[numu] += fac;
/*<  940    continue >*/
L940:
	    ;
	}
/*<         nroldu = numu >*/
	nroldu = numu;
/*<  950  continue >*/
/* L950: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpgrsp_ */

