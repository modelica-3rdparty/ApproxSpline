/* fpgrpa.f -- translated by f2c (version 20061008).
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
/* Subroutine */ int fpgrpa_(integer *ifsu, integer *ifsv, integer *ifbu, 
	integer *ifbv, integer *idim, integer *ipar, doublereal *u, integer *
	mu, doublereal *v, integer *mv, doublereal *z__, integer *mz, 
	doublereal *tu, integer *nu, doublereal *tv, integer *nv, doublereal *
	p, doublereal *c__, integer *nc, doublereal *fp, doublereal *fpu, 
	doublereal *fpv, integer *mm, integer *mvnu, doublereal *spu, 
	doublereal *spv, doublereal *right, doublereal *q, doublereal *au, 
	doublereal *au1, doublereal *av, doublereal *av1, doublereal *bu, 
	doublereal *bv, integer *nru, integer *nrv)
{
    /* System generated locals */
    integer spu_dim1, spu_offset, spv_dim1, spv_offset, au_dim1, au_offset, 
	    au1_dim1, au1_offset, av_dim1, av_offset, av1_dim1, av1_offset, 
	    bu_dim1, bu_offset, bv_dim1, bv_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal h__[5];
    static integer i__, j, k, l, i1, i2, k1, k2, l1, l2, k0, id, ii, n33, it, 
	    iz, jz, nu4, nv4, nu7, nu8, nv7, nv8;
    static doublereal fac, arg;
    static integer nmd;
    static doublereal one;
    static integer muu, mvv, nuu, nvv;
    static doublereal half;
    static integer ncof;
    static doublereal term;
    static integer numu, numv, numu1, numv1;
    static doublereal value;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), fpbacp_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), fpdisc_(doublereal *, integer *, integer *,
	     doublereal *, integer *), fpbspl_(doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer number;
    extern /* Subroutine */ int fptrpe_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer nroldu, nroldv;
    extern /* Subroutine */ int fptrnp_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/*  .. */
/*  ..scalar arguments.. */
/*<       real p,fp >*/
/*<       integer ifsu,ifsv,ifbu,ifbv,idim,mu,mv,mz,nu,nv,nc,mm,mvnu >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer ipar(2),nru(mu),nrv(mv) >*/
/*  ..local scalars.. */
/*<       real arg,fac,term,one,half,value >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real h(5) >*/
/*  ..subroutine references.. */
/*    fpback,fpbspl,fpdisc,fpbacp,fptrnp,fptrpe */
/*  .. */
/*  let */
/*               |   (spu)    |            |   (spv)    | */
/*        (au) = | ---------- |     (av) = | ---------- | */
/*               | (1/p) (bu) |            | (1/p) (bv) | */

/*                                | z  ' 0 | */
/*                            q = | ------ | */
/*                                | 0  ' 0 | */

/*  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline */
/*                coefficients. */
/*       z      : the mu x mv matrix which contains the function values. */
/*       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices */
/*                according to the least-squares problems in the u-,resp. */
/*                v-direction. */
/*       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices */
/*                containing the discontinuity jumps of the derivatives */
/*                of the b-splines in the u-,resp.v-variable at the knots */
/*  the b-spline coefficients of the smoothing spline are then calculated */
/*  as the least-squares solution of the following over-determined linear */
/*  system of equations */

/*    (1)  (av) c (au)' = q */

/*  subject to the constraints */

/*    (2)  c(nu-3+i,j) = c(i,j), i=1,2,3 ; j=1,2,...,nv-4 */
/*            if(ipar(1).ne.0) */

/*    (3)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4 */
/*            if(ipar(2).ne.0) */

/*  set constants */
/*<       one = 1 >*/
    /* Parameter adjustments */
    --ipar;
    --nru;
    spu_dim1 = *mu;
    spu_offset = 1 + spu_dim1;
    spu -= spu_offset;
    --u;
    --nrv;
    spv_dim1 = *mv;
    spv_offset = 1 + spv_dim1;
    spv -= spv_offset;
    --v;
    --z__;
    bu_dim1 = *nu;
    bu_offset = 1 + bu_dim1;
    bu -= bu_offset;
    au1_dim1 = *nu;
    au1_offset = 1 + au1_dim1;
    au1 -= au1_offset;
    au_dim1 = *nu;
    au_offset = 1 + au_dim1;
    au -= au_offset;
    --fpu;
    --tu;
    bv_dim1 = *nv;
    bv_offset = 1 + bv_dim1;
    bv -= bv_offset;
    av1_dim1 = *nv;
    av1_offset = 1 + av1_dim1;
    av1 -= av1_offset;
    av_dim1 = *nv;
    av_offset = 1 + av_dim1;
    av -= av_offset;
    --fpv;
    --tv;
    --c__;
    --right;
    --q;

    /* Function Body */
    one = 1.;
/*<       half = 0.5 >*/
    half = .5;
/*  initialization */
/*<       nu4 = nu-4 >*/
    nu4 = *nu - 4;
/*<       nu7 = nu-7 >*/
    nu7 = *nu - 7;
/*<       nu8 = nu-8 >*/
    nu8 = *nu - 8;
/*<       nv4 = nv-4 >*/
    nv4 = *nv - 4;
/*<       nv7 = nv-7 >*/
    nv7 = *nv - 7;
/*<       nv8 = nv-8 >*/
    nv8 = *nv - 8;
/*<       muu = mu >*/
    muu = *mu;
/*<       if(ipar(1).ne.0) muu = mu-1 >*/
    if (ipar[1] != 0) {
	muu = *mu - 1;
    }
/*<       mvv = mv >*/
    mvv = *mv;
/*<       if(ipar(2).ne.0) mvv = mv-1 >*/
    if (ipar[2] != 0) {
	mvv = *mv - 1;
    }
/*  it depends on the value of the flags ifsu,ifsv,ifbu  and ibvand */
/*  on the value of p whether the matrices (spu), (spv), (bu) and (bv) */
/*  still must be determined. */
/*<       if(ifsu.ne.0) go to 50 >*/
    if (*ifsu != 0) {
	goto L50;
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
/*<       do 40 it=1,muu >*/
    i__1 = muu;
    for (it = 1; it <= i__1; ++it) {
/*<         arg = u(it) >*/
	arg = u[it];
/*<   10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20 >*/
L10:
	if (arg < tu[l1] || l == nu4) {
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
/*<   20    call fpbspl(tu,nu,3,arg,l,h) >*/
L20:
	fpbspl_(&tu[1], nu, &c__3, &arg, &l, h__);
/*<         do 30 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           spu(it,i) = h(i) >*/
	    spu[it + i__ * spu_dim1] = h__[i__ - 1];
/*<   30    continue >*/
/* L30: */
	}
/*<         nru(it) = number >*/
	nru[it] = number;
/*<   40  continue >*/
/* L40: */
    }
/*<       ifsu = 1 >*/
    *ifsu = 1;
/*  calculate the non-zero elements of the matrix (spv) which is the ob- */
/*  servation matrix according to the least-squares spline approximation */
/*  problem in the v-direction. */
/*<   50  if(ifsv.ne.0) go to 100 >*/
L50:
    if (*ifsv != 0) {
	goto L100;
    }
/*<       l = 4 >*/
    l = 4;
/*<       l1 = 5 >*/
    l1 = 5;
/*<       number = 0 >*/
    number = 0;
/*<       do 90 it=1,mvv >*/
    i__1 = mvv;
    for (it = 1; it <= i__1; ++it) {
/*<         arg = v(it) >*/
	arg = v[it];
/*<   60    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 70 >*/
L60:
	if (arg < tv[l1] || l == nv4) {
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
/*<   70    call fpbspl(tv,nv,3,arg,l,h) >*/
L70:
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
/*<         do 80 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           spv(it,i) = h(i) >*/
	    spv[it + i__ * spv_dim1] = h__[i__ - 1];
/*<   80    continue >*/
/* L80: */
	}
/*<         nrv(it) = number >*/
	nrv[it] = number;
/*<   90  continue >*/
/* L90: */
    }
/*<       ifsv = 1 >*/
    *ifsv = 1;
/*<  100  if(p.le.0.) go to  150 >*/
L100:
    if (*p <= 0.) {
	goto L150;
    }
/*  calculate the non-zero elements of the matrix (bu). */
/*<       if(ifbu.ne.0 .or. nu8.eq.0) go to 110 >*/
    if (*ifbu != 0 || nu8 == 0) {
	goto L110;
    }
/*<       call fpdisc(tu,nu,5,bu,nu) >*/
    fpdisc_(&tu[1], nu, &c__5, &bu[bu_offset], nu);
/*<       ifbu = 1 >*/
    *ifbu = 1;
/*  calculate the non-zero elements of the matrix (bv). */
/*<  110  if(ifbv.ne.0 .or. nv8.eq.0) go to 150 >*/
L110:
    if (*ifbv != 0 || nv8 == 0) {
	goto L150;
    }
/*<       call fpdisc(tv,nv,5,bv,nv) >*/
    fpdisc_(&tv[1], nv, &c__5, &bv[bv_offset], nv);
/*<       ifbv = 1 >*/
    *ifbv = 1;
/*  substituting (2)  and (3) into (1), we obtain the overdetermined */
/*  system */
/*         (4)  (avv) (cr) (auu)' = (qq) */
/*  from which the nuu*nvv remaining coefficients */
/*         c(i,j) , i=1,...,nu-4-3*ipar(1) ; j=1,...,nv-4-3*ipar(2) , */
/*  the elements of (cr), are then determined in the least-squares sense. */
/*  we first determine the matrices (auu) and (qq). then we reduce the */
/*  matrix (auu) to upper triangular form (ru) using givens rotations. */
/*  we apply the same transformations to the rows of matrix qq to obtain */
/*  the (mv) x nuu matrix g. */
/*  we store matrix (ru) into au (and au1 if ipar(1)=1) and g into q. */
/*<  150  if(ipar(1).ne.0) go to 160 >*/
L150:
    if (ipar[1] != 0) {
	goto L160;
    }
/*<       nuu = nu4 >*/
    nuu = nu4;
/*<       call fptrnp(mu,mv,idim,nu,nru,spu,p,bu,z,au,q,right) >*/
    fptrnp_(mu, mv, idim, nu, &nru[1], &spu[spu_offset], p, &bu[bu_offset], &
	    z__[1], &au[au_offset], &q[1], &right[1]);
/*<       go to 180 >*/
    goto L180;
/*<  160  nuu = nu7 >*/
L160:
    nuu = nu7;
/*<       call fptrpe(mu,mv,idim,nu,nru,spu,p,bu,z,au,au1,q,right) >*/
    fptrpe_(mu, mv, idim, nu, &nru[1], &spu[spu_offset], p, &bu[bu_offset], &
	    z__[1], &au[au_offset], &au1[au1_offset], &q[1], &right[1]);
/*  we determine the matrix (avv) and then we reduce this matrix to */
/*  upper triangular form (rv) using givens rotations. */
/*  we apply the same transformations to the columns of matrix */
/*  g to obtain the (nvv) x (nuu) matrix h. */
/*  we store matrix (rv) into av (and av1 if ipar(2)=1) and h into c. */
/*<  180  if(ipar(2).ne.0) go to 190 >*/
L180:
    if (ipar[2] != 0) {
	goto L190;
    }
/*<       nvv = nv4 >*/
    nvv = nv4;
/*<       call fptrnp(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,c,right) >*/
    fptrnp_(mv, &nuu, idim, nv, &nrv[1], &spv[spv_offset], p, &bv[bv_offset], 
	    &q[1], &av[av_offset], &c__[1], &right[1]);
/*<       go to 200 >*/
    goto L200;
/*<  190  nvv = nv7 >*/
L190:
    nvv = nv7;
/*<       call fptrpe(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,av1,c,right) >*/
    fptrpe_(mv, &nuu, idim, nv, &nrv[1], &spv[spv_offset], p, &bv[bv_offset], 
	    &q[1], &av[av_offset], &av1[av1_offset], &c__[1], &right[1]);
/*  backward substitution to obtain the b-spline coefficients as the */
/*  solution of the linear system    (rv) (cr) (ru)' = h. */
/*  first step: solve the system  (rv) (c1) = h. */
/*<  200  ncof = nuu*nvv >*/
L200:
    ncof = nuu * nvv;
/*<       k = 1 >*/
    k = 1;
/*<       if(ipar(2).ne.0) go to 240 >*/
    if (ipar[2] != 0) {
	goto L240;
    }
/*<       do 220 ii=1,idim >*/
    i__1 = *idim;
    for (ii = 1; ii <= i__1; ++ii) {
/*<       do 220 i=1,nuu >*/
	i__2 = nuu;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<          call fpback(av,c(k),nvv,5,c(k),nv) >*/
	    fpback_(&av[av_offset], &c__[k], &nvv, &c__5, &c__[k], nv);
/*<          k = k+nvv >*/
	    k += nvv;
/*<  220  continue >*/
/* L220: */
	}
    }
/*<       go to 300 >*/
    goto L300;
/*<  240  do 260 ii=1,idim >*/
L240:
    i__2 = *idim;
    for (ii = 1; ii <= i__2; ++ii) {
/*<       do 260 i=1,nuu >*/
	i__1 = nuu;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<          call fpbacp(av,av1,c(k),nvv,4,c(k),5,nv) >*/
	    fpbacp_(&av[av_offset], &av1[av1_offset], &c__[k], &nvv, &c__4, &
		    c__[k], &c__5, nv);
/*<          k = k+nvv >*/
	    k += nvv;
/*<  260  continue >*/
/* L260: */
	}
    }
/*  second step: solve the system  (cr) (ru)' = (c1). */
/*<  300  if(ipar(1).ne.0) go to 400 >*/
L300:
    if (ipar[1] != 0) {
	goto L400;
    }
/*<       do 360 ii=1,idim >*/
    i__1 = *idim;
    for (ii = 1; ii <= i__1; ++ii) {
/*<       k = (ii-1)*ncof >*/
	k = (ii - 1) * ncof;
/*<       do 360 j=1,nvv >*/
	i__2 = nvv;
	for (j = 1; j <= i__2; ++j) {
/*<         k = k+1 >*/
	    ++k;
/*<         l = k >*/
	    l = k;
/*<         do 320 i=1,nuu >*/
	    i__3 = nuu;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<           right(i) = c(l) >*/
		right[i__] = c__[l];
/*<           l = l+nvv >*/
		l += nvv;
/*<  320    continue >*/
/* L320: */
	    }
/*<         call fpback(au,right,nuu,5,right,nu) >*/
	    fpback_(&au[au_offset], &right[1], &nuu, &c__5, &right[1], nu);
/*<         l = k >*/
	    l = k;
/*<         do 340 i=1,nuu >*/
	    i__3 = nuu;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<            c(l) = right(i) >*/
		c__[l] = right[i__];
/*<            l = l+nvv >*/
		l += nvv;
/*<  340    continue >*/
/* L340: */
	    }
/*<  360  continue >*/
/* L360: */
	}
    }
/*<       go to 500 >*/
    goto L500;
/*<  400  do 460 ii=1,idim >*/
L400:
    i__2 = *idim;
    for (ii = 1; ii <= i__2; ++ii) {
/*<       k = (ii-1)*ncof >*/
	k = (ii - 1) * ncof;
/*<       do 460 j=1,nvv >*/
	i__1 = nvv;
	for (j = 1; j <= i__1; ++j) {
/*<         k = k+1 >*/
	    ++k;
/*<         l = k >*/
	    l = k;
/*<         do 420 i=1,nuu >*/
	    i__3 = nuu;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<           right(i) = c(l) >*/
		right[i__] = c__[l];
/*<           l = l+nvv >*/
		l += nvv;
/*<  420    continue >*/
/* L420: */
	    }
/*<         call fpbacp(au,au1,right,nuu,4,right,5,nu) >*/
	    fpbacp_(&au[au_offset], &au1[au1_offset], &right[1], &nuu, &c__4, 
		    &right[1], &c__5, nu);
/*<         l = k >*/
	    l = k;
/*<         do 440 i=1,nuu >*/
	    i__3 = nuu;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<            c(l) = right(i) >*/
		c__[l] = right[i__];
/*<            l = l+nvv >*/
		l += nvv;
/*<  440    continue >*/
/* L440: */
	    }
/*<  460  continue >*/
/* L460: */
	}
    }
/*  calculate from the conditions (2)-(3), the remaining b-spline */
/*  coefficients. */
/*<  500  if(ipar(2).eq.0) go to 600 >*/
L500:
    if (ipar[2] == 0) {
	goto L600;
    }
/*<       i = 0 >*/
    i__ = 0;
/*<       j = 0 >*/
    j = 0;
/*<       do 560 id=1,idim >*/
    i__1 = *idim;
    for (id = 1; id <= i__1; ++id) {
/*<       do 560 l=1,nuu >*/
	i__2 = nuu;
	for (l = 1; l <= i__2; ++l) {
/*<          ii = i >*/
	    ii = i__;
/*<          do 520 k=1,nvv >*/
	    i__3 = nvv;
	    for (k = 1; k <= i__3; ++k) {
/*<             i = i+1 >*/
		++i__;
/*<             j = j+1 >*/
		++j;
/*<             q(i) = c(j) >*/
		q[i__] = c__[j];
/*<  520     continue >*/
/* L520: */
	    }
/*<          do 540 k=1,3 >*/
	    for (k = 1; k <= 3; ++k) {
/*<             ii = ii+1 >*/
		++ii;
/*<             i = i+1 >*/
		++i__;
/*<             q(i) = q(ii) >*/
		q[i__] = q[ii];
/*<  540     continue >*/
/* L540: */
	    }
/*<  560  continue >*/
/* L560: */
	}
    }
/*<       ncof = nv4*nuu >*/
    ncof = nv4 * nuu;
/*<       nmd = ncof*idim >*/
    nmd = ncof * *idim;
/*<       do 580 i=1,nmd >*/
    i__2 = nmd;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<          c(i) = q(i) >*/
	c__[i__] = q[i__];
/*<  580  continue >*/
/* L580: */
    }
/*<  600  if(ipar(1).eq.0) go to 700 >*/
L600:
    if (ipar[1] == 0) {
	goto L700;
    }
/*<       i = 0 >*/
    i__ = 0;
/*<       j = 0 >*/
    j = 0;
/*<       n33 = 3*nv4 >*/
    n33 = nv4 * 3;
/*<       do 660 id=1,idim >*/
    i__2 = *idim;
    for (id = 1; id <= i__2; ++id) {
/*<          ii = i >*/
	ii = i__;
/*<          do 620 k=1,ncof >*/
	i__1 = ncof;
	for (k = 1; k <= i__1; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             j = j+1 >*/
	    ++j;
/*<             q(i) = c(j) >*/
	    q[i__] = c__[j];
/*<  620     continue >*/
/* L620: */
	}
/*<          do 640 k=1,n33 >*/
	i__1 = n33;
	for (k = 1; k <= i__1; ++k) {
/*<             ii = ii+1 >*/
	    ++ii;
/*<             i = i+1 >*/
	    ++i__;
/*<             q(i) = q(ii) >*/
	    q[i__] = q[ii];
/*<  640     continue >*/
/* L640: */
	}
/*<  660  continue >*/
/* L660: */
    }
/*<       ncof = nv4*nu4 >*/
    ncof = nv4 * nu4;
/*<       nmd = ncof*idim >*/
    nmd = ncof * *idim;
/*<       do 680 i=1,nmd >*/
    i__2 = nmd;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<          c(i) = q(i) >*/
	c__[i__] = q[i__];
/*<  680  continue >*/
/* L680: */
    }
/*  calculate the quantities */
/*    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv */
/*    fp = sumi=1,mu(sumj=1,mv(res(i,j))) */
/*    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7 */
/*                  tu(r+3) <= u(i) <= tu(r+4) */
/*    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7 */
/*                  tv(r+3) <= v(j) <= tv(r+4) */
/*<  700  fp = 0. >*/
L700:
    *fp = 0.;
/*<       do 720 i=1,nu >*/
    i__2 = *nu;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         fpu(i) = 0. >*/
	fpu[i__] = 0.;
/*<  720  continue >*/
/* L720: */
    }
/*<       do 740 i=1,nv >*/
    i__2 = *nv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/*<         fpv(i) = 0. >*/
	fpv[i__] = 0.;
/*<  740  continue >*/
/* L740: */
    }
/*<       nroldu = 0 >*/
    nroldu = 0;
/*  main loop for the different grid points. */
/*<       do 860 i1=1,muu >*/
    i__2 = muu;
    for (i1 = 1; i1 <= i__2; ++i1) {
/*<         numu = nru(i1) >*/
	numu = nru[i1];
/*<         numu1 = numu+1 >*/
	numu1 = numu + 1;
/*<         nroldv = 0 >*/
	nroldv = 0;
/*<         iz = (i1-1)*mv >*/
	iz = (i1 - 1) * *mv;
/*<         do 840 i2=1,mvv >*/
	i__1 = mvv;
	for (i2 = 1; i2 <= i__1; ++i2) {
/*<           numv = nrv(i2) >*/
	    numv = nrv[i2];
/*<           numv1 = numv+1 >*/
	    numv1 = numv + 1;
/*<           iz = iz+1 >*/
	    ++iz;
/*  evaluate s(u,v) at the current grid point by making the sum of the */
/*  cross products of the non-zero b-splines at (u,v), multiplied with */
/*  the appropiate b-spline coefficients. */
/*<           term = 0. >*/
	    term = 0.;
/*<           k0 = numu*nv4+numv >*/
	    k0 = numu * nv4 + numv;
/*<           jz = iz >*/
	    jz = iz;
/*<           do 800 id=1,idim >*/
	    i__3 = *idim;
	    for (id = 1; id <= i__3; ++id) {
/*<           k1 = k0 >*/
		k1 = k0;
/*<           value = 0. >*/
		value = 0.;
/*<           do 780 l1=1,4 >*/
		for (l1 = 1; l1 <= 4; ++l1) {
/*<             k2 = k1 >*/
		    k2 = k1;
/*<             fac = spu(i1,l1) >*/
		    fac = spu[i1 + l1 * spu_dim1];
/*<             do 760 l2=1,4 >*/
		    for (l2 = 1; l2 <= 4; ++l2) {
/*<               k2 = k2+1 >*/
			++k2;
/*<               value = value+fac*spv(i2,l2)*c(k2) >*/
			value += fac * spv[i2 + l2 * spv_dim1] * c__[k2];
/*<  760        continue >*/
/* L760: */
		    }
/*<             k1 = k1+nv4 >*/
		    k1 += nv4;
/*<  780      continue >*/
/* L780: */
		}
/*  calculate the squared residual at the current grid point. */
/*<           term = term+(z(jz)-value)**2 >*/
/* Computing 2nd power */
		d__1 = z__[jz] - value;
		term += d__1 * d__1;
/*<           jz = jz+mz >*/
		jz += *mz;
/*<           k0 = k0+ncof >*/
		k0 += ncof;
/*<  800      continue >*/
/* L800: */
	    }
/*  adjust the different parameters. */
/*<           fp = fp+term >*/
	    *fp += term;
/*<           fpu(numu1) = fpu(numu1)+term >*/
	    fpu[numu1] += term;
/*<           fpv(numv1) = fpv(numv1)+term >*/
	    fpv[numv1] += term;
/*<           fac = term*half >*/
	    fac = term * half;
/*<           if(numv.eq.nroldv) go to 820 >*/
	    if (numv == nroldv) {
		goto L820;
	    }
/*<           fpv(numv1) = fpv(numv1)-fac >*/
	    fpv[numv1] -= fac;
/*<           fpv(numv) = fpv(numv)+fac >*/
	    fpv[numv] += fac;
/*<  820      nroldv = numv >*/
L820:
	    nroldv = numv;
/*<           if(numu.eq.nroldu) go to 840 >*/
	    if (numu == nroldu) {
		goto L840;
	    }
/*<           fpu(numu1) = fpu(numu1)-fac >*/
	    fpu[numu1] -= fac;
/*<           fpu(numu) = fpu(numu)+fac >*/
	    fpu[numu] += fac;
/*<  840    continue >*/
L840:
	    ;
	}
/*<         nroldu = numu >*/
	nroldu = numu;
/*<  860  continue >*/
/* L860: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpgrpa_ */

