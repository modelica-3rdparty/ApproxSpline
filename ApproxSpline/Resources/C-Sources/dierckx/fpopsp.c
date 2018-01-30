/* fpopsp.f -- translated by f2c (version 20061008).
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

static integer c__0 = 0;
static integer c__1 = 1;

/*<    >*/
/* Subroutine */ int fpopsp_(integer *ifsu, integer *ifsv, integer *ifbu, 
	integer *ifbv, doublereal *u, integer *mu, doublereal *v, integer *mv,
	 doublereal *r__, integer *mr, doublereal *r0, doublereal *r1, 
	doublereal *dr, integer *iopt, integer *ider, doublereal *tu, integer 
	*nu, doublereal *tv, integer *nv, integer *nuest, integer *nvest, 
	doublereal *p, doublereal *step, doublereal *c__, integer *nc, 
	doublereal *fp, doublereal *fpu, doublereal *fpv, integer *nru, 
	integer *nrv, doublereal *wrk, integer *lwrk)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal a[36]	/* was [6][6] */, g[6];
    static integer i__, j, l, i1, l1, l2, mm, lq, nr[6];
    static doublereal sq;
    static integer id0, la0, la1, lb0, lb1, lc0, lc1, id1;
    static doublereal sq0, sq1;
    static integer lau, lbu, lbv, lcs, lri;
    static doublereal drr[6], sqq;
    static integer lsu, lsv;
    static doublereal sum[6];
    static integer lav1, lav2, iop0, iop1, mvnu;
    static doublereal step1, step2, delta[6], three;
    static integer number;
    extern /* Subroutine */ int fpgrsp_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    fpsysy_(doublereal *, integer *, doublereal *);

/*  given the set of function values r(i,j) defined on the rectangular */
/*  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopsp determines a */
/*  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu */
/*  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this */
/*  spline sp(u,v) will be periodic in the variable v and will satisfy */
/*  the following constraints */

/*     s(tu(1),v) = dr(1) , tv(4) <=v<= tv(nv-3) */

/*     s(tu(nu),v) = dr(4) , tv(4) <=v<= tv(nv-3) */

/*  and (if iopt(2) = 1) */

/*     d s(tu(1),v) */
/*     ------------ =  dr(2)*cos(v)+dr(3)*sin(v) , tv(4) <=v<= tv(nv-3) */
/*     d u */

/*  and (if iopt(3) = 1) */

/*     d s(tu(nu),v) */
/*     ------------- =  dr(5)*cos(v)+dr(6)*sin(v) , tv(4) <=v<= tv(nv-3) */
/*     d u */

/*  where the parameters dr(i) correspond to the derivative values at the */
/*  poles as defined in subroutine spgrid. */

/*  the b-spline coefficients of sp(u,v) are determined as the least- */
/*  squares solution  of an overdetermined linear system which depends */
/*  on the value of p and on the values dr(i),i=1,...,6. the correspond- */
/*  ing sum of squared residuals sq is a simple quadratic function in */
/*  the variables dr(i). these may or may not be provided. the values */
/*  dr(i) which are not given will be determined so as to minimize the */
/*  resulting sum of squared residuals sq. in that case the user must */
/*  provide some initial guess dr(i) and some estimate (dr(i)-step, */
/*  dr(i)+step) of the range of possible values for these latter. */

/*  sp(u,v) also depends on the parameter p (p>0) in such a way that */
/*    - if p tends to infinity, sp(u,v) becomes the least-squares spline */
/*      with given knots, satisfying the constraints. */
/*    - if p tends to zero, sp(u,v) becomes the least-squares polynomial, */
/*      satisfying the constraints. */
/*    - the function  f(p)=sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2) */
/*      is continuous and strictly decreasing for p>0. */

/*  ..scalar arguments.. */
/*<    >*/
/*<       real r0,r1,p,fp >*/
/*  ..array arguments.. */
/*<       integer ider(4),nru(mu),nrv(mv),iopt(3) >*/
/*<    >*/
/*  ..local scalars.. */
/*<       real res,sq,sqq,sq0,sq1,step1,step2,three >*/
/*<    >*/
/*  ..local arrays.. */
/*<       integer nr(6) >*/
/*<       real delta(6),drr(6),sum(6),a(6,6),g(6) >*/
/*  ..function references.. */
/*<       integer max0 >*/
/*  ..subroutine references.. */
/*    fpgrsp,fpsysy */
/*  .. */
/*  set constant */
/*<       three = 3 >*/
    /* Parameter adjustments */
    --nru;
    --u;
    --nrv;
    --v;
    --r__;
    --dr;
    --iopt;
    --ider;
    --fpu;
    --tu;
    --fpv;
    --tv;
    --step;
    --c__;
    --wrk;

    /* Function Body */
    three = 3.;
/*  we partition the working space */
/*<       lsu = 1 >*/
    lsu = 1;
/*<       lsv = lsu+4*mu >*/
    lsv = lsu + (*mu << 2);
/*<       lri = lsv+4*mv >*/
    lri = lsv + (*mv << 2);
/*<       mm = max0(nuest,mv+nvest) >*/
/* Computing MAX */
    i__1 = *nuest, i__2 = *mv + *nvest;
    mm = max(i__1,i__2);
/*<       lq = lri+mm >*/
    lq = lri + mm;
/*<       mvnu = nuest*(mv+nvest-8) >*/
    mvnu = *nuest * (*mv + *nvest - 8);
/*<       lau = lq+mvnu >*/
    lau = lq + mvnu;
/*<       lav1 = lau+5*nuest >*/
    lav1 = lau + *nuest * 5;
/*<       lav2 = lav1+6*nvest >*/
    lav2 = lav1 + *nvest * 6;
/*<       lbu = lav2+4*nvest >*/
    lbu = lav2 + (*nvest << 2);
/*<       lbv = lbu+5*nuest >*/
    lbv = lbu + *nuest * 5;
/*<       la0 = lbv+5*nvest >*/
    la0 = lbv + *nvest * 5;
/*<       la1 = la0+2*mv >*/
    la1 = la0 + (*mv << 1);
/*<       lb0 = la1+2*mv >*/
    lb0 = la1 + (*mv << 1);
/*<       lb1 = lb0+2*nvest >*/
    lb1 = lb0 + (*nvest << 1);
/*<       lc0 = lb1+2*nvest >*/
    lc0 = lb1 + (*nvest << 1);
/*<       lc1 = lc0+nvest >*/
    lc1 = lc0 + *nvest;
/*<       lcs = lc1+nvest >*/
    lcs = lc1 + *nvest;
/*  we calculate the smoothing spline sp(u,v) according to the input */
/*  values dr(i),i=1,...,6. */
/*<       iop0 = iopt(2) >*/
    iop0 = iopt[2];
/*<       iop1 = iopt(3) >*/
    iop1 = iopt[3];
/*<       id0 = ider(1) >*/
    id0 = ider[1];
/*<       id1 = ider(3) >*/
    id1 = ider[3];
/*<    >*/
    fpgrsp_(ifsu, ifsv, ifbu, ifbv, &c__0, &u[1], mu, &v[1], mv, &r__[1], mr, 
	    &dr[1], &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc, &sq,
	     fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv], &wrk[lri]
	    , &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[lav2], &wrk[lbu], &wrk[
	    lbv], &wrk[la0], &wrk[la1], &wrk[lb0], &wrk[lb1], &wrk[lc0], &wrk[
	    lc1], &wrk[lcs], &nru[1], &nrv[1]);
/*<       sq0 = 0. >*/
    sq0 = 0.;
/*<       sq1 = 0. >*/
    sq1 = 0.;
/*<       if(id0.eq.0) sq0 = (r0-dr(1))**2 >*/
    if (id0 == 0) {
/* Computing 2nd power */
	d__1 = *r0 - dr[1];
	sq0 = d__1 * d__1;
    }
/*<       if(id1.eq.0) sq1 = (r1-dr(4))**2 >*/
    if (id1 == 0) {
/* Computing 2nd power */
	d__1 = *r1 - dr[4];
	sq1 = d__1 * d__1;
    }
/*<       sq = sq+sq0+sq1 >*/
    sq = sq + sq0 + sq1;
/* in case all derivative values dr(i) are given (step<=0) or in case */
/* we have spline interpolation, we accept this spline as a solution. */
/*<       if(sq.le.0.) return >*/
    if (sq <= 0.) {
	return 0;
    }
/*<       if(step(1).le.0. .and. step(2).le.0.) return >*/
    if (step[1] <= 0. && step[2] <= 0.) {
	return 0;
    }
/*<       do 10 i=1,6 >*/
    for (i__ = 1; i__ <= 6; ++i__) {
/*<         drr(i) = dr(i) >*/
	drr[i__ - 1] = dr[i__];
/*<   10  continue >*/
/* L10: */
    }
/* number denotes the number of derivative values dr(i) that still must */
/* be optimized. let us denote these parameters by g(j),j=1,...,number. */
/*<       number = 0 >*/
    number = 0;
/*<       if(id0.gt.0) go to 20 >*/
    if (id0 > 0) {
	goto L20;
    }
/*<       number = 1 >*/
    number = 1;
/*<       nr(1) = 1 >*/
    nr[0] = 1;
/*<       delta(1) = step(1) >*/
    delta[0] = step[1];
/*<   20  if(iop0.eq.0) go to 30 >*/
L20:
    if (iop0 == 0) {
	goto L30;
    }
/*<       if(ider(2).ne.0) go to 30 >*/
    if (ider[2] != 0) {
	goto L30;
    }
/*<       step2 = step(1)*three/(tu(5)-tu(4)) >*/
    step2 = step[1] * three / (tu[5] - tu[4]);
/*<       nr(number+1) = 2 >*/
    nr[number] = 2;
/*<       nr(number+2) = 3 >*/
    nr[number + 1] = 3;
/*<       delta(number+1) = step2 >*/
    delta[number] = step2;
/*<       delta(number+2) = step2 >*/
    delta[number + 1] = step2;
/*<       number = number+2 >*/
    number += 2;
/*<   30  if(id1.gt.0) go to 40 >*/
L30:
    if (id1 > 0) {
	goto L40;
    }
/*<       number = number+1 >*/
    ++number;
/*<       nr(number) = 4 >*/
    nr[number - 1] = 4;
/*<       delta(number) = step(2) >*/
    delta[number - 1] = step[2];
/*<   40  if(iop1.eq.0) go to 50 >*/
L40:
    if (iop1 == 0) {
	goto L50;
    }
/*<       if(ider(4).ne.0) go to 50 >*/
    if (ider[4] != 0) {
	goto L50;
    }
/*<       step2 = step(2)*three/(tu(nu)-tu(nu-4)) >*/
    step2 = step[2] * three / (tu[*nu] - tu[*nu - 4]);
/*<       nr(number+1) = 5 >*/
    nr[number] = 5;
/*<       nr(number+2) = 6 >*/
    nr[number + 1] = 6;
/*<       delta(number+1) = step2 >*/
    delta[number] = step2;
/*<       delta(number+2) = step2 >*/
    delta[number + 1] = step2;
/*<       number = number+2 >*/
    number += 2;
/*<   50  if(number.eq.0) return >*/
L50:
    if (number == 0) {
	return 0;
    }
/* the sum of squared residulas sq is a quadratic polynomial in the */
/* parameters g(j). we determine the unknown coefficients of this */
/* polymomial by calculating (number+1)*(number+2)/2 different splines */
/* according to specific values for g(j). */
/*<       do 60 i=1,number >*/
    i__1 = number;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          l = nr(i) >*/
	l = nr[i__ - 1];
/*<          step1 = delta(i) >*/
	step1 = delta[i__ - 1];
/*<          drr(l) = dr(l)+step1 >*/
	drr[l - 1] = dr[l] + step1;
/*<    >*/
	fpgrsp_(ifsu, ifsv, ifbu, ifbv, &c__1, &u[1], mu, &v[1], mv, &r__[1], 
		mr, drr, &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc,
		 &sum[i__ - 1], fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &
		wrk[lsv], &wrk[lri], &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[
		lav2], &wrk[lbu], &wrk[lbv], &wrk[la0], &wrk[la1], &wrk[lb0], 
		&wrk[lb1], &wrk[lc0], &wrk[lc1], &wrk[lcs], &nru[1], &nrv[1]);
/*<          if(id0.eq.0) sq0 = (r0-drr(1))**2 >*/
	if (id0 == 0) {
/* Computing 2nd power */
	    d__1 = *r0 - drr[0];
	    sq0 = d__1 * d__1;
	}
/*<          if(id1.eq.0) sq1 = (r1-drr(4))**2 >*/
	if (id1 == 0) {
/* Computing 2nd power */
	    d__1 = *r1 - drr[3];
	    sq1 = d__1 * d__1;
	}
/*<          sum(i) = sum(i)+sq0+sq1 >*/
	sum[i__ - 1] = sum[i__ - 1] + sq0 + sq1;
/*<          drr(l) = dr(l)-step1 >*/
	drr[l - 1] = dr[l] - step1;
/*<    >*/
	fpgrsp_(ifsu, ifsv, ifbu, ifbv, &c__1, &u[1], mu, &v[1], mv, &r__[1], 
		mr, drr, &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc,
		 &sqq, fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv],
		 &wrk[lri], &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[lav2], &wrk[
		lbu], &wrk[lbv], &wrk[la0], &wrk[la1], &wrk[lb0], &wrk[lb1], &
		wrk[lc0], &wrk[lc1], &wrk[lcs], &nru[1], &nrv[1]);
/*<          if(id0.eq.0) sq0 = (r0-drr(1))**2 >*/
	if (id0 == 0) {
/* Computing 2nd power */
	    d__1 = *r0 - drr[0];
	    sq0 = d__1 * d__1;
	}
/*<          if(id1.eq.0) sq1 = (r1-drr(4))**2 >*/
	if (id1 == 0) {
/* Computing 2nd power */
	    d__1 = *r1 - drr[3];
	    sq1 = d__1 * d__1;
	}
/*<          sqq = sqq+sq0+sq1 >*/
	sqq = sqq + sq0 + sq1;
/*<          drr(l) = dr(l) >*/
	drr[l - 1] = dr[l];
/*<          a(i,i) = (sum(i)+sqq-sq-sq)/step1**2 >*/
/* Computing 2nd power */
	d__1 = step1;
	a[i__ + i__ * 6 - 7] = (sum[i__ - 1] + sqq - sq - sq) / (d__1 * d__1);
/*<          if(a(i,i).le.0.) go to 110 >*/
	if (a[i__ + i__ * 6 - 7] <= 0.) {
	    goto L110;
	}
/*<          g(i) = (sqq-sum(i))/(step1+step1) >*/
	g[i__ - 1] = (sqq - sum[i__ - 1]) / (step1 + step1);
/*<   60  continue >*/
/* L60: */
    }
/*<       if(number.eq.1) go to 90 >*/
    if (number == 1) {
	goto L90;
    }
/*<       do 80 i=2,number >*/
    i__1 = number;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          l1 = nr(i) >*/
	l1 = nr[i__ - 1];
/*<          step1 = delta(i) >*/
	step1 = delta[i__ - 1];
/*<          drr(l1) = dr(l1)+step1 >*/
	drr[l1 - 1] = dr[l1] + step1;
/*<          i1 = i-1 >*/
	i1 = i__ - 1;
/*<          do 70 j=1,i1 >*/
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
/*<             l2 = nr(j) >*/
	    l2 = nr[j - 1];
/*<             step2 = delta(j) >*/
	    step2 = delta[j - 1];
/*<             drr(l2) = dr(l2)+step2 >*/
	    drr[l2 - 1] = dr[l2] + step2;
/*<    >*/
	    fpgrsp_(ifsu, ifsv, ifbu, ifbv, &c__1, &u[1], mu, &v[1], mv, &r__[
		    1], mr, drr, &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &
		    c__[1], nc, &sqq, fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[
		    lsu], &wrk[lsv], &wrk[lri], &wrk[lq], &wrk[lau], &wrk[
		    lav1], &wrk[lav2], &wrk[lbu], &wrk[lbv], &wrk[la0], &wrk[
		    la1], &wrk[lb0], &wrk[lb1], &wrk[lc0], &wrk[lc1], &wrk[
		    lcs], &nru[1], &nrv[1]);
/*<             if(id0.eq.0) sq0 = (r0-drr(1))**2 >*/
	    if (id0 == 0) {
/* Computing 2nd power */
		d__1 = *r0 - drr[0];
		sq0 = d__1 * d__1;
	    }
/*<             if(id1.eq.0) sq1 = (r1-drr(4))**2 >*/
	    if (id1 == 0) {
/* Computing 2nd power */
		d__1 = *r1 - drr[3];
		sq1 = d__1 * d__1;
	    }
/*<             sqq = sqq+sq0+sq1 >*/
	    sqq = sqq + sq0 + sq1;
/*<             a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2) >*/
	    a[i__ + j * 6 - 7] = (sq + sqq - sum[i__ - 1] - sum[j - 1]) / (
		    step1 * step2);
/*<             drr(l2) = dr(l2) >*/
	    drr[l2 - 1] = dr[l2];
/*<   70     continue >*/
/* L70: */
	}
/*<          drr(l1) = dr(l1) >*/
	drr[l1 - 1] = dr[l1];
/*<   80  continue >*/
/* L80: */
    }
/* the optimal values g(j) are found as the solution of the system */
/* d (sq) / d (g(j)) = 0 , j=1,...,number. */
/*<   90  call fpsysy(a,number,g) >*/
L90:
    fpsysy_(a, &number, g);
/*<       do 100 i=1,number >*/
    i__1 = number;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          l = nr(i) >*/
	l = nr[i__ - 1];
/*<          dr(l) = dr(l)+g(i) >*/
	dr[l] += g[i__ - 1];
/*<  100  continue >*/
/* L100: */
    }
/* we determine the spline sp(u,v) according to the optimal values g(j). */
/*<  1 >*/
L110:
    fpgrsp_(ifsu, ifsv, ifbu, ifbv, &c__0, &u[1], mu, &v[1], mv, &r__[1], mr, 
	    &dr[1], &iop0, &iop1, &tu[1], nu, &tv[1], nv, p, &c__[1], nc, &sq,
	     fp, &fpu[1], &fpv[1], &mm, &mvnu, &wrk[lsu], &wrk[lsv], &wrk[lri]
	    , &wrk[lq], &wrk[lau], &wrk[lav1], &wrk[lav2], &wrk[lbu], &wrk[
	    lbv], &wrk[la0], &wrk[la1], &wrk[lb0], &wrk[lb1], &wrk[lc0], &wrk[
	    lc1], &wrk[lcs], &nru[1], &nrv[1]);
/*<       if(id0.eq.0) sq0 = (r0-dr(1))**2 >*/
    if (id0 == 0) {
/* Computing 2nd power */
	d__1 = *r0 - dr[1];
	sq0 = d__1 * d__1;
    }
/*<       if(id1.eq.0) sq1 = (r1-dr(4))**2 >*/
    if (id1 == 0) {
/* Computing 2nd power */
	d__1 = *r1 - dr[4];
	sq1 = d__1 * d__1;
    }
/*<       sq = sq+sq0+sq1 >*/
    sq = sq + sq0 + sq1;
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpopsp_ */

