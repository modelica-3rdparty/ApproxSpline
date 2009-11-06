/* fppola.f -- translated by f2c (version 20061008).
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

/*<    >*/
/* Subroutine */ int fppola_(integer *iopt1, integer *iopt2, integer *iopt3, 
	integer *m, doublereal *u, doublereal *v, doublereal *z__, doublereal 
	*w, D_fp rad, doublereal *s, integer *nuest, integer *nvest, 
	doublereal *eta, doublereal *tol, integer *maxit, integer *ib1, 
	integer *ib3, integer *nc, integer *ncc, integer *intest, integer *
	nrest, integer *nu, doublereal *tu, integer *nv, doublereal *tv, 
	doublereal *c__, doublereal *fp, doublereal *sup, doublereal *fpint, 
	doublereal *coord, doublereal *f, doublereal *ff, doublereal *row, 
	doublereal *cs, doublereal *cosi, doublereal *a, doublereal *q, 
	doublereal *bu, doublereal *bv, doublereal *spu, doublereal *spv, 
	doublereal *h__, integer *index, integer *nummer, doublereal *wrk, 
	integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, bu_dim1, bu_offset, bv_dim1, 
	    bv_offset, spu_dim1, spu_offset, spv_dim1, spv_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal p, r__, c1, c2, c3, c4, f1, f2, f3;
    static integer i1, i2, i3, j1, j2, l1;
    static doublereal p1, p2, p3;
    static integer l2, l3, l4;
    static doublereal u2, u3;
    static integer la;
    static doublereal co;
    static integer ii, lf, il;
    static doublereal pi;
    static integer in;
    static doublereal si;
    static integer lh, ll;
    static doublereal wi, rn;
    static integer lu;
    static doublereal sq, zi;
    static integer lv;
    static doublereal hu[4], uu, hv[4], pi2;
    static integer nr1, nu4, nv4;
    static doublereal acc, fac, arg, one, hui, huj, eps, ten;
    static integer jlu;
    static doublereal piv;
    static integer num, nrr;
    static doublereal fac1, fac2, two;
    static integer nvv, nuu, ich1, ich3;
    static doublereal con1, con4, con9;
    static integer num1;
    static doublereal half;
    static integer ncof;
    static doublereal dmax__;
    static integer ipar, nreg, rank, iter;
    static doublereal fpms, pinv;
    static integer irot, jrot, ipar1, iband, ncoff;
    static doublereal sigma, three, fpmax, ratio;
    static integer numin, nvmin, nrint;
    static doublereal store;
    static integer iband3, iband4, lwest, iband1;
    extern /* Subroutine */ int fpback_(doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *), fpdisc_(doublereal *, 
	    integer *, integer *, doublereal *, integer *), fporde_(
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *), fprank_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *)
	    ;
    extern doublereal fprati_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *), fprota_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fpgivs_(doublereal *, 
	    doublereal *, doublereal *, doublereal *), fprppo_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);

/*  ..scalar arguments.. */
/*<    >*/
/*<       real s,eta,tol,fp,sup >*/
/*  ..array arguments.. */
/*<       integer index(nrest),nummer(m) >*/
/*<    >*/
/*  ..user supplied function.. */
/*<       real rad >*/
/*  ..local scalars.. */
/*<    >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real hu(4),hv(4) >*/
/*  ..function references.. */
/*<       real abs,atan,cos,fprati,sin,sqrt >*/
/*<       integer min0 >*/
/*  ..subroutine references.. */
/*    fporde,fpbspl,fpback,fpgivs,fprota,fprank,fpdisc,fprppo */
/*  .. */
/*  set constants */
/*<       one = 1 >*/
    /* Parameter adjustments */
    --nummer;
    spv_dim1 = *m;
    spv_offset = 1 + spv_dim1;
    spv -= spv_offset;
    spu_dim1 = *m;
    spu_offset = 1 + spu_dim1;
    spu -= spu_offset;
    --w;
    --z__;
    --v;
    --u;
    bu_dim1 = *nuest;
    bu_offset = 1 + bu_dim1;
    bu -= bu_offset;
    --tu;
    bv_dim1 = *nvest;
    bv_offset = 1 + bv_dim1;
    bv -= bv_offset;
    cosi -= 6;
    --cs;
    --row;
    --tv;
    --h__;
    --ff;
    --c__;
    q_dim1 = *ncc;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    a_dim1 = *ncc;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --f;
    --coord;
    --fpint;
    --index;
    --wrk;

    /* Function Body */
    one = 1.;
/*<       two = 2 >*/
    two = 2.;
/*<       three = 3 >*/
    three = 3.;
/*<       ten = 10 >*/
    ten = 10.;
/*<       half = 0.5e0 >*/
    half = .5;
/*<       con1 = 0.1e0 >*/
    con1 = .1;
/*<       con9 = 0.9e0 >*/
    con9 = .9;
/*<       con4 = 0.4e-01 >*/
    con4 = .04;
/*<       pi = atan(one)*4 >*/
    pi = atan(one) * 4;
/*<       pi2 = pi+pi >*/
    pi2 = pi + pi;
/*<       ipar = iopt2*(iopt2+3)/2 >*/
    ipar = *iopt2 * (*iopt2 + 3) / 2;
/*<       ipar1 = ipar+1 >*/
    ipar1 = ipar + 1;
/*<       eps = sqrt(eta) >*/
    eps = sqrt(*eta);
/*<       if(iopt1.lt.0) go to 90 >*/
    if (*iopt1 < 0) {
	goto L90;
    }
/*<       numin = 9 >*/
    numin = 9;
/*<       nvmin = 9+iopt2*(iopt2+1) >*/
    nvmin = *iopt2 * (*iopt2 + 1) + 9;
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
/*<       acc = tol*s >*/
    acc = *tol * *s;
/*<       if(iopt1.eq.0) go to 10 >*/
    if (*iopt1 == 0) {
	goto L10;
    }
/*<       if(s.lt.sup) if(nv-nvmin) 70,90,90 >*/
    if (*s < *sup) {
	if (*nv - nvmin >= 0) {
	    goto L90;
	} else {
	    goto L70;
	}
    }
/*  if iopt1 = 0 we begin by computing the weighted least-squares */
/*  polymomial of the form */
/*     s(u,v) = f(1)*(1-u**3)+f(2)*u**3+f(3)*(u**2-u**3)+f(4)*(u-u**3) */
/*  where f(4) = 0 if iopt2> 0 , f(3) = 0 if iopt2 > 1 and */
/*        f(2) = 0 if iopt3> 0. */
/*  the corresponding weighted sum of squared residuals gives the upper */
/*  bound sup for the smoothing factor s. */
/*<   10  sup = 0. >*/
L10:
    *sup = 0.;
/*<       do 20 i=1,4 >*/
    for (i__ = 1; i__ <= 4; ++i__) {
/*<          f(i) = 0. >*/
	f[i__] = 0.;
/*<          do 20 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             a(i,j) = 0. >*/
	    a[i__ + j * a_dim1] = 0.;
/*<  20   continue >*/
/* L20: */
	}
    }
/*<       do 50 i=1,m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          wi = w(i) >*/
	wi = w[i__];
/*<          zi = z(i)*wi >*/
	zi = z__[i__] * wi;
/*<          uu = u(i) >*/
	uu = u[i__];
/*<          u2 = uu*uu >*/
	u2 = uu * uu;
/*<          u3 = uu*u2 >*/
	u3 = uu * u2;
/*<          h(1) = (one-u3)*wi >*/
	h__[1] = (one - u3) * wi;
/*<          h(2) = u3*wi >*/
	h__[2] = u3 * wi;
/*<          h(3) = u2*(one-uu)*wi >*/
	h__[3] = u2 * (one - uu) * wi;
/*<          h(4) = uu*(one-u2)*wi >*/
	h__[4] = uu * (one - u2) * wi;
/*<          if(iopt3.ne.0) h(2) = 0. >*/
	if (*iopt3 != 0) {
	    h__[2] = 0.;
	}
/*<          if(iopt2.gt.1) h(3) = 0. >*/
	if (*iopt2 > 1) {
	    h__[3] = 0.;
	}
/*<          if(iopt2.gt.0) h(4) = 0. >*/
	if (*iopt2 > 0) {
	    h__[4] = 0.;
	}
/*<          do 40 j=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<             piv = h(j) >*/
	    piv = h__[j];
/*<             if(piv.eq.0.) go to 40 >*/
	    if (piv == 0.) {
		goto L40;
	    }
/*<             call fpgivs(piv,a(j,1),co,si) >*/
	    fpgivs_(&piv, &a[j + a_dim1], &co, &si);
/*<             call fprota(co,si,zi,f(j)) >*/
	    fprota_(&co, &si, &zi, &f[j]);
/*<             if(j.eq.4) go to 40 >*/
	    if (j == 4) {
		goto L40;
	    }
/*<             j1 = j+1 >*/
	    j1 = j + 1;
/*<             j2 = 1 >*/
	    j2 = 1;
/*<             do 30 l=j1,4 >*/
	    for (l = j1; l <= 4; ++l) {
/*<                j2 = j2+1 >*/
		++j2;
/*<                call fprota(co,si,h(l),a(j,j2)) >*/
		fprota_(&co, &si, &h__[l], &a[j + j2 * a_dim1]);
/*<   30        continue >*/
/* L30: */
	    }
/*<   40     continue >*/
L40:
	    ;
	}
/*<          sup = sup+zi*zi >*/
	*sup += zi * zi;
/*<   50  continue >*/
/* L50: */
    }
/*<       if(a(4,1).ne.0.) f(4) = f(4)/a(4,1) >*/
    if (a[a_dim1 + 4] != 0.) {
	f[4] /= a[a_dim1 + 4];
    }
/*<       if(a(3,1).ne.0.) f(3) = (f(3)-a(3,2)*f(4))/a(3,1) >*/
    if (a[a_dim1 + 3] != 0.) {
	f[3] = (f[3] - a[(a_dim1 << 1) + 3] * f[4]) / a[a_dim1 + 3];
    }
/*<       if(a(2,1).ne.0.) f(2) = (f(2)-a(2,2)*f(3)-a(2,3)*f(4))/a(2,1) >*/
    if (a[a_dim1 + 2] != 0.) {
	f[2] = (f[2] - a[(a_dim1 << 1) + 2] * f[3] - a[a_dim1 * 3 + 2] * f[4])
		 / a[a_dim1 + 2];
    }
/*<    >*/
    if (a[a_dim1 + 1] != 0.) {
	f[1] = (f[1] - a[(a_dim1 << 1) + 1] * f[2] - a[a_dim1 * 3 + 1] * f[3] 
		- a[(a_dim1 << 2) + 1] * f[4]) / a[a_dim1 + 1];
    }
/*  find the b-spline representation of this least-squares polynomial */
/*<       c1 = f(1) >*/
    c1 = f[1];
/*<       c4 = f(2) >*/
    c4 = f[2];
/*<       c2 = f(4)/three+c1 >*/
    c2 = f[4] / three + c1;
/*<       c3 = (f(3)+two*f(4))/three+c1 >*/
    c3 = (f[3] + two * f[4]) / three + c1;
/*<       nu = 8 >*/
    *nu = 8;
/*<       nv = 8 >*/
    *nv = 8;
/*<       do 60 i=1,4 >*/
    for (i__ = 1; i__ <= 4; ++i__) {
/*<          c(i) = c1 >*/
	c__[i__] = c1;
/*<          c(i+4) = c2 >*/
	c__[i__ + 4] = c2;
/*<          c(i+8) = c3 >*/
	c__[i__ + 8] = c3;
/*<          c(i+12) = c4 >*/
	c__[i__ + 12] = c4;
/*<          tu(i) = 0. >*/
	tu[i__] = 0.;
/*<          tu(i+4) = one >*/
	tu[i__ + 4] = one;
/*<          rn = 2*i-9 >*/
	rn = (doublereal) ((i__ << 1) - 9);
/*<          tv(i) = rn*pi >*/
	tv[i__] = rn * pi;
/*<          rn = 2*i-1 >*/
	rn = (doublereal) ((i__ << 1) - 1);
/*<          tv(i+4) = rn*pi >*/
	tv[i__ + 4] = rn * pi;
/*<   60  continue >*/
/* L60: */
    }
/*<       fp = sup >*/
    *fp = *sup;
/*  test whether the least-squares polynomial is an acceptable solution */
/*<       fpms = sup-s >*/
    fpms = *sup - *s;
/*<       if(fpms.lt.acc) go to 960 >*/
    if (fpms < acc) {
	goto L960;
    }
/*  test whether we cannot further increase the number of knots. */
/*<   70  if(nuest.lt.numin .or. nvest.lt.nvmin) go to 950 >*/
L70:
    if (*nuest < numin || *nvest < nvmin) {
	goto L950;
    }
/*  find the initial set of interior knots of the spline in case iopt1=0. */
/*<       nu = numin >*/
    *nu = numin;
/*<       nv = nvmin >*/
    *nv = nvmin;
/*<       tu(5) = half >*/
    tu[5] = half;
/*<       nvv = nv-8 >*/
    nvv = *nv - 8;
/*<       rn = nvv+1 >*/
    rn = (doublereal) (nvv + 1);
/*<       fac = pi2/rn >*/
    fac = pi2 / rn;
/*<       do 80 i=1,nvv >*/
    i__1 = nvv;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          rn = i >*/
	rn = (doublereal) i__;
/*<          tv(i+4) = rn*fac-pi >*/
	tv[i__ + 4] = rn * fac - pi;
/*<   80  continue >*/
/* L80: */
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 1 : computation of least-squares bicubic splines.              c */
/*  ******************************************************              c */
/*  if iopt1<0 we compute the least-squares bicubic spline according    c */
/*  to the given set of knots.                                          c */
/*  if iopt1>=0 we compute least-squares bicubic splines with in-       c */
/*  creasing numbers of knots until the corresponding sum f(p=inf)<=s.  c */
/*  the initial set of knots then depends on the value of iopt1         c */
/*    if iopt1=0 we start with one interior knot in the u-direction     c */
/*              (0.5) and 1+iopt2*(iopt2+1) in the v-direction.         c */
/*    if iopt1>0 we start with the set of knots found at the last       c */
/*              call of the routine.                                    c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
/*<   90  do 570 iter=1,m >*/
L90:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(u,v). */
/*<          l1 = 4 >*/
	l1 = 4;
/*<          l2 = l1 >*/
	l2 = l1;
/*<          l3 = nv-3 >*/
	l3 = *nv - 3;
/*<          l4 = l3 >*/
	l4 = l3;
/*<          tv(l2) = -pi >*/
	tv[l2] = -pi;
/*<          tv(l3) = pi >*/
	tv[l3] = pi;
/*<          do 120 i=1,3 >*/
	for (i__ = 1; i__ <= 3; ++i__) {
/*<             l1 = l1+1 >*/
	    ++l1;
/*<             l2 = l2-1 >*/
	    --l2;
/*<             l3 = l3+1 >*/
	    ++l3;
/*<             l4 = l4-1 >*/
	    --l4;
/*<             tv(l2) = tv(l4)-pi2 >*/
	    tv[l2] = tv[l4] - pi2;
/*<             tv(l3) = tv(l1)+pi2 >*/
	    tv[l3] = tv[l1] + pi2;
/*<  120     continue >*/
/* L120: */
	}
/*<         l = nu >*/
	l = *nu;
/*<         do 130 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           tu(i) = 0. >*/
	    tu[i__] = 0.;
/*<           tu(l) = one >*/
	    tu[l] = one;
/*<           l = l-1 >*/
	    --l;
/*<  130    continue >*/
/* L130: */
	}
/*  find nrint, the total number of knot intervals and nreg, the number */
/*  of panels in which the approximation domain is subdivided by the */
/*  intersection of knots. */
/*<         nuu = nu-7 >*/
	nuu = *nu - 7;
/*<         nvv = nv-7 >*/
	nvv = *nv - 7;
/*<         nrr = nvv/2 >*/
	nrr = nvv / 2;
/*<         nr1 = nrr+1 >*/
	nr1 = nrr + 1;
/*<         nrint = nuu+nvv >*/
	nrint = nuu + nvv;
/*<         nreg = nuu*nvv >*/
	nreg = nuu * nvv;
/*  arrange the data points according to the panel they belong to. */
/*<         call fporde(u,v,m,3,3,tu,nu,tv,nv,nummer,index,nreg) >*/
	fporde_(&u[1], &v[1], m, &c__3, &c__3, &tu[1], nu, &tv[1], nv, &
		nummer[1], &index[1], &nreg);
/*<         if(iopt2.eq.0) go to 195 >*/
	if (*iopt2 == 0) {
	    goto L195;
	}
/*  find the b-spline coefficients cosi of the cubic spline */
/*  approximations for cr(v)=rad(v)*cos(v) and sr(v) = rad(v)*sin(v) */
/*  if iopt2=1, and additionally also for cr(v)**2,sr(v)**2 and */
/*  2*cr(v)*sr(v) if iopt2=2 */
/*<         do 140 i=1,nvv >*/
	i__2 = nvv;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<            do 135 j=1,ipar >*/
	    i__3 = ipar;
	    for (j = 1; j <= i__3; ++j) {
/*<               cosi(j,i) = 0. >*/
		cosi[j + i__ * 5] = 0.;
/*<  135       continue >*/
/* L135: */
	    }
/*<            do 140 j=1,nvv >*/
	    i__3 = nvv;
	    for (j = 1; j <= i__3; ++j) {
/*<               a(i,j) = 0. >*/
		a[i__ + j * a_dim1] = 0.;
/*<  140    continue >*/
/* L140: */
	    }
	}
/*  the coefficients cosi are obtained from interpolation conditions */
/*  at the knots tv(i),i=4,5,...nv-4. */
/*<         do 175 i=1,nvv >*/
	i__3 = nvv;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<            l2 = i+3 >*/
	    l2 = i__ + 3;
/*<            arg = tv(l2) >*/
	    arg = tv[l2];
/*<            call fpbspl(tv,nv,3,arg,l2,hv) >*/
	    fpbspl_(&tv[1], nv, &c__3, &arg, &l2, hv);
/*<            do 145 j=1,nvv >*/
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
/*<               row(j) = 0. >*/
		row[j] = 0.;
/*<  145       continue >*/
/* L145: */
	    }
/*<            ll = i >*/
	    ll = i__;
/*<            do 150 j=1,3 >*/
	    for (j = 1; j <= 3; ++j) {
/*<               if(ll.gt.nvv) ll= 1 >*/
		if (ll > nvv) {
		    ll = 1;
		}
/*<               row(ll) = row(ll)+hv(j) >*/
		row[ll] += hv[j - 1];
/*<               ll = ll+1 >*/
		++ll;
/*<  150       continue >*/
/* L150: */
	    }
/*<            co = cos(arg) >*/
	    co = cos(arg);
/*<            si = sin(arg) >*/
	    si = sin(arg);
/*<            r = rad(arg) >*/
	    r__ = (*rad)(&arg);
/*<            cs(1) = co*r >*/
	    cs[1] = co * r__;
/*<            cs(2) = si*r >*/
	    cs[2] = si * r__;
/*<            if(iopt2.eq.1) go to 155 >*/
	    if (*iopt2 == 1) {
		goto L155;
	    }
/*<            cs(3) = cs(1)*cs(1) >*/
	    cs[3] = cs[1] * cs[1];
/*<            cs(4) = cs(2)*cs(2) >*/
	    cs[4] = cs[2] * cs[2];
/*<            cs(5) = cs(1)*cs(2) >*/
	    cs[5] = cs[1] * cs[2];
/*<  155       do 170 j=1,nvv >*/
L155:
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
/*<               piv = row(j) >*/
		piv = row[j];
/*<               if(piv.eq.0.) go to 170 >*/
		if (piv == 0.) {
		    goto L170;
		}
/*<               call fpgivs(piv,a(j,1),co,si) >*/
		fpgivs_(&piv, &a[j + a_dim1], &co, &si);
/*<               do 160 l=1,ipar >*/
		i__4 = ipar;
		for (l = 1; l <= i__4; ++l) {
/*<                  call fprota(co,si,cs(l),cosi(l,j)) >*/
		    fprota_(&co, &si, &cs[l], &cosi[l + j * 5]);
/*<  160          continue >*/
/* L160: */
		}
/*<               if(j.eq.nvv) go to 175 >*/
		if (j == nvv) {
		    goto L175;
		}
/*<               j1 = j+1 >*/
		j1 = j + 1;
/*<               j2 = 1 >*/
		j2 = 1;
/*<               do 165 l=j1,nvv >*/
		i__4 = nvv;
		for (l = j1; l <= i__4; ++l) {
/*<                  j2 = j2+1 >*/
		    ++j2;
/*<                  call fprota(co,si,row(l),a(j,j2)) >*/
		    fprota_(&co, &si, &row[l], &a[j + j2 * a_dim1]);
/*<  165          continue >*/
/* L165: */
		}
/*<  170       continue >*/
L170:
		;
	    }
/*<  175    continue >*/
L175:
	    ;
	}
/*<          do 190 l=1,ipar >*/
	i__3 = ipar;
	for (l = 1; l <= i__3; ++l) {
/*<             do 180 j=1,nvv >*/
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
/*<                cs(j) = cosi(l,j) >*/
		cs[j] = cosi[l + j * 5];
/*<  180        continue >*/
/* L180: */
	    }
/*<             call fpback(a,cs,nvv,nvv,cs,ncc) >*/
	    fpback_(&a[a_offset], &cs[1], &nvv, &nvv, &cs[1], ncc);
/*<             do 185 j=1,nvv >*/
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
/*<                cosi(l,j) = cs(j) >*/
		cosi[l + j * 5] = cs[j];
/*<  185        continue >*/
/* L185: */
	    }
/*<  190     continue >*/
/* L190: */
	}
/*  find ncof, the dimension of the spline and ncoff, the number */
/*  of coefficients in the standard b-spline representation. */
/*<  195    nu4 = nu-4 >*/
L195:
	nu4 = *nu - 4;
/*<         nv4 = nv-4 >*/
	nv4 = *nv - 4;
/*<         ncoff = nu4*nv4 >*/
	ncoff = nu4 * nv4;
/*<         ncof = ipar1+nvv*(nu4-1-iopt2-iopt3) >*/
	ncof = ipar1 + nvv * (nu4 - 1 - *iopt2 - *iopt3);
/*  find the bandwidth of the observation matrix a. */
/*<         iband = 4*nvv >*/
	iband = nvv << 2;
/*<         if(nuu-iopt2-iopt3.le.1) iband = ncof >*/
	if (nuu - *iopt2 - *iopt3 <= 1) {
	    iband = ncof;
	}
/*<         iband1 = iband-1 >*/
	iband1 = iband - 1;
/*  initialize the observation matrix a. */
/*<         do 200 i=1,ncof >*/
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           f(i) = 0. >*/
	    f[i__] = 0.;
/*<           do 200 j=1,iband >*/
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
/*<             a(i,j) = 0. >*/
		a[i__ + j * a_dim1] = 0.;
/*<  200    continue >*/
/* L200: */
	    }
	}
/*  initialize the sum of squared residuals. */
/*<         fp = 0. >*/
	*fp = 0.;
/*<         ratio = one+tu(6)/tu(5) >*/
	ratio = one + tu[6] / tu[5];
/*  fetch the data points in the new order. main loop for the */
/*  different panels. */
/*<         do 380 num=1,nreg >*/
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
/*  fix certain constants for the current panel; jrot records the column */
/*  number of the first non-zero element in a row of the observation */
/*  matrix according to a data point of the panel. */
/*<           num1 = num-1 >*/
	    num1 = num - 1;
/*<           lu = num1/nvv >*/
	    lu = num1 / nvv;
/*<           l1 = lu+4 >*/
	    l1 = lu + 4;
/*<           lv = num1-lu*nvv+1 >*/
	    lv = num1 - lu * nvv + 1;
/*<           l2 = lv+3 >*/
	    l2 = lv + 3;
/*<           jrot = 0 >*/
	    jrot = 0;
/*<           if(lu.gt.iopt2) jrot = ipar1+(lu-iopt2-1)*nvv >*/
	    if (lu > *iopt2) {
		jrot = ipar1 + (lu - *iopt2 - 1) * nvv;
	    }
/*<           lu = lu+1 >*/
	    ++lu;
/*  test whether there are still data points in the current panel. */
/*<           in = index(num) >*/
	    in = index[num];
/*<  210      if(in.eq.0) go to 380 >*/
L210:
	    if (in == 0) {
		goto L380;
	    }
/*  fetch a new data point. */
/*<           wi = w(in) >*/
	    wi = w[in];
/*<           zi = z(in)*wi >*/
	    zi = z__[in] * wi;
/*  evaluate for the u-direction, the 4 non-zero b-splines at u(in) */
/*<           call fpbspl(tu,nu,3,u(in),l1,hu) >*/
	    fpbspl_(&tu[1], nu, &c__3, &u[in], &l1, hu);
/*  evaluate for the v-direction, the 4 non-zero b-splines at v(in) */
/*<           call fpbspl(tv,nv,3,v(in),l2,hv) >*/
	    fpbspl_(&tv[1], nv, &c__3, &v[in], &l2, hv);
/*  store the value of these b-splines in spu and spv resp. */
/*<           do 220 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<             spu(in,i) = hu(i) >*/
		spu[in + i__ * spu_dim1] = hu[i__ - 1];
/*<             spv(in,i) = hv(i) >*/
		spv[in + i__ * spv_dim1] = hv[i__ - 1];
/*<  220      continue >*/
/* L220: */
	    }
/*  initialize the new row of observation matrix. */
/*<           do 240 i=1,iband >*/
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             h(i) = 0. >*/
		h__[i__] = 0.;
/*<  240      continue >*/
/* L240: */
	    }
/*  calculate the non-zero elements of the new row by making the cross */
/*  products of the non-zero b-splines in u- and v-direction and */
/*  by taking into account the conditions of the splines. */
/*<           do 250 i=1,nvv >*/
	    i__3 = nvv;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<              row(i) = 0. >*/
		row[i__] = 0.;
/*<  250      continue >*/
/* L250: */
	    }
/*  take into account the periodicity condition of the bicubic splines. */
/*<           ll = lv >*/
	    ll = lv;
/*<           do 260 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<              if(ll.gt.nvv) ll=1 >*/
		if (ll > nvv) {
		    ll = 1;
		}
/*<              row(ll) = row(ll)+hv(i) >*/
		row[ll] += hv[i__ - 1];
/*<              ll = ll+1 >*/
		++ll;
/*<  260      continue >*/
/* L260: */
	    }
/*  take into account the other conditions of the splines. */
/*<           if(iopt2.eq.0 .or. lu.gt.iopt2+1) go to 280 >*/
	    if (*iopt2 == 0 || lu > *iopt2 + 1) {
		goto L280;
	    }
/*<           do 270 l=1,ipar >*/
	    i__3 = ipar;
	    for (l = 1; l <= i__3; ++l) {
/*<              cs(l) = 0. >*/
		cs[l] = 0.;
/*<              do 270 i=1,nvv >*/
		i__4 = nvv;
		for (i__ = 1; i__ <= i__4; ++i__) {
/*<                 cs(l) = cs(l)+row(i)*cosi(l,i) >*/
		    cs[l] += row[i__] * cosi[l + i__ * 5];
/*<  270     continue >*/
/* L270: */
		}
	    }
/*  fill in the non-zero elements of the new row. */
/*<  280     j1 = 0 >*/
L280:
	    j1 = 0;
/*<          do 330 j =1,4 >*/
	    for (j = 1; j <= 4; ++j) {
/*<             jlu = j+lu >*/
		jlu = j + lu;
/*<             huj = hu(j) >*/
		huj = hu[j - 1];
/*<             if(jlu.gt.iopt2+2) go to 320 >*/
		if (jlu > *iopt2 + 2) {
		    goto L320;
		}
/*<             go to (290,290,300,310),jlu >*/
		switch (jlu) {
		    case 1:  goto L290;
		    case 2:  goto L290;
		    case 3:  goto L300;
		    case 4:  goto L310;
		}
/*<  290        h(1) = huj >*/
L290:
		h__[1] = huj;
/*<             j1 = 1 >*/
		j1 = 1;
/*<             go to 330 >*/
		goto L330;
/*<  300        h(1) = h(1)+huj >*/
L300:
		h__[1] += huj;
/*<             h(2) = huj*cs(1) >*/
		h__[2] = huj * cs[1];
/*<             h(3) = huj*cs(2) >*/
		h__[3] = huj * cs[2];
/*<             j1 = 3 >*/
		j1 = 3;
/*<             go to 330 >*/
		goto L330;
/*<  310        h(1) = h(1)+huj >*/
L310:
		h__[1] += huj;
/*<             h(2) = h(2)+huj*ratio*cs(1) >*/
		h__[2] += huj * ratio * cs[1];
/*<             h(3) = h(3)+huj*ratio*cs(2) >*/
		h__[3] += huj * ratio * cs[2];
/*<             h(4) = huj*cs(3) >*/
		h__[4] = huj * cs[3];
/*<             h(5) = huj*cs(4) >*/
		h__[5] = huj * cs[4];
/*<             h(6) = huj*cs(5) >*/
		h__[6] = huj * cs[5];
/*<             j1 = 6 >*/
		j1 = 6;
/*<             go to 330 >*/
		goto L330;
/*<  320        if(jlu.gt.nu4 .and. iopt3.ne.0) go to 330 >*/
L320:
		if (jlu > nu4 && *iopt3 != 0) {
		    goto L330;
		}
/*<             do 325 i=1,nvv >*/
		i__4 = nvv;
		for (i__ = 1; i__ <= i__4; ++i__) {
/*<                j1 = j1+1 >*/
		    ++j1;
/*<                h(j1) = row(i)*huj >*/
		    h__[j1] = row[i__] * huj;
/*<  325        continue >*/
/* L325: */
		}
/*<  330      continue >*/
L330:
		;
	    }
/*<           do 335 i=1,iband >*/
	    i__4 = iband;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<             h(i) = h(i)*wi >*/
		h__[i__] *= wi;
/*<  335      continue >*/
/* L335: */
	    }
/*  rotate the row into triangle by givens transformations. */
/*<           irot = jrot >*/
	    irot = jrot;
/*<           do 350 i=1,iband >*/
	    i__4 = iband;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/*<             irot = irot+1 >*/
		++irot;
/*<             piv = h(i) >*/
		piv = h__[i__];
/*<             if(piv.eq.0.) go to 350 >*/
		if (piv == 0.) {
		    goto L350;
		}
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a(irot,1),co,si) >*/
		fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the right hand side. */
/*<             call fprota(co,si,zi,f(irot)) >*/
		fprota_(&co, &si, &zi, &f[irot]);
/*<             if(i.eq.iband) go to 360 >*/
		if (i__ == iband) {
		    goto L360;
		}
/*  apply that transformation to the left hand side. */
/*<             i2 = 1 >*/
		i2 = 1;
/*<             i3 = i+1 >*/
		i3 = i__ + 1;
/*<             do 340 j=i3,iband >*/
		i__3 = iband;
		for (j = i3; j <= i__3; ++j) {
/*<               i2 = i2+1 >*/
		    ++i2;
/*<               call fprota(co,si,h(j),a(irot,i2)) >*/
		    fprota_(&co, &si, &h__[j], &a[irot + i2 * a_dim1]);
/*<  340        continue >*/
/* L340: */
		}
/*<  350      continue >*/
L350:
		;
	    }
/*  add the contribution of the row to the sum of squares of residual */
/*  right hand sides. */
/*<  360      fp = fp+zi**2 >*/
L360:
/* Computing 2nd power */
	    d__1 = zi;
	    *fp += d__1 * d__1;
/*  find the number of the next data point in the panel. */
/*<  370      in = nummer(in) >*/
/* L370: */
	    in = nummer[in];
/*<           go to 210 >*/
	    goto L210;
/*<  380    continue >*/
L380:
	    ;
	}
/*  find dmax, the maximum value for the diagonal elements in the reduced */
/*  triangle. */
/*<         dmax = 0. >*/
	dmax__ = 0.;
/*<         do 390 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(a(i,1).le.dmax) go to 390 >*/
	    if (a[i__ + a_dim1] <= dmax__) {
		goto L390;
	    }
/*<           dmax = a(i,1) >*/
	    dmax__ = a[i__ + a_dim1];
/*<  390    continue >*/
L390:
	    ;
	}
/*  check whether the observation matrix is rank deficient. */
/*<         sigma = eps*dmax >*/
	sigma = eps * dmax__;
/*<         do 400 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(a(i,1).le.sigma) go to 410 >*/
	    if (a[i__ + a_dim1] <= sigma) {
		goto L410;
	    }
/*<  400    continue >*/
/* L400: */
	}
/*  backward substitution in case of full rank. */
/*<         call fpback(a,f,ncof,iband,c,ncc) >*/
	fpback_(&a[a_offset], &f[1], &ncof, &iband, &c__[1], ncc);
/*<         rank = ncof >*/
	rank = ncof;
/*<         do 405 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           q(i,1) = a(i,1)/dmax >*/
	    q[i__ + q_dim1] = a[i__ + a_dim1] / dmax__;
/*<  405    continue >*/
/* L405: */
	}
/*<         go to 430 >*/
	goto L430;
/*  in case of rank deficiency, find the minimum norm solution. */
/*<  410    lwest = ncof*iband+ncof+iband >*/
L410:
	lwest = ncof * iband + ncof + iband;
/*<         if(lwrk.lt.lwest) go to 925 >*/
	if (*lwrk < lwest) {
	    goto L925;
	}
/*<         lf = 1 >*/
	lf = 1;
/*<         lh = lf+ncof >*/
	lh = lf + ncof;
/*<         la = lh+iband >*/
	la = lh + iband;
/*<         do 420 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           ff(i) = f(i) >*/
	    ff[i__] = f[i__];
/*<           do 420 j=1,iband >*/
	    i__4 = iband;
	    for (j = 1; j <= i__4; ++j) {
/*<             q(i,j) = a(i,j) >*/
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/*<  420    continue >*/
/* L420: */
	    }
	}
/*<    >*/
	fprank_(&q[q_offset], &ff[1], &ncof, &iband, ncc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
/*<         do 425 i=1,ncof >*/
	i__4 = ncof;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<           q(i,1) = q(i,1)/dmax >*/
	    q[i__ + q_dim1] /= dmax__;
/*<  425    continue >*/
/* L425: */
	}
/*  add to the sum of squared residuals, the contribution of reducing */
/*  the rank. */
/*<         fp = fp+sq >*/
	*fp += sq;
/*  find the coefficients in the standard b-spline representation of */
/*  the spline. */
/*<  430    call fprppo(nu,nv,iopt2,iopt3,cosi,ratio,c,ff,ncoff) >*/
L430:
	fprppo_(nu, nv, iopt2, iopt3, &cosi[6], &ratio, &c__[1], &ff[1], &
		ncoff);
/*  test whether the least-squares spline is an acceptable solution. */
/*<         if(iopt1.lt.0) if(fp) 970,970,980 >*/
	if (*iopt1 < 0) {
	    if (*fp <= 0.) {
		goto L970;
	    } else {
		goto L980;
	    }
	}
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).le.acc) if(fp) 970,970,980 >*/
	if (abs(fpms) <= acc) {
	    if (*fp <= 0.) {
		goto L970;
	    } else {
		goto L980;
	    }
	}
/*  if f(p=inf) < s, accept the choice of knots. */
/*<         if(fpms.lt.0.) go to 580 >*/
	if (fpms < 0.) {
	    goto L580;
	}
/*  test whether we cannot further increase the number of knots */
/*<         if(m.lt.ncof) go to 935 >*/
	if (*m < ncof) {
	    goto L935;
	}
/*  search where to add a new knot. */
/*  find for each interval the sum of squared residuals fpint for the */
/*  data points having the coordinate belonging to that knot interval. */
/*  calculate also coord which is the same sum, weighted by the position */
/*  of the data points considered. */
/*<  440    do 450 i=1,nrint >*/
/* L440: */
	i__4 = nrint;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<           fpint(i) = 0. >*/
	    fpint[i__] = 0.;
/*<           coord(i) = 0. >*/
	    coord[i__] = 0.;
/*<  450    continue >*/
/* L450: */
	}
/*<         do 490 num=1,nreg >*/
	i__4 = nreg;
	for (num = 1; num <= i__4; ++num) {
/*<           num1 = num-1 >*/
	    num1 = num - 1;
/*<           lu = num1/nvv >*/
	    lu = num1 / nvv;
/*<           l1 = lu+1 >*/
	    l1 = lu + 1;
/*<           lv = num1-lu*nvv >*/
	    lv = num1 - lu * nvv;
/*<           l2 = lv+1+nuu >*/
	    l2 = lv + 1 + nuu;
/*<           jrot = lu*nv4+lv >*/
	    jrot = lu * nv4 + lv;
/*<           in = index(num) >*/
	    in = index[num];
/*<  460      if(in.eq.0) go to 490 >*/
L460:
	    if (in == 0) {
		goto L490;
	    }
/*<           store = 0. >*/
	    store = 0.;
/*<           i1 = jrot >*/
	    i1 = jrot;
/*<           do 480 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<             hui = spu(in,i) >*/
		hui = spu[in + i__ * spu_dim1];
/*<             j1 = i1 >*/
		j1 = i1;
/*<             do 470 j=1,4 >*/
		for (j = 1; j <= 4; ++j) {
/*<               j1 = j1+1 >*/
		    ++j1;
/*<               store = store+hui*spv(in,j)*c(j1) >*/
		    store += hui * spv[in + j * spv_dim1] * c__[j1];
/*<  470        continue >*/
/* L470: */
		}
/*<             i1 = i1+nv4 >*/
		i1 += nv4;
/*<  480      continue >*/
/* L480: */
	    }
/*<           store = (w(in)*(z(in)-store))**2 >*/
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    store = d__1 * d__1;
/*<           fpint(l1) = fpint(l1)+store >*/
	    fpint[l1] += store;
/*<           coord(l1) = coord(l1)+store*u(in) >*/
	    coord[l1] += store * u[in];
/*<           fpint(l2) = fpint(l2)+store >*/
	    fpint[l2] += store;
/*<           coord(l2) = coord(l2)+store*v(in) >*/
	    coord[l2] += store * v[in];
/*<           in = nummer(in) >*/
	    in = nummer[in];
/*<           go to 460 >*/
	    goto L460;
/*<  490    continue >*/
L490:
	    ;
	}
/* bring together the information concerning knot panels which are */
/* symmetric with respect to the origin. */
/*<         do 495 i=1,nrr >*/
	i__4 = nrr;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<           l1 = nuu+i >*/
	    l1 = nuu + i__;
/*<           l2 = l1+nrr >*/
	    l2 = l1 + nrr;
/*<           fpint(l1) = fpint(l1)+fpint(l2) >*/
	    fpint[l1] += fpint[l2];
/*<           coord(l1) = coord(l1)+coord(l2)-pi*fpint(l2) >*/
	    coord[l1] = coord[l1] + coord[l2] - pi * fpint[l2];
/*<  495    continue >*/
/* L495: */
	}
/*  find the interval for which fpint is maximal on the condition that */
/*  there still can be added a knot. */
/*<         l1 = 1 >*/
	l1 = 1;
/*<         l2 = nuu+nrr >*/
	l2 = nuu + nrr;
/*<         if(nuest.lt.nu+1) l1=nuu+1 >*/
	if (*nuest < *nu + 1) {
	    l1 = nuu + 1;
	}
/*<         if(nvest.lt.nv+2) l2=nuu >*/
	if (*nvest < *nv + 2) {
	    l2 = nuu;
	}
/*  test whether we cannot further increase the number of knots. */
/*<         if(l1.gt.l2) go to 950 >*/
	if (l1 > l2) {
	    goto L950;
	}
/*<  500    fpmax = 0. >*/
L500:
	fpmax = 0.;
/*<         l = 0 >*/
	l = 0;
/*<         do 510 i=l1,l2 >*/
	i__4 = l2;
	for (i__ = l1; i__ <= i__4; ++i__) {
/*<           if(fpmax.ge.fpint(i)) go to 510 >*/
	    if (fpmax >= fpint[i__]) {
		goto L510;
	    }
/*<           l = i >*/
	    l = i__;
/*<           fpmax = fpint(i) >*/
	    fpmax = fpint[i__];
/*<  510    continue >*/
L510:
	    ;
	}
/*<         if(l.eq.0) go to 930 >*/
	if (l == 0) {
	    goto L930;
	}
/*  calculate the position of the new knot. */
/*<         arg = coord(l)/fpint(l) >*/
	arg = coord[l] / fpint[l];
/*  test in what direction the new knot is going to be added. */
/*<         if(l.gt.nuu) go to 530 >*/
	if (l > nuu) {
	    goto L530;
	}
/*  addition in the u-direction */
/*<         l4 = l+4 >*/
	l4 = l + 4;
/*<         fpint(l) = 0. >*/
	fpint[l] = 0.;
/*<         fac1 = tu(l4)-arg >*/
	fac1 = tu[l4] - arg;
/*<         fac2 = arg-tu(l4-1) >*/
	fac2 = arg - tu[l4 - 1];
/*<         if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500 >*/
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
/*<         j = nu >*/
	j = *nu;
/*<         do 520 i=l4,nu >*/
	i__4 = *nu;
	for (i__ = l4; i__ <= i__4; ++i__) {
/*<           tu(j+1) = tu(j) >*/
	    tu[j + 1] = tu[j];
/*<           j = j-1 >*/
	    --j;
/*<  520    continue >*/
/* L520: */
	}
/*<         tu(l4) = arg >*/
	tu[l4] = arg;
/*<         nu = nu+1 >*/
	++(*nu);
/*<         go to 570 >*/
	goto L570;
/*  addition in the v-direction */
/*<  530    l4 = l+4-nuu >*/
L530:
	l4 = l + 4 - nuu;
/*<         fpint(l) = 0. >*/
	fpint[l] = 0.;
/*<         fac1 = tv(l4)-arg >*/
	fac1 = tv[l4] - arg;
/*<         fac2 = arg-tv(l4-1) >*/
	fac2 = arg - tv[l4 - 1];
/*<         if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500 >*/
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
/*<         ll = nrr+4 >*/
	ll = nrr + 4;
/*<         j = ll >*/
	j = ll;
/*<         do 550 i=l4,ll >*/
	i__4 = ll;
	for (i__ = l4; i__ <= i__4; ++i__) {
/*<           tv(j+1) = tv(j) >*/
	    tv[j + 1] = tv[j];
/*<           j = j-1 >*/
	    --j;
/*<  550    continue >*/
/* L550: */
	}
/*<         tv(l4) = arg >*/
	tv[l4] = arg;
/*<         nv = nv+2 >*/
	*nv += 2;
/*<         nrr = nrr+1 >*/
	++nrr;
/*<         do 560 i=5,ll >*/
	i__4 = ll;
	for (i__ = 5; i__ <= i__4; ++i__) {
/*<           j = i+nrr >*/
	    j = i__ + nrr;
/*<           tv(j) = tv(i)+pi >*/
	    tv[j] = tv[i__] + pi;
/*<  560    continue >*/
/* L560: */
	}
/*  restart the computations with the new set of knots. */
/*<  570  continue >*/
L570:
	;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing bicubic spline.               c */
/* ******************************************************               c */
/* we have determined the number of knots and their position. we now    c */
/* compute the coefficients of the smoothing spline sp(u,v).            c */
/* the observation matrix a is extended by the rows of a matrix, expres-c */
/* sing that sp(u,v) must be a constant function in the variable        c */
/* v and a cubic polynomial in the variable u. the corresponding        c */
/* weights of these additional rows are set to 1/(p). iteratively       c */
/* we than have to determine the value of p such that f(p) = sum((w(i)* c */
/* (z(i)-sp(u(i),v(i))))**2)  be = s.                                   c */
/* we already know that the least-squares polynomial corresponds to p=0,c */
/* and that the least-squares bicubic spline corresponds to p=infin.    c */
/* the iteration process makes use of rational interpolation. since f(p)c */
/* is a convex and strictly decreasing function of p, it can be approx- c */
/* imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c */
/* three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c */
/* f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c */
/* of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tu(l),l=5,...,nu-4. */
/*<  580  call fpdisc(tu,nu,5,bu,nuest) >*/
L580:
    fpdisc_(&tu[1], nu, &c__5, &bu[bu_offset], nuest);
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tv(l),l=5,...,nv-4. */
/*<       call fpdisc(tv,nv,5,bv,nvest) >*/
    fpdisc_(&tv[1], nv, &c__5, &bv[bv_offset], nvest);
/*  initial value for p. */
/*<       p1 = 0. >*/
    p1 = 0.;
/*<       f1 = sup-s >*/
    f1 = *sup - *s;
/*<       p3 = -one >*/
    p3 = -one;
/*<       f3 = fpms >*/
    f3 = fpms;
/*<       p = 0. >*/
    p = 0.;
/*<       do 590 i=1,ncof >*/
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         p = p+a(i,1) >*/
	p += a[i__ + a_dim1];
/*<  590  continue >*/
/* L590: */
    }
/*<       rn = ncof >*/
    rn = (doublereal) ncof;
/*<       p = rn/p >*/
    p = rn / p;
/*  find the bandwidth of the extended observation matrix. */
/*<       iband4 = iband+ipar1 >*/
    iband4 = iband + ipar1;
/*<       if(iband4.gt.ncof) iband4 = ncof >*/
    if (iband4 > ncof) {
	iband4 = ncof;
    }
/*<       iband3 = iband4 -1 >*/
    iband3 = iband4 - 1;
/*<       ich1 = 0 >*/
    ich1 = 0;
/*<       ich3 = 0 >*/
    ich3 = 0;
/*<       nuu = nu4-iopt3-1 >*/
    nuu = nu4 - *iopt3 - 1;
/*  iteration process to find the root of f(p)=s. */
/*<       do 920 iter=1,maxit >*/
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*<         pinv = one/p >*/
	pinv = one / p;
/*  store the triangularized observation matrix into q. */
/*<         do 630 i=1,ncof >*/
	i__4 = ncof;
	for (i__ = 1; i__ <= i__4; ++i__) {
/*<           ff(i) = f(i) >*/
	    ff[i__] = f[i__];
/*<           do 620 j=1,iband4 >*/
	    i__2 = iband4;
	    for (j = 1; j <= i__2; ++j) {
/*<             q(i,j) = 0. >*/
		q[i__ + j * q_dim1] = 0.;
/*<  620      continue >*/
/* L620: */
	    }
/*<           do 630 j=1,iband >*/
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
/*<             q(i,j) = a(i,j) >*/
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/*<  630    continue >*/
/* L630: */
	    }
	}
/*  extend the observation matrix with the rows of a matrix, expressing */
/*  that for u=constant sp(u,v) must be a constant function. */
/*<         do 720 i=5,nv4 >*/
	i__2 = nv4;
	for (i__ = 5; i__ <= i__2; ++i__) {
/*<           ii = i-4 >*/
	    ii = i__ - 4;
/*<           do 635 l=1,nvv >*/
	    i__4 = nvv;
	    for (l = 1; l <= i__4; ++l) {
/*<              row(l) = 0. >*/
		row[l] = 0.;
/*<  635      continue >*/
/* L635: */
	    }
/*<           ll = ii >*/
	    ll = ii;
/*<           do 640  l=1,5 >*/
	    for (l = 1; l <= 5; ++l) {
/*<              if(ll.gt.nvv) ll=1 >*/
		if (ll > nvv) {
		    ll = 1;
		}
/*<              row(ll) = row(ll)+bv(ii,l) >*/
		row[ll] += bv[ii + l * bv_dim1];
/*<              ll = ll+1 >*/
		++ll;
/*<  640      continue >*/
/* L640: */
	    }
/*<           do 720 j=1,nuu >*/
	    i__4 = nuu;
	    for (j = 1; j <= i__4; ++j) {
/*  initialize the new row. */
/*<             do 645 l=1,iband >*/
		i__3 = iband;
		for (l = 1; l <= i__3; ++l) {
/*<               h(l) = 0. >*/
		    h__[l] = 0.;
/*<  645        continue >*/
/* L645: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
/*<             if(j.gt.iopt2) go to 665 >*/
		if (j > *iopt2) {
		    goto L665;
		}
/*<             if(j.eq.2) go to 655 >*/
		if (j == 2) {
		    goto L655;
		}
/*<             do 650 k=1,2 >*/
		for (k = 1; k <= 2; ++k) {
/*<                cs(k) = 0. >*/
		    cs[k] = 0.;
/*<                do 650 l=1,nvv >*/
		    i__3 = nvv;
		    for (l = 1; l <= i__3; ++l) {
/*<                   cs(k) = cs(k)+cosi(k,l)*row(l) >*/
			cs[k] += cosi[k + l * 5] * row[l];
/*<  650        continue >*/
/* L650: */
		    }
		}
/*<             h(1) = cs(1) >*/
		h__[1] = cs[1];
/*<             h(2) = cs(2) >*/
		h__[2] = cs[2];
/*<             jrot = 2 >*/
		jrot = 2;
/*<             go to 675 >*/
		goto L675;
/*<  655        do 660 k=3,5 >*/
L655:
		for (k = 3; k <= 5; ++k) {
/*<                cs(k) = 0. >*/
		    cs[k] = 0.;
/*<                do 660 l=1,nvv >*/
		    i__3 = nvv;
		    for (l = 1; l <= i__3; ++l) {
/*<                   cs(k) = cs(k)+cosi(k,l)*row(l) >*/
			cs[k] += cosi[k + l * 5] * row[l];
/*<  660        continue >*/
/* L660: */
		    }
		}
/*<             h(1) = cs(1)*ratio >*/
		h__[1] = cs[1] * ratio;
/*<             h(2) = cs(2)*ratio >*/
		h__[2] = cs[2] * ratio;
/*<             h(3) = cs(3) >*/
		h__[3] = cs[3];
/*<             h(4) = cs(4) >*/
		h__[4] = cs[4];
/*<             h(5) = cs(5) >*/
		h__[5] = cs[5];
/*<             jrot = 2 >*/
		jrot = 2;
/*<             go to 675 >*/
		goto L675;
/*<  665        do 670 l=1,nvv >*/
L665:
		i__3 = nvv;
		for (l = 1; l <= i__3; ++l) {
/*<                h(l) = row(l) >*/
		    h__[l] = row[l];
/*<  670        continue >*/
/* L670: */
		}
/*<             jrot = ipar1+1+(j-iopt2-1)*nvv >*/
		jrot = ipar1 + 1 + (j - *iopt2 - 1) * nvv;
/*<  675        do 677 l=1,iband >*/
L675:
		i__3 = iband;
		for (l = 1; l <= i__3; ++l) {
/*<               h(l) = h(l)*pinv >*/
		    h__[l] *= pinv;
/*<  677        continue >*/
/* L677: */
		}
/*<             zi = 0. >*/
		zi = 0.;
/*  rotate the new row into triangle by givens transformations. */
/*<             do 710 irot=jrot,ncof >*/
		i__3 = ncof;
		for (irot = jrot; irot <= i__3; ++irot) {
/*<               piv = h(1) >*/
		    piv = h__[1];
/*<               i2 = min0(iband1,ncof-irot) >*/
/* Computing MIN */
		    i__5 = iband1, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
/*<               if(piv.eq.0.) if(i2) 720,720,690 >*/
		    if (piv == 0.) {
			if (i2 <= 0) {
			    goto L720;
			} else {
			    goto L690;
			}
		    }
/*  calculate the parameters of the givens transformation. */
/*<               call fpgivs(piv,q(irot,1),co,si) >*/
		    fpgivs_(&piv, &q[irot + q_dim1], &co, &si);
/*  apply that givens transformation to the right hand side. */
/*<               call fprota(co,si,zi,ff(irot)) >*/
		    fprota_(&co, &si, &zi, &ff[irot]);
/*<               if(i2.eq.0) go to 720 >*/
		    if (i2 == 0) {
			goto L720;
		    }
/*  apply that givens transformation to the left hand side. */
/*<               do 680 l=1,i2 >*/
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
/*<                 l1 = l+1 >*/
			l1 = l + 1;
/*<                 call fprota(co,si,h(l1),q(irot,l1)) >*/
			fprota_(&co, &si, &h__[l1], &q[irot + l1 * q_dim1]);
/*<  680          continue >*/
/* L680: */
		    }
/*<  690          do 700 l=1,i2 >*/
L690:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
/*<                 h(l) = h(l+1) >*/
			h__[l] = h__[l + 1];
/*<  700          continue >*/
/* L700: */
		    }
/*<               h(i2+1) = 0. >*/
		    h__[i2 + 1] = 0.;
/*<  710        continue >*/
/* L710: */
		}
/*<  720    continue >*/
L720:
		;
	    }
	}
/*  extend the observation matrix with the rows of a matrix expressing */
/*  that for v=constant. sp(u,v) must be a cubic polynomial. */
/*<         do 810 i=5,nu4 >*/
	i__4 = nu4;
	for (i__ = 5; i__ <= i__4; ++i__) {
/*<           ii = i-4 >*/
	    ii = i__ - 4;
/*<           do 810 j=1,nvv >*/
	    i__2 = nvv;
	    for (j = 1; j <= i__2; ++j) {
/*  initialize the new row */
/*<             do 730 l=1,iband4 >*/
		i__3 = iband4;
		for (l = 1; l <= i__3; ++l) {
/*<               h(l) = 0. >*/
		    h__[l] = 0.;
/*<  730        continue >*/
/* L730: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
/*<             j1 = 1 >*/
		j1 = 1;
/*<             do 760 l=1,5 >*/
		for (l = 1; l <= 5; ++l) {
/*<                il = ii+l-1 >*/
		    il = ii + l - 1;
/*<                if(il.eq.nu4 .and. iopt3.ne.0) go to 760 >*/
		    if (il == nu4 && *iopt3 != 0) {
			goto L760;
		    }
/*<                if(il.gt.iopt2+1) go to 750 >*/
		    if (il > *iopt2 + 1) {
			goto L750;
		    }
/*<                go to (735,740,745),il >*/
		    switch (il) {
			case 1:  goto L735;
			case 2:  goto L740;
			case 3:  goto L745;
		    }
/*<  735           h(1) = bu(ii,l) >*/
L735:
		    h__[1] = bu[ii + l * bu_dim1];
/*<                j1 = j+1 >*/
		    j1 = j + 1;
/*<                go to 760 >*/
		    goto L760;
/*<  740           h(1) = h(1)+bu(ii,l) >*/
L740:
		    h__[1] += bu[ii + l * bu_dim1];
/*<                h(2) = bu(ii,l)*cosi(1,j) >*/
		    h__[2] = bu[ii + l * bu_dim1] * cosi[j * 5 + 1];
/*<                h(3) = bu(ii,l)*cosi(2,j) >*/
		    h__[3] = bu[ii + l * bu_dim1] * cosi[j * 5 + 2];
/*<                j1 = j+3 >*/
		    j1 = j + 3;
/*<                go to 760 >*/
		    goto L760;
/*<  745           h(1) = h(1)+bu(ii,l) >*/
L745:
		    h__[1] += bu[ii + l * bu_dim1];
/*<                h(2) = bu(ii,l)*cosi(1,j)*ratio >*/
		    h__[2] = bu[ii + l * bu_dim1] * cosi[j * 5 + 1] * ratio;
/*<                h(3) = bu(ii,l)*cosi(2,j)*ratio >*/
		    h__[3] = bu[ii + l * bu_dim1] * cosi[j * 5 + 2] * ratio;
/*<                h(4) = bu(ii,l)*cosi(3,j) >*/
		    h__[4] = bu[ii + l * bu_dim1] * cosi[j * 5 + 3];
/*<                h(5) = bu(ii,l)*cosi(4,j) >*/
		    h__[5] = bu[ii + l * bu_dim1] * cosi[j * 5 + 4];
/*<                h(6) = bu(ii,l)*cosi(5,j) >*/
		    h__[6] = bu[ii + l * bu_dim1] * cosi[j * 5 + 5];
/*<                j1 = j+6 >*/
		    j1 = j + 6;
/*<                go to 760 >*/
		    goto L760;
/*<  750           h(j1) = bu(ii,l) >*/
L750:
		    h__[j1] = bu[ii + l * bu_dim1];
/*<                j1 = j1+nvv >*/
		    j1 += nvv;
/*<  760        continue >*/
L760:
		    ;
		}
/*<             do 765 l=1,iband4 >*/
		i__3 = iband4;
		for (l = 1; l <= i__3; ++l) {
/*<               h(l) = h(l)*pinv >*/
		    h__[l] *= pinv;
/*<  765        continue >*/
/* L765: */
		}
/*<             zi = 0. >*/
		zi = 0.;
/*<             jrot = 1 >*/
		jrot = 1;
/*<             if(ii.gt.iopt2+1) jrot = ipar1+(ii-iopt2-2)*nvv+j >*/
		if (ii > *iopt2 + 1) {
		    jrot = ipar1 + (ii - *iopt2 - 2) * nvv + j;
		}
/*  rotate the new row into triangle by givens transformations. */
/*<             do 800 irot=jrot,ncof >*/
		i__3 = ncof;
		for (irot = jrot; irot <= i__3; ++irot) {
/*<               piv = h(1) >*/
		    piv = h__[1];
/*<               i2 = min0(iband3,ncof-irot) >*/
/* Computing MIN */
		    i__5 = iband3, i__6 = ncof - irot;
		    i2 = min(i__5,i__6);
/*<               if(piv.eq.0.) if(i2) 810,810,780 >*/
		    if (piv == 0.) {
			if (i2 <= 0) {
			    goto L810;
			} else {
			    goto L780;
			}
		    }
/*  calculate the parameters of the givens transformation. */
/*<               call fpgivs(piv,q(irot,1),co,si) >*/
		    fpgivs_(&piv, &q[irot + q_dim1], &co, &si);
/*  apply that givens transformation to the right hand side. */
/*<               call fprota(co,si,zi,ff(irot)) >*/
		    fprota_(&co, &si, &zi, &ff[irot]);
/*<               if(i2.eq.0) go to 810 >*/
		    if (i2 == 0) {
			goto L810;
		    }
/*  apply that givens transformation to the left hand side. */
/*<               do 770 l=1,i2 >*/
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
/*<                 l1 = l+1 >*/
			l1 = l + 1;
/*<                 call fprota(co,si,h(l1),q(irot,l1)) >*/
			fprota_(&co, &si, &h__[l1], &q[irot + l1 * q_dim1]);
/*<  770          continue >*/
/* L770: */
		    }
/*<  780          do 790 l=1,i2 >*/
L780:
		    i__5 = i2;
		    for (l = 1; l <= i__5; ++l) {
/*<                 h(l) = h(l+1) >*/
			h__[l] = h__[l + 1];
/*<  790          continue >*/
/* L790: */
		    }
/*<               h(i2+1) = 0. >*/
		    h__[i2 + 1] = 0.;
/*<  800        continue >*/
/* L800: */
		}
/*<  810    continue >*/
L810:
		;
	    }
	}
/*  find dmax, the maximum value for the diagonal elements in the */
/*  reduced triangle. */
/*<         dmax = 0. >*/
	dmax__ = 0.;
/*<         do 820 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(q(i,1).le.dmax) go to 820 >*/
	    if (q[i__ + q_dim1] <= dmax__) {
		goto L820;
	    }
/*<           dmax = q(i,1) >*/
	    dmax__ = q[i__ + q_dim1];
/*<  820    continue >*/
L820:
	    ;
	}
/*  check whether the matrix is rank deficient. */
/*<         sigma = eps*dmax >*/
	sigma = eps * dmax__;
/*<         do 830 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(q(i,1).le.sigma) go to 840 >*/
	    if (q[i__ + q_dim1] <= sigma) {
		goto L840;
	    }
/*<  830    continue >*/
/* L830: */
	}
/*  backward substitution in case of full rank. */
/*<         call fpback(q,ff,ncof,iband4,c,ncc) >*/
	fpback_(&q[q_offset], &ff[1], &ncof, &iband4, &c__[1], ncc);
/*<         rank = ncof >*/
	rank = ncof;
/*<         go to 845 >*/
	goto L845;
/*  in case of rank deficiency, find the minimum norm solution. */
/*<  840    lwest = ncof*iband4+ncof+iband4 >*/
L840:
	lwest = ncof * iband4 + ncof + iband4;
/*<         if(lwrk.lt.lwest) go to 925 >*/
	if (*lwrk < lwest) {
	    goto L925;
	}
/*<         lf = 1 >*/
	lf = 1;
/*<         lh = lf+ncof >*/
	lh = lf + ncof;
/*<         la = lh+iband4 >*/
	la = lh + iband4;
/*<    >*/
	fprank_(&q[q_offset], &ff[1], &ncof, &iband4, ncc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
/*<  845    do 850 i=1,ncof >*/
L845:
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<            q(i,1) = q(i,1)/dmax >*/
	    q[i__ + q_dim1] /= dmax__;
/*<  850    continue >*/
/* L850: */
	}
/*  find the coefficients in the standard b-spline representation of */
/*  the polar spline. */
/*<         call fprppo(nu,nv,iopt2,iopt3,cosi,ratio,c,ff,ncoff) >*/
	fprppo_(nu, nv, iopt2, iopt3, &cosi[6], &ratio, &c__[1], &ff[1], &
		ncoff);
/*  compute f(p). */
/*<         fp = 0. >*/
	*fp = 0.;
/*<         do 890 num = 1,nreg >*/
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
/*<           num1 = num-1 >*/
	    num1 = num - 1;
/*<           lu = num1/nvv >*/
	    lu = num1 / nvv;
/*<           lv = num1-lu*nvv >*/
	    lv = num1 - lu * nvv;
/*<           jrot = lu*nv4+lv >*/
	    jrot = lu * nv4 + lv;
/*<           in = index(num) >*/
	    in = index[num];
/*<  860      if(in.eq.0) go to 890 >*/
L860:
	    if (in == 0) {
		goto L890;
	    }
/*<           store = 0. >*/
	    store = 0.;
/*<           i1 = jrot >*/
	    i1 = jrot;
/*<           do 880 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<             hui = spu(in,i) >*/
		hui = spu[in + i__ * spu_dim1];
/*<             j1 = i1 >*/
		j1 = i1;
/*<             do 870 j=1,4 >*/
		for (j = 1; j <= 4; ++j) {
/*<               j1 = j1+1 >*/
		    ++j1;
/*<               store = store+hui*spv(in,j)*c(j1) >*/
		    store += hui * spv[in + j * spv_dim1] * c__[j1];
/*<  870        continue >*/
/* L870: */
		}
/*<             i1 = i1+nv4 >*/
		i1 += nv4;
/*<  880      continue >*/
/* L880: */
	    }
/*<           fp = fp+(w(in)*(z(in)-store))**2 >*/
/* Computing 2nd power */
	    d__1 = w[in] * (z__[in] - store);
	    *fp += d__1 * d__1;
/*<           in = nummer(in) >*/
	    in = nummer[in];
/*<           go to 860 >*/
	    goto L860;
/*<  890    continue >*/
L890:
	    ;
	}
/*  test whether the approximation sp(u,v) is an acceptable solution */
/*<         fpms = fp-s >*/
	fpms = *fp - *s;
/*<         if(abs(fpms).le.acc) go to 980 >*/
	if (abs(fpms) <= acc) {
	    goto L980;
	}
/*  test whether the maximum allowable number of iterations has been */
/*  reached. */
/*<         if(iter.eq.maxit) go to 940 >*/
	if (iter == *maxit) {
	    goto L940;
	}
/*  carry out one more step of the iteration process. */
/*<         p2 = p >*/
	p2 = p;
/*<         f2 = fpms >*/
	f2 = fpms;
/*<         if(ich3.ne.0) go to 900 >*/
	if (ich3 != 0) {
	    goto L900;
	}
/*<         if((f2-f3).gt.acc) go to 895 >*/
	if (f2 - f3 > acc) {
	    goto L895;
	}
/*  our initial choice of p is too large. */
/*<         p3 = p2 >*/
	p3 = p2;
/*<         f3 = f2 >*/
	f3 = f2;
/*<         p = p*con4 >*/
	p *= con4;
/*<         if(p.le.p1) p = p1*con9 + p2*con1 >*/
	if (p <= p1) {
	    p = p1 * con9 + p2 * con1;
	}
/*<         go to 920 >*/
	goto L920;
/*<  895    if(f2.lt.0.) ich3 = 1 >*/
L895:
	if (f2 < 0.) {
	    ich3 = 1;
	}
/*<  900    if(ich1.ne.0) go to 910 >*/
L900:
	if (ich1 != 0) {
	    goto L910;
	}
/*<         if((f1-f2).gt.acc) go to 905 >*/
	if (f1 - f2 > acc) {
	    goto L905;
	}
/*  our initial choice of p is too small */
/*<         p1 = p2 >*/
	p1 = p2;
/*<         f1 = f2 >*/
	f1 = f2;
/*<         p = p/con4 >*/
	p /= con4;
/*<         if(p3.lt.0.) go to 920 >*/
	if (p3 < 0.) {
	    goto L920;
	}
/*<         if(p.ge.p3) p = p2*con1 +p3*con9 >*/
	if (p >= p3) {
	    p = p2 * con1 + p3 * con9;
	}
/*<         go to 920 >*/
	goto L920;
/*<  905    if(f2.gt.0.) ich1 = 1 >*/
L905:
	if (f2 > 0.) {
	    ich1 = 1;
	}
/*  test whether the iteration process proceeds as theoretically */
/*  expected. */
/*<  910    if(f2.ge.f1 .or. f2.le.f3) go to 945 >*/
L910:
	if (f2 >= f1 || f2 <= f3) {
	    goto L945;
	}
/*  find the new value of p. */
/*<         p = fprati(p1,f1,p2,f2,p3,f3) >*/
	p = fprati_(&p1, &f1, &p2, &f2, &p3, &f3);
/*<  920  continue >*/
L920:
	;
    }
/*  error codes and messages. */
/*<  925  ier = lwest >*/
L925:
    *ier = lwest;
/*<       go to 990 >*/
    goto L990;
/*<  930  ier = 5 >*/
L930:
    *ier = 5;
/*<       go to 990 >*/
    goto L990;
/*<  935  ier = 4 >*/
L935:
    *ier = 4;
/*<       go to 990 >*/
    goto L990;
/*<  940  ier = 3 >*/
L940:
    *ier = 3;
/*<       go to 990 >*/
    goto L990;
/*<  945  ier = 2 >*/
L945:
    *ier = 2;
/*<       go to 990 >*/
    goto L990;
/*<  950  ier = 1 >*/
L950:
    *ier = 1;
/*<       go to 990 >*/
    goto L990;
/*<  960  ier = -2 >*/
L960:
    *ier = -2;
/*<       go to 990 >*/
    goto L990;
/*<  970  ier = -1 >*/
L970:
    *ier = -1;
/*<       fp = 0. >*/
    *fp = 0.;
/*<  980  if(ncof.ne.rank) ier = -rank >*/
L980:
    if (ncof != rank) {
	*ier = -rank;
    }
/*<  990  return >*/
L990:
    return 0;
/*<       end >*/
} /* fppola_ */

