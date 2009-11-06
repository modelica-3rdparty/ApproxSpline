/* fpsphe.f -- translated by f2c (version 20061008).
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
/* Subroutine */ int fpsphe_(integer *iopt, integer *m, doublereal *teta, 
	doublereal *phi, doublereal *r__, doublereal *w, doublereal *s, 
	integer *ntest, integer *npest, doublereal *eta, doublereal *tol, 
	integer *maxit, integer *ib1, integer *ib3, integer *nc, integer *ncc,
	 integer *intest, integer *nrest, integer *nt, doublereal *tt, 
	integer *np, doublereal *tp, doublereal *c__, doublereal *fp, 
	doublereal *sup, doublereal *fpint, doublereal *coord, doublereal *f, 
	doublereal *ff, doublereal *row, doublereal *coco, doublereal *cosi, 
	doublereal *a, doublereal *q, doublereal *bt, doublereal *bp, 
	doublereal *spt, doublereal *spp, doublereal *h__, integer *index, 
	integer *nummer, doublereal *wrk, integer *lwrk, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, bt_dim1, bt_offset, bp_dim1, 
	    bp_offset, spt_dim1, spt_offset, spp_dim1, spp_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), cos(doublereal), sin(
	    doublereal);

    /* Local variables */
    static integer i__, j, l;
    static doublereal p, c1, d1, d2, f1, f2, f3;
    static integer i1, i2, i3, j1, j2, l1, l2;
    static doublereal p1, p2, p3;
    static integer l3, l4;
    static doublereal aa;
    static integer la;
    static doublereal cn, co, fn;
    static integer ii;
    static doublereal pi;
    static integer ij;
    static doublereal ri, si;
    static integer il, in;
    static doublereal wi, rn;
    static integer lf;
    static doublereal sq;
    static integer lh, ll, lp, lt;
    static doublereal ht[4], hp[4], pi2;
    static integer nr1, np4, nt4, nt6;
    static doublereal acc, arg, one, hti, htj, ten, eps;
    static integer jlt, npp;
    static doublereal piv;
    static integer num, nrr, ntt;
    static doublereal fac1, fac2;
    static integer ich1, ich3;
    static doublereal con1, con4, con9;
    static integer num1;
    static doublereal facc, half, facs;
    static integer ncof;
    static doublereal dmax__;
    static integer nreg, rank, iter;
    static doublereal fpms, pinv;
    static integer irot, jrot, iband, ncoff;
    static doublereal sigma, fpmax;
    static integer nrint;
    static doublereal store;
    static integer iband1, lwest, iband3, iband4;
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
	    doublereal *, doublereal *, doublereal *), fprpsp_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *);

/*  .. */
/*  ..scalar arguments.. */
/*<    >*/
/*<       real s,eta,tol,fp,sup >*/
/*  ..array arguments.. */
/*<    >*/
/*<       integer index(nrest),nummer(m) >*/
/*  ..local scalars.. */
/*<    >*/
/*<    >*/
/*  ..local arrays.. */
/*<       real ht(4),hp(4) >*/
/*  ..function references.. */
/*<       real abs,atan,fprati,sqrt,cos,sin >*/
/*<       integer min0 >*/
/*  ..subroutine references.. */
/*   fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota,fprpsp */
/*  .. */
/*  set constants */
/*<       one = 0.1e+01 >*/
    /* Parameter adjustments */
    --nummer;
    spp_dim1 = *m;
    spp_offset = 1 + spp_dim1;
    spp -= spp_offset;
    spt_dim1 = *m;
    spt_offset = 1 + spt_dim1;
    spt -= spt_offset;
    --w;
    --r__;
    --phi;
    --teta;
    bt_dim1 = *ntest;
    bt_offset = 1 + bt_dim1;
    bt -= bt_offset;
    --tt;
    bp_dim1 = *npest;
    bp_offset = 1 + bp_dim1;
    bp -= bp_offset;
    --cosi;
    --coco;
    --row;
    --tp;
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
/*<       con1 = 0.1e0 >*/
    con1 = .1;
/*<       con9 = 0.9e0 >*/
    con9 = .9;
/*<       con4 = 0.4e-01 >*/
    con4 = .04;
/*<       half = 0.5e0 >*/
    half = .5;
/*<       ten = 0.1e+02 >*/
    ten = 10.;
/*<       pi = atan(one)*4 >*/
    pi = atan(one) * 4;
/*<       pi2 = pi+pi >*/
    pi2 = pi + pi;
/*<       eps = sqrt(eta) >*/
    eps = sqrt(*eta);
/*<       if(iopt.lt.0) go to 70 >*/
    if (*iopt < 0) {
	goto L70;
    }
/*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
/*<       acc = tol*s >*/
    acc = *tol * *s;
/*<       if(iopt.eq.0) go to 10 >*/
    if (*iopt == 0) {
	goto L10;
    }
/*<       if(s.lt.sup) if(np-11) 60,70,70 >*/
    if (*s < *sup) {
	if (*np - 11 >= 0) {
	    goto L70;
	} else {
	    goto L60;
	}
    }
/*  if iopt=0 we begin by computing the weighted least-squares polynomial */
/*  of the form */
/*     s(teta,phi) = c1*f1(teta) + cn*fn(teta) */
/*  where f1(teta) and fn(teta) are the cubic polynomials satisfying */
/*     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta). */
/*  the corresponding weighted sum of squared residuals gives the upper */
/*  bound sup for the smoothing factor s. */
/*<   10  sup = 0. >*/
L10:
    *sup = 0.;
/*<       d1 = 0. >*/
    d1 = 0.;
/*<       d2 = 0. >*/
    d2 = 0.;
/*<       c1 = 0. >*/
    c1 = 0.;
/*<       cn = 0. >*/
    cn = 0.;
/*<       fac1 = pi*(one + half) >*/
    fac1 = pi * (one + half);
/*<       fac2 = (one + one)/pi**3 >*/
/* Computing 3rd power */
    d__1 = pi;
    fac2 = (one + one) / (d__1 * (d__1 * d__1));
/*<       aa = 0. >*/
    aa = 0.;
/*<       do 40 i=1,m >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          wi = w(i) >*/
	wi = w[i__];
/*<          ri = r(i)*wi >*/
	ri = r__[i__] * wi;
/*<          arg = teta(i) >*/
	arg = teta[i__];
/*<          fn = fac2*arg*arg*(fac1-arg) >*/
	fn = fac2 * arg * arg * (fac1 - arg);
/*<          f1 = (one-fn)*wi >*/
	f1 = (one - fn) * wi;
/*<          fn = fn*wi >*/
	fn *= wi;
/*<          if(fn.eq.0.) go to 20 >*/
	if (fn == 0.) {
	    goto L20;
	}
/*<          call fpgivs(fn,d1,co,si) >*/
	fpgivs_(&fn, &d1, &co, &si);
/*<          call fprota(co,si,f1,aa) >*/
	fprota_(&co, &si, &f1, &aa);
/*<          call fprota(co,si,ri,cn) >*/
	fprota_(&co, &si, &ri, &cn);
/*<  20      if(f1.eq.0.) go to 30 >*/
L20:
	if (f1 == 0.) {
	    goto L30;
	}
/*<          call fpgivs(f1,d2,co,si) >*/
	fpgivs_(&f1, &d2, &co, &si);
/*<          call fprota(co,si,ri,c1) >*/
	fprota_(&co, &si, &ri, &c1);
/*<  30      sup = sup+ri*ri >*/
L30:
	*sup += ri * ri;
/*<  40   continue >*/
/* L40: */
    }
/*<       if(d2.ne.0.) c1 = c1/d2 >*/
    if (d2 != 0.) {
	c1 /= d2;
    }
/*<       if(d1.ne.0.) cn = (cn-aa*c1)/d1 >*/
    if (d1 != 0.) {
	cn = (cn - aa * c1) / d1;
    }
/*  find the b-spline representation of this least-squares polynomial */
/*<       nt = 8 >*/
    *nt = 8;
/*<       np = 8 >*/
    *np = 8;
/*<       do 50 i=1,4 >*/
    for (i__ = 1; i__ <= 4; ++i__) {
/*<          c(i) = c1 >*/
	c__[i__] = c1;
/*<          c(i+4) = c1 >*/
	c__[i__ + 4] = c1;
/*<          c(i+8) = cn >*/
	c__[i__ + 8] = cn;
/*<          c(i+12) = cn >*/
	c__[i__ + 12] = cn;
/*<          tt(i) = 0. >*/
	tt[i__] = 0.;
/*<          tt(i+4) = pi >*/
	tt[i__ + 4] = pi;
/*<          tp(i) = 0. >*/
	tp[i__] = 0.;
/*<          tp(i+4) = pi2 >*/
	tp[i__ + 4] = pi2;
/*<   50  continue >*/
/* L50: */
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
/*<   60  if(npest.lt.11 .or. ntest.lt.9) go to 950 >*/
L60:
    if (*npest < 11 || *ntest < 9) {
	goto L950;
    }
/*  find the initial set of interior knots of the spherical spline in */
/*  case iopt = 0. */
/*<       np = 11 >*/
    *np = 11;
/*<       tp(5) = pi*half >*/
    tp[5] = pi * half;
/*<       tp(6) = pi >*/
    tp[6] = pi;
/*<       tp(7) = tp(5)+pi >*/
    tp[7] = tp[5] + pi;
/*<       nt = 9 >*/
    *nt = 9;
/*<       tt(5) = tp(5) >*/
    tt[5] = tp[5];
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  part 1 : computation of least-squares spherical splines.            c */
/*  ********************************************************            c */
/*  if iopt < 0 we compute the least-squares spherical spline according c */
/*  to the given set of knots.                                          c */
/*  if iopt >=0 we compute least-squares spherical splines with increas-c */
/*  ing numbers of knots until the corresponding sum f(p=inf)<=s.       c */
/*  the initial set of knots then depends on the value of iopt:         c */
/*    if iopt=0 we start with one interior knot in the teta-direction   c */
/*              (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). c */
/*    if iopt>0 we start with the set of knots found at the last call   c */
/*              of the routine.                                         c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  main loop for the different sets of knots. m is a save upper bound */
/*  for the number of trials. */
/*<   70  do 570 iter=1,m >*/
L70:
    i__1 = *m;
    for (iter = 1; iter <= i__1; ++iter) {
/*  find the position of the additional knots which are needed for the */
/*  b-spline representation of s(teta,phi). */
/*<          l1 = 4 >*/
	l1 = 4;
/*<          l2 = l1 >*/
	l2 = l1;
/*<          l3 = np-3 >*/
	l3 = *np - 3;
/*<          l4 = l3 >*/
	l4 = l3;
/*<          tp(l2) = 0. >*/
	tp[l2] = 0.;
/*<          tp(l3) = pi2 >*/
	tp[l3] = pi2;
/*<          do 80 i=1,3 >*/
	for (i__ = 1; i__ <= 3; ++i__) {
/*<             l1 = l1+1 >*/
	    ++l1;
/*<             l2 = l2-1 >*/
	    --l2;
/*<             l3 = l3+1 >*/
	    ++l3;
/*<             l4 = l4-1 >*/
	    --l4;
/*<             tp(l2) = tp(l4)-pi2 >*/
	    tp[l2] = tp[l4] - pi2;
/*<             tp(l3) = tp(l1)+pi2 >*/
	    tp[l3] = tp[l1] + pi2;
/*<   80     continue >*/
/* L80: */
	}
/*<         l = nt >*/
	l = *nt;
/*<         do 90 i=1,4 >*/
	for (i__ = 1; i__ <= 4; ++i__) {
/*<           tt(i) = 0. >*/
	    tt[i__] = 0.;
/*<           tt(l) = pi >*/
	    tt[l] = pi;
/*<           l = l-1 >*/
	    --l;
/*<   90    continue >*/
/* L90: */
	}
/*  find nrint, the total number of knot intervals and nreg, the number */
/*  of panels in which the approximation domain is subdivided by the */
/*  intersection of knots. */
/*<         ntt = nt-7 >*/
	ntt = *nt - 7;
/*<         npp = np-7 >*/
	npp = *np - 7;
/*<         nrr = npp/2 >*/
	nrr = npp / 2;
/*<         nr1 = nrr+1 >*/
	nr1 = nrr + 1;
/*<         nrint = ntt+npp >*/
	nrint = ntt + npp;
/*<         nreg = ntt*npp >*/
	nreg = ntt * npp;
/*  arrange the data points according to the panel they belong to. */
/*<         call fporde(teta,phi,m,3,3,tt,nt,tp,np,nummer,index,nreg) >*/
	fporde_(&teta[1], &phi[1], m, &c__3, &c__3, &tt[1], nt, &tp[1], np, &
		nummer[1], &index[1], &nreg);
/*  find the b-spline coefficients coco and cosi of the cubic spline */
/*  approximations sc(phi) and ss(phi) for cos(phi) and sin(phi). */
/*<         do 100 i=1,npp >*/
	i__2 = npp;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<            coco(i) = 0. >*/
	    coco[i__] = 0.;
/*<            cosi(i) = 0. >*/
	    cosi[i__] = 0.;
/*<            do 100 j=1,npp >*/
	    i__3 = npp;
	    for (j = 1; j <= i__3; ++j) {
/*<               a(i,j) = 0. >*/
		a[i__ + j * a_dim1] = 0.;
/*<  100    continue >*/
/* L100: */
	    }
	}
/*  the coefficients coco and cosi are obtained from the conditions */
/*  sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4. */
/*<         do 150 i=1,npp >*/
	i__3 = npp;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<            l2 = i+3 >*/
	    l2 = i__ + 3;
/*<            arg = tp(l2) >*/
	    arg = tp[l2];
/*<            call fpbspl(tp,np,3,arg,l2,hp) >*/
	    fpbspl_(&tp[1], np, &c__3, &arg, &l2, hp);
/*<            do 110 j=1,npp >*/
	    i__2 = npp;
	    for (j = 1; j <= i__2; ++j) {
/*<               row(j) = 0. >*/
		row[j] = 0.;
/*<  110       continue >*/
/* L110: */
	    }
/*<            ll = i >*/
	    ll = i__;
/*<            do 120 j=1,3 >*/
	    for (j = 1; j <= 3; ++j) {
/*<               if(ll.gt.npp) ll= 1 >*/
		if (ll > npp) {
		    ll = 1;
		}
/*<               row(ll) = row(ll)+hp(j) >*/
		row[ll] += hp[j - 1];
/*<               ll = ll+1 >*/
		++ll;
/*<  120       continue >*/
/* L120: */
	    }
/*<            facc = cos(arg) >*/
	    facc = cos(arg);
/*<            facs = sin(arg) >*/
	    facs = sin(arg);
/*<            do 140 j=1,npp >*/
	    i__2 = npp;
	    for (j = 1; j <= i__2; ++j) {
/*<               piv = row(j) >*/
		piv = row[j];
/*<               if(piv.eq.0.) go to 140 >*/
		if (piv == 0.) {
		    goto L140;
		}
/*<               call fpgivs(piv,a(j,1),co,si) >*/
		fpgivs_(&piv, &a[j + a_dim1], &co, &si);
/*<               call fprota(co,si,facc,coco(j)) >*/
		fprota_(&co, &si, &facc, &coco[j]);
/*<               call fprota(co,si,facs,cosi(j)) >*/
		fprota_(&co, &si, &facs, &cosi[j]);
/*<               if(j.eq.npp) go to 150 >*/
		if (j == npp) {
		    goto L150;
		}
/*<               j1 = j+1 >*/
		j1 = j + 1;
/*<               i2 = 1 >*/
		i2 = 1;
/*<               do 130 l=j1,npp >*/
		i__4 = npp;
		for (l = j1; l <= i__4; ++l) {
/*<                  i2 = i2+1 >*/
		    ++i2;
/*<                  call fprota(co,si,row(l),a(j,i2)) >*/
		    fprota_(&co, &si, &row[l], &a[j + i2 * a_dim1]);
/*<  130          continue >*/
/* L130: */
		}
/*<  140       continue >*/
L140:
		;
	    }
/*<  150    continue >*/
L150:
	    ;
	}
/*<         call fpback(a,coco,npp,npp,coco,ncc) >*/
	fpback_(&a[a_offset], &coco[1], &npp, &npp, &coco[1], ncc);
/*<         call fpback(a,cosi,npp,npp,cosi,ncc) >*/
	fpback_(&a[a_offset], &cosi[1], &npp, &npp, &cosi[1], ncc);
/*  find ncof, the dimension of the spherical spline and ncoff, the */
/*  number of coefficients in the standard b-spline representation. */
/*<         nt4 = nt-4 >*/
	nt4 = *nt - 4;
/*<         np4 = np-4 >*/
	np4 = *np - 4;
/*<         ncoff = nt4*np4 >*/
	ncoff = nt4 * np4;
/*<         ncof = 6+npp*(ntt-1) >*/
	ncof = npp * (ntt - 1) + 6;
/*  find the bandwidth of the observation matrix a. */
/*<         iband = 4*npp >*/
	iband = npp << 2;
/*<         if(ntt.eq.4) iband = 3*(npp+1) >*/
	if (ntt == 4) {
	    iband = (npp + 1) * 3;
	}
/*<         if(ntt.lt.4) iband = ncof >*/
	if (ntt < 4) {
	    iband = ncof;
	}
/*<         iband1 = iband-1 >*/
	iband1 = iband - 1;
/*  initialize the observation matrix a. */
/*<         do 160 i=1,ncof >*/
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           f(i) = 0. >*/
	    f[i__] = 0.;
/*<           do 160 j=1,iband >*/
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
/*<             a(i,j) = 0. >*/
		a[i__ + j * a_dim1] = 0.;
/*<  160    continue >*/
/* L160: */
	    }
	}
/*  initialize the sum of squared residuals. */
/*<         fp = 0. >*/
	*fp = 0.;
/*  fetch the data points in the new order. main loop for the */
/*  different panels. */
/*<         do 340 num=1,nreg >*/
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
/*  fix certain constants for the current panel; jrot records the column */
/*  number of the first non-zero element in a row of the observation */
/*  matrix according to a data point of the panel. */
/*<           num1 = num-1 >*/
	    num1 = num - 1;
/*<           lt = num1/npp >*/
	    lt = num1 / npp;
/*<           l1 = lt+4 >*/
	    l1 = lt + 4;
/*<           lp = num1-lt*npp+1 >*/
	    lp = num1 - lt * npp + 1;
/*<           l2 = lp+3 >*/
	    l2 = lp + 3;
/*<           lt = lt+1 >*/
	    ++lt;
/*<           jrot = 0 >*/
	    jrot = 0;
/*<           if(lt.gt.2) jrot = 3+(lt-3)*npp >*/
	    if (lt > 2) {
		jrot = (lt - 3) * npp + 3;
	    }
/*  test whether there are still data points in the current panel. */
/*<           in = index(num) >*/
	    in = index[num];
/*<  170      if(in.eq.0) go to 340 >*/
L170:
	    if (in == 0) {
		goto L340;
	    }
/*  fetch a new data point. */
/*<           wi = w(in) >*/
	    wi = w[in];
/*<           ri = r(in)*wi >*/
	    ri = r__[in] * wi;
/*  evaluate for the teta-direction, the 4 non-zero b-splines at teta(in) */
/*<           call fpbspl(tt,nt,3,teta(in),l1,ht) >*/
	    fpbspl_(&tt[1], nt, &c__3, &teta[in], &l1, ht);
/*  evaluate for the phi-direction, the 4 non-zero b-splines at phi(in) */
/*<           call fpbspl(tp,np,3,phi(in),l2,hp) >*/
	    fpbspl_(&tp[1], np, &c__3, &phi[in], &l2, hp);
/*  store the value of these b-splines in spt and spp resp. */
/*<           do 180 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<             spp(in,i) = hp(i) >*/
		spp[in + i__ * spp_dim1] = hp[i__ - 1];
/*<             spt(in,i) = ht(i) >*/
		spt[in + i__ * spt_dim1] = ht[i__ - 1];
/*<  180      continue >*/
/* L180: */
	    }
/*  initialize the new row of observation matrix. */
/*<           do 190 i=1,iband >*/
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             h(i) = 0. >*/
		h__[i__] = 0.;
/*<  190      continue >*/
/* L190: */
	    }
/*  calculate the non-zero elements of the new row by making the cross */
/*  products of the non-zero b-splines in teta- and phi-direction and */
/*  by taking into account the conditions of the spherical splines. */
/*<           do 200 i=1,npp >*/
	    i__3 = npp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<              row(i) = 0. >*/
		row[i__] = 0.;
/*<  200      continue >*/
/* L200: */
	    }
/*  take into account the condition (3) of the spherical splines. */
/*<           ll = lp >*/
	    ll = lp;
/*<           do 210 i=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<              if(ll.gt.npp) ll=1 >*/
		if (ll > npp) {
		    ll = 1;
		}
/*<              row(ll) = row(ll)+hp(i) >*/
		row[ll] += hp[i__ - 1];
/*<              ll = ll+1 >*/
		++ll;
/*<  210      continue >*/
/* L210: */
	    }
/*  take into account the other conditions of the spherical splines. */
/*<           if(lt.gt.2 .and. lt.lt.(ntt-1)) go to 230 >*/
	    if (lt > 2 && lt < ntt - 1) {
		goto L230;
	    }
/*<           facc = 0. >*/
	    facc = 0.;
/*<           facs = 0. >*/
	    facs = 0.;
/*<           do 220 i=1,npp >*/
	    i__3 = npp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<              facc = facc+row(i)*coco(i) >*/
		facc += row[i__] * coco[i__];
/*<              facs = facs+row(i)*cosi(i) >*/
		facs += row[i__] * cosi[i__];
/*<  220     continue >*/
/* L220: */
	    }
/*  fill in the non-zero elements of the new row. */
/*<  230     j1 = 0 >*/
L230:
	    j1 = 0;
/*<          do 280 j =1,4 >*/
	    for (j = 1; j <= 4; ++j) {
/*<             jlt = j+lt >*/
		jlt = j + lt;
/*<             htj = ht(j) >*/
		htj = ht[j - 1];
/*<             if(jlt.gt.2 .and. jlt.le.nt4) go to 240 >*/
		if (jlt > 2 && jlt <= nt4) {
		    goto L240;
		}
/*<             j1 = j1+1 >*/
		++j1;
/*<             h(j1) = h(j1)+htj >*/
		h__[j1] += htj;
/*<             go to 280 >*/
		goto L280;
/*<  240        if(jlt.eq.3 .or. jlt.eq.nt4) go to 260 >*/
L240:
		if (jlt == 3 || jlt == nt4) {
		    goto L260;
		}
/*<             do 250 i=1,npp >*/
		i__3 = npp;
		for (i__ = 1; i__ <= i__3; ++i__) {
/*<                j1 = j1+1 >*/
		    ++j1;
/*<                h(j1) = row(i)*htj >*/
		    h__[j1] = row[i__] * htj;
/*<  250        continue >*/
/* L250: */
		}
/*<             go to 280 >*/
		goto L280;
/*<  260        if(jlt.eq.3) go to 270 >*/
L260:
		if (jlt == 3) {
		    goto L270;
		}
/*<             h(j1+1) = facc*htj >*/
		h__[j1 + 1] = facc * htj;
/*<             h(j1+2) = facs*htj >*/
		h__[j1 + 2] = facs * htj;
/*<             h(j1+3) = htj >*/
		h__[j1 + 3] = htj;
/*<             j1 = j1+2 >*/
		j1 += 2;
/*<             go to 280 >*/
		goto L280;
/*<  270        h(1) = h(1)+htj >*/
L270:
		h__[1] += htj;
/*<             h(2) = facc*htj >*/
		h__[2] = facc * htj;
/*<             h(3) = facs*htj >*/
		h__[3] = facs * htj;
/*<             j1 = 3 >*/
		j1 = 3;
/*<  280      continue >*/
L280:
		;
	    }
/*<           do 290 i=1,iband >*/
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             h(i) = h(i)*wi >*/
		h__[i__] *= wi;
/*<  290      continue >*/
/* L290: */
	    }
/*  rotate the row into triangle by givens transformations. */
/*<           irot = jrot >*/
	    irot = jrot;
/*<           do 310 i=1,iband >*/
	    i__3 = iband;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             irot = irot+1 >*/
		++irot;
/*<             piv = h(i) >*/
		piv = h__[i__];
/*<             if(piv.eq.0.) go to 310 >*/
		if (piv == 0.) {
		    goto L310;
		}
/*  calculate the parameters of the givens transformation. */
/*<             call fpgivs(piv,a(irot,1),co,si) >*/
		fpgivs_(&piv, &a[irot + a_dim1], &co, &si);
/*  apply that transformation to the right hand side. */
/*<             call fprota(co,si,ri,f(irot)) >*/
		fprota_(&co, &si, &ri, &f[irot]);
/*<             if(i.eq.iband) go to 320 >*/
		if (i__ == iband) {
		    goto L320;
		}
/*  apply that transformation to the left hand side. */
/*<             i2 = 1 >*/
		i2 = 1;
/*<             i3 = i+1 >*/
		i3 = i__ + 1;
/*<             do 300 j=i3,iband >*/
		i__4 = iband;
		for (j = i3; j <= i__4; ++j) {
/*<               i2 = i2+1 >*/
		    ++i2;
/*<               call fprota(co,si,h(j),a(irot,i2)) >*/
		    fprota_(&co, &si, &h__[j], &a[irot + i2 * a_dim1]);
/*<  300        continue >*/
/* L300: */
		}
/*<  310      continue >*/
L310:
		;
	    }
/*  add the contribution of the row to the sum of squares of residual */
/*  right hand sides. */
/*<  320      fp = fp+ri**2 >*/
L320:
/* Computing 2nd power */
	    d__1 = ri;
	    *fp += d__1 * d__1;
/*  find the number of the next data point in the panel. */
/*<  330      in = nummer(in) >*/
/* L330: */
	    in = nummer[in];
/*<           go to 170 >*/
	    goto L170;
/*<  340    continue >*/
L340:
	    ;
	}
/*  find dmax, the maximum value for the diagonal elements in the reduced */
/*  triangle. */
/*<         dmax = 0. >*/
	dmax__ = 0.;
/*<         do 350 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(a(i,1).le.dmax) go to 350 >*/
	    if (a[i__ + a_dim1] <= dmax__) {
		goto L350;
	    }
/*<           dmax = a(i,1) >*/
	    dmax__ = a[i__ + a_dim1];
/*<  350    continue >*/
L350:
	    ;
	}
/*  check whether the observation matrix is rank deficient. */
/*<         sigma = eps*dmax >*/
	sigma = eps * dmax__;
/*<         do 360 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           if(a(i,1).le.sigma) go to 370 >*/
	    if (a[i__ + a_dim1] <= sigma) {
		goto L370;
	    }
/*<  360    continue >*/
/* L360: */
	}
/*  backward substitution in case of full rank. */
/*<         call fpback(a,f,ncof,iband,c,ncc) >*/
	fpback_(&a[a_offset], &f[1], &ncof, &iband, &c__[1], ncc);
/*<         rank = ncof >*/
	rank = ncof;
/*<         do 365 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           q(i,1) = a(i,1)/dmax >*/
	    q[i__ + q_dim1] = a[i__ + a_dim1] / dmax__;
/*<  365    continue >*/
/* L365: */
	}
/*<         go to 390 >*/
	goto L390;
/*  in case of rank deficiency, find the minimum norm solution. */
/*<  370    lwest = ncof*iband+ncof+iband >*/
L370:
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
/*<         do 380 i=1,ncof >*/
	i__2 = ncof;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           ff(i) = f(i) >*/
	    ff[i__] = f[i__];
/*<           do 380 j=1,iband >*/
	    i__3 = iband;
	    for (j = 1; j <= i__3; ++j) {
/*<             q(i,j) = a(i,j) >*/
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/*<  380    continue >*/
/* L380: */
	    }
	}
/*<    >*/
	fprank_(&q[q_offset], &ff[1], &ncof, &iband, ncc, &sigma, &c__[1], &
		sq, &rank, &wrk[la], &wrk[lf], &wrk[lh]);
/*<         do 385 i=1,ncof >*/
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           q(i,1) = q(i,1)/dmax >*/
	    q[i__ + q_dim1] /= dmax__;
/*<  385    continue >*/
/* L385: */
	}
/*  add to the sum of squared residuals, the contribution of reducing */
/*  the rank. */
/*<         fp = fp+sq >*/
	*fp += sq;
/*  find the coefficients in the standard b-spline representation of */
/*  the spherical spline. */
/*<  390    call fprpsp(nt,np,coco,cosi,c,ff,ncoff) >*/
L390:
	fprpsp_(nt, np, &coco[1], &cosi[1], &c__[1], &ff[1], &ncoff);
/*  test whether the least-squares spline is an acceptable solution. */
/*<         if(iopt.lt.0) if(fp) 970,970,980 >*/
	if (*iopt < 0) {
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
/*  test whether we cannot further increase the number of knots. */
/*<         if(ncof.gt.m) go to 935 >*/
	if (ncof > *m) {
	    goto L935;
	}
/*  search where to add a new knot. */
/*  find for each interval the sum of squared residuals fpint for the */
/*  data points having the coordinate belonging to that knot interval. */
/*  calculate also coord which is the same sum, weighted by the position */
/*  of the data points considered. */
/*<  440    do 450 i=1,nrint >*/
/* L440: */
	i__3 = nrint;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           fpint(i) = 0. >*/
	    fpint[i__] = 0.;
/*<           coord(i) = 0. >*/
	    coord[i__] = 0.;
/*<  450    continue >*/
/* L450: */
	}
/*<         do 490 num=1,nreg >*/
	i__3 = nreg;
	for (num = 1; num <= i__3; ++num) {
/*<           num1 = num-1 >*/
	    num1 = num - 1;
/*<           lt = num1/npp >*/
	    lt = num1 / npp;
/*<           l1 = lt+1 >*/
	    l1 = lt + 1;
/*<           lp = num1-lt*npp >*/
	    lp = num1 - lt * npp;
/*<           l2 = lp+1+ntt >*/
	    l2 = lp + 1 + ntt;
/*<           jrot = lt*np4+lp >*/
	    jrot = lt * np4 + lp;
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
/*<             hti = spt(in,i) >*/
		hti = spt[in + i__ * spt_dim1];
/*<             j1 = i1 >*/
		j1 = i1;
/*<             do 470 j=1,4 >*/
		for (j = 1; j <= 4; ++j) {
/*<               j1 = j1+1 >*/
		    ++j1;
/*<               store = store+hti*spp(in,j)*c(j1) >*/
		    store += hti * spp[in + j * spp_dim1] * c__[j1];
/*<  470        continue >*/
/* L470: */
		}
/*<             i1 = i1+np4 >*/
		i1 += np4;
/*<  480      continue >*/
/* L480: */
	    }
/*<           store = (w(in)*(r(in)-store))**2 >*/
/* Computing 2nd power */
	    d__1 = w[in] * (r__[in] - store);
	    store = d__1 * d__1;
/*<           fpint(l1) = fpint(l1)+store >*/
	    fpint[l1] += store;
/*<           coord(l1) = coord(l1)+store*teta(in) >*/
	    coord[l1] += store * teta[in];
/*<           fpint(l2) = fpint(l2)+store >*/
	    fpint[l2] += store;
/*<           coord(l2) = coord(l2)+store*phi(in) >*/
	    coord[l2] += store * phi[in];
/*<           in = nummer(in) >*/
	    in = nummer[in];
/*<           go to 460 >*/
	    goto L460;
/*<  490    continue >*/
L490:
	    ;
	}
/*  find the interval for which fpint is maximal on the condition that */
/*  there still can be added a knot. */
/*<         l1 = 1 >*/
	l1 = 1;
/*<         l2 = nrint >*/
	l2 = nrint;
/*<         if(ntest.lt.nt+1) l1=ntt+1 >*/
	if (*ntest < *nt + 1) {
	    l1 = ntt + 1;
	}
/*<         if(npest.lt.np+2) l2=ntt >*/
	if (*npest < *np + 2) {
	    l2 = ntt;
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
	i__3 = l2;
	for (i__ = l1; i__ <= i__3; ++i__) {
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
/*<         if(l.gt.ntt) go to 530 >*/
	if (l > ntt) {
	    goto L530;
	}
/*  addition in the teta-direction */
/*<         l4 = l+4 >*/
	l4 = l + 4;
/*<         fpint(l) = 0. >*/
	fpint[l] = 0.;
/*<         fac1 = tt(l4)-arg >*/
	fac1 = tt[l4] - arg;
/*<         fac2 = arg-tt(l4-1) >*/
	fac2 = arg - tt[l4 - 1];
/*<         if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500 >*/
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
/*<         j = nt >*/
	j = *nt;
/*<         do 520 i=l4,nt >*/
	i__3 = *nt;
	for (i__ = l4; i__ <= i__3; ++i__) {
/*<           tt(j+1) = tt(j) >*/
	    tt[j + 1] = tt[j];
/*<           j = j-1 >*/
	    --j;
/*<  520    continue >*/
/* L520: */
	}
/*<         tt(l4) = arg >*/
	tt[l4] = arg;
/*<         nt = nt+1 >*/
	++(*nt);
/*<         go to 570 >*/
	goto L570;
/*  addition in the phi-direction */
/*<  530    l4 = l+4-ntt >*/
L530:
	l4 = l + 4 - ntt;
/*<         if(arg.lt.pi) go to 540 >*/
	if (arg < pi) {
	    goto L540;
	}
/*<         arg = arg-pi >*/
	arg -= pi;
/*<         l4 = l4-nrr >*/
	l4 -= nrr;
/*<  540    fpint(l) = 0. >*/
L540:
	fpint[l] = 0.;
/*<         fac1 = tp(l4)-arg >*/
	fac1 = tp[l4] - arg;
/*<         fac2 = arg-tp(l4-1) >*/
	fac2 = arg - tp[l4 - 1];
/*<         if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500 >*/
	if (fac1 > ten * fac2 || fac2 > ten * fac1) {
	    goto L500;
	}
/*<         ll = nrr+4 >*/
	ll = nrr + 4;
/*<         j = ll >*/
	j = ll;
/*<         do 550 i=l4,ll >*/
	i__3 = ll;
	for (i__ = l4; i__ <= i__3; ++i__) {
/*<           tp(j+1) = tp(j) >*/
	    tp[j + 1] = tp[j];
/*<           j = j-1 >*/
	    --j;
/*<  550    continue >*/
/* L550: */
	}
/*<         tp(l4) = arg >*/
	tp[l4] = arg;
/*<         np = np+2 >*/
	*np += 2;
/*<         nrr = nrr+1 >*/
	++nrr;
/*<         do 560 i=5,ll >*/
	i__3 = ll;
	for (i__ = 5; i__ <= i__3; ++i__) {
/*<           j = i+nrr >*/
	    j = i__ + nrr;
/*<           tp(j) = tp(i)+pi >*/
	    tp[j] = tp[i__] + pi;
/*<  560    continue >*/
/* L560: */
	}
/*  restart the computations with the new set of knots. */
/*<  570  continue >*/
L570:
	;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* part 2: determination of the smoothing spherical spline.             c */
/* ********************************************************             c */
/* we have determined the number of knots and their position. we now    c */
/* compute the coefficients of the smoothing spline sp(teta,phi).       c */
/* the observation matrix a is extended by the rows of a matrix, expres-c */
/* sing that sp(teta,phi) must be a constant function in the variable   c */
/* phi and a cubic polynomial in the variable teta. the corresponding   c */
/* weights of these additional rows are set to 1/(p). iteratively       c */
/* we than have to determine the value of p such that f(p) = sum((w(i)* c */
/* (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              c */
/* we already know that the least-squares polynomial corresponds to p=0,c */
/* and that the least-squares spherical spline corresponds to p=infin.  c */
/* the iteration process makes use of rational interpolation. since f(p)c */
/* is a convex and strictly decreasing function of p, it can be approx- c */
/* imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c */
/* three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c */
/* f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c */
/* of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tt(l),l=5,...,nt-4. */
/*<  580  call fpdisc(tt,nt,5,bt,ntest) >*/
L580:
    fpdisc_(&tt[1], nt, &c__5, &bt[bt_offset], ntest);
/*  evaluate the discontinuity jumps of the 3-th order derivative of */
/*  the b-splines at the knots tp(l),l=5,...,np-4. */
/*<       call fpdisc(tp,np,5,bp,npest) >*/
    fpdisc_(&tp[1], np, &c__5, &bp[bp_offset], npest);
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
/*<       do 585 i=1,ncof >*/
    i__1 = ncof;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         p = p+a(i,1) >*/
	p += a[i__ + a_dim1];
/*<  585  continue >*/
/* L585: */
    }
/*<       rn = ncof >*/
    rn = (doublereal) ncof;
/*<       p = rn/p >*/
    p = rn / p;
/*  find the bandwidth of the extended observation matrix. */
/*<       iband4 = iband+3 >*/
    iband4 = iband + 3;
/*<       if(ntt.le.4) iband4 = ncof >*/
    if (ntt <= 4) {
	iband4 = ncof;
    }
/*<       iband3 = iband4 -1 >*/
    iband3 = iband4 - 1;
/*<       ich1 = 0 >*/
    ich1 = 0;
/*<       ich3 = 0 >*/
    ich3 = 0;
/*  iteration process to find the root of f(p)=s. */
/*<       do 920 iter=1,maxit >*/
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/*<         pinv = one/p >*/
	pinv = one / p;
/*  store the triangularized observation matrix into q. */
/*<         do 600 i=1,ncof >*/
	i__3 = ncof;
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<           ff(i) = f(i) >*/
	    ff[i__] = f[i__];
/*<           do 590 j=1,iband4 >*/
	    i__2 = iband4;
	    for (j = 1; j <= i__2; ++j) {
/*<             q(i,j) = 0. >*/
		q[i__ + j * q_dim1] = 0.;
/*<  590      continue >*/
/* L590: */
	    }
/*<           do 600 j=1,iband >*/
	    i__2 = iband;
	    for (j = 1; j <= i__2; ++j) {
/*<             q(i,j) = a(i,j) >*/
		q[i__ + j * q_dim1] = a[i__ + j * a_dim1];
/*<  600    continue >*/
/* L600: */
	    }
	}
/*  extend the observation matrix with the rows of a matrix, expressing */
/*  that for teta=cst. sp(teta,phi) must be a constant function. */
/*<         nt6 = nt-6 >*/
	nt6 = *nt - 6;
/*<         do 720 i=5,np4 >*/
	i__2 = np4;
	for (i__ = 5; i__ <= i__2; ++i__) {
/*<           ii = i-4 >*/
	    ii = i__ - 4;
/*<           do 610 l=1,npp >*/
	    i__3 = npp;
	    for (l = 1; l <= i__3; ++l) {
/*<              row(l) = 0. >*/
		row[l] = 0.;
/*<  610      continue >*/
/* L610: */
	    }
/*<           ll = ii >*/
	    ll = ii;
/*<           do 620  l=1,5 >*/
	    for (l = 1; l <= 5; ++l) {
/*<              if(ll.gt.npp) ll=1 >*/
		if (ll > npp) {
		    ll = 1;
		}
/*<              row(ll) = row(ll)+bp(ii,l) >*/
		row[ll] += bp[ii + l * bp_dim1];
/*<              ll = ll+1 >*/
		++ll;
/*<  620      continue >*/
/* L620: */
	    }
/*<           facc = 0. >*/
	    facc = 0.;
/*<           facs = 0. >*/
	    facs = 0.;
/*<           do 630 l=1,npp >*/
	    i__3 = npp;
	    for (l = 1; l <= i__3; ++l) {
/*<              facc = facc+row(l)*coco(l) >*/
		facc += row[l] * coco[l];
/*<              facs = facs+row(l)*cosi(l) >*/
		facs += row[l] * cosi[l];
/*<  630      continue >*/
/* L630: */
	    }
/*<           do 720 j=1,nt6 >*/
	    i__3 = nt6;
	    for (j = 1; j <= i__3; ++j) {
/*  initialize the new row. */
/*<             do 640 l=1,iband >*/
		i__4 = iband;
		for (l = 1; l <= i__4; ++l) {
/*<               h(l) = 0. >*/
		    h__[l] = 0.;
/*<  640        continue >*/
/* L640: */
		}
/*  fill in the non-zero elements of the row. jrot records the column */
/*  number of the first non-zero element in the row. */
/*<             jrot = 4+(j-2)*npp >*/
		jrot = (j - 2) * npp + 4;
/*<             if(j.gt.1 .and. j.lt.nt6) go to 650 >*/
		if (j > 1 && j < nt6) {
		    goto L650;
		}
/*<             h(1) = facc >*/
		h__[1] = facc;
/*<             h(2) = facs >*/
		h__[2] = facs;
/*<             if(j.eq.1) jrot = 2 >*/
		if (j == 1) {
		    jrot = 2;
		}
/*<             go to 670 >*/
		goto L670;
/*<  650        do 660 l=1,npp >*/
L650:
		i__4 = npp;
		for (l = 1; l <= i__4; ++l) {
/*<                h(l)=row(l) >*/
		    h__[l] = row[l];
/*<  660        continue >*/
/* L660: */
		}
/*<  670        do 675 l=1,iband >*/
L670:
		i__4 = iband;
		for (l = 1; l <= i__4; ++l) {
/*<                h(l) = h(l)*pinv >*/
		    h__[l] *= pinv;
/*<  675        continue >*/
/* L675: */
		}
/*<             ri = 0. >*/
		ri = 0.;
/*  rotate the new row into triangle by givens transformations. */
/*<             do 710 irot=jrot,ncof >*/
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
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
/*<               call fprota(co,si,ri,ff(irot)) >*/
		    fprota_(&co, &si, &ri, &ff[irot]);
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
/*  that for phi=cst. sp(teta,phi) must be a cubic polynomial. */
/*<         do 810 i=5,nt4 >*/
	i__3 = nt4;
	for (i__ = 5; i__ <= i__3; ++i__) {
/*<           ii = i-4 >*/
	    ii = i__ - 4;
/*<           do 810 j=1,npp >*/
	    i__2 = npp;
	    for (j = 1; j <= i__2; ++j) {
/*  initialize the new row */
/*<             do 730 l=1,iband4 >*/
		i__4 = iband4;
		for (l = 1; l <= i__4; ++l) {
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
/*<                il = ii+l >*/
		    il = ii + l;
/*<                ij = npp >*/
		    ij = npp;
/*<                if(il.ne.3 .and. il.ne.nt4) go to 750 >*/
		    if (il != 3 && il != nt4) {
			goto L750;
		    }
/*<                j1 = j1+3-j >*/
		    j1 = j1 + 3 - j;
/*<                j2 = j1-2 >*/
		    j2 = j1 - 2;
/*<                ij = 0 >*/
		    ij = 0;
/*<                if(il.ne.3) go to 740 >*/
		    if (il != 3) {
			goto L740;
		    }
/*<                j1 = 1 >*/
		    j1 = 1;
/*<                j2 = 2 >*/
		    j2 = 2;
/*<                ij = j+2 >*/
		    ij = j + 2;
/*<  740           h(j2) = bt(ii,l)*coco(j) >*/
L740:
		    h__[j2] = bt[ii + l * bt_dim1] * coco[j];
/*<                h(j2+1) = bt(ii,l)*cosi(j) >*/
		    h__[j2 + 1] = bt[ii + l * bt_dim1] * cosi[j];
/*<  750           h(j1) = h(j1)+bt(ii,l) >*/
L750:
		    h__[j1] += bt[ii + l * bt_dim1];
/*<                j1 = j1+ij >*/
		    j1 += ij;
/*<  760        continue >*/
/* L760: */
		}
/*<             do 765 l=1,iband4 >*/
		i__4 = iband4;
		for (l = 1; l <= i__4; ++l) {
/*<                h(l) = h(l)*pinv >*/
		    h__[l] *= pinv;
/*<  765        continue >*/
/* L765: */
		}
/*<             ri = 0. >*/
		ri = 0.;
/*<             jrot = 1 >*/
		jrot = 1;
/*<             if(ii.gt.2) jrot = 3+j+(ii-3)*npp >*/
		if (ii > 2) {
		    jrot = j + 3 + (ii - 3) * npp;
		}
/*  rotate the new row into triangle by givens transformations. */
/*<             do 800 irot=jrot,ncof >*/
		i__4 = ncof;
		for (irot = jrot; irot <= i__4; ++irot) {
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
/*<               call fprota(co,si,ri,ff(irot)) >*/
		    fprota_(&co, &si, &ri, &ff[irot]);
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
/*  the spherical spline. */
/*<         call fprpsp(nt,np,coco,cosi,c,ff,ncoff) >*/
	fprpsp_(nt, np, &coco[1], &cosi[1], &c__[1], &ff[1], &ncoff);
/*  compute f(p). */
/*<         fp = 0. >*/
	*fp = 0.;
/*<         do 890 num = 1,nreg >*/
	i__2 = nreg;
	for (num = 1; num <= i__2; ++num) {
/*<           num1 = num-1 >*/
	    num1 = num - 1;
/*<           lt = num1/npp >*/
	    lt = num1 / npp;
/*<           lp = num1-lt*npp >*/
	    lp = num1 - lt * npp;
/*<           jrot = lt*np4+lp >*/
	    jrot = lt * np4 + lp;
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
/*<             hti = spt(in,i) >*/
		hti = spt[in + i__ * spt_dim1];
/*<             j1 = i1 >*/
		j1 = i1;
/*<             do 870 j=1,4 >*/
		for (j = 1; j <= 4; ++j) {
/*<               j1 = j1+1 >*/
		    ++j1;
/*<               store = store+hti*spp(in,j)*c(j1) >*/
		    store += hti * spp[in + j * spp_dim1] * c__[j1];
/*<  870        continue >*/
/* L870: */
		}
/*<             i1 = i1+np4 >*/
		i1 += np4;
/*<  880      continue >*/
/* L880: */
	    }
/*<           fp = fp+(w(in)*(r(in)-store))**2 >*/
/* Computing 2nd power */
	    d__1 = w[in] * (r__[in] - store);
	    *fp += d__1 * d__1;
/*<           in = nummer(in) >*/
	    in = nummer[in];
/*<           go to 860 >*/
	    goto L860;
/*<  890    continue >*/
L890:
	    ;
	}
/*  test whether the approximation sp(teta,phi) is an acceptable solution */
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
} /* fpsphe_ */

