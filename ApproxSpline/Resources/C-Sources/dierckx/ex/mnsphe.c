/* mnsphe.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__9 = 9;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c             mnsphe : sphere test program                           cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0021latitude-longitude values of the data points (degrees)\002)";
    static char fmt_905[] = "(\0020\002,4(3x,\002teta   phi\002,3x))";
    static char fmt_910[] = "(8f6.0)";
    static char fmt_915[] = "(\002 \002,4(3x,f4.0,2x,f4.0,3x))";
    static char fmt_920[] = "(\0020least-squares spline approximation on the sphere.\002)";
    static char fmt_925[] = "(\0020smoothing spline on the sphere.\002)";
    static char fmt_930[] = "(\002 smoothing factor s=\002,f9.0)";
    static char fmt_935[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i3)";
    static char fmt_940[] = "(1x,\002total number of knots in the teta-direction =\002,i3)";
    static char fmt_945[] = "(1x,\002position of the knots \002)";
    static char fmt_950[] = "(5x,8f8.4)";
    static char fmt_955[] = "(1x,\002total number of knots in the phi-direction =\002,i3)";
    static char fmt_960[] = "(\0020b-spline coefficients \002)";
    static char fmt_965[] = "(5x,8f9.4)";
    static char fmt_970[] = "(\002      phi\002,9f7.3)";
    static char fmt_975[] = "(\002  teta\002)";
    static char fmt_980[] = "(\002 \002,f6.3,2x,9f7.3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal);
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real teta[192];
    static integer iopt, iwrk[300], kwrk, lwrk1, lwrk2;
    static real c__[300], f[81];
    static integer i__, j, l, m;
    static real p[9], r__[192], s, t[9], w[192], scale;
    static integer npest, i1, i2, l1, l2, ntest;
    static real ai;
    static integer nc;
    static real fp, pi;
    static integer is, np, nt;
    static real tp[30], tt[30];
    extern /* Subroutine */ int sphere_(integer *, integer *, real *, real *, real *, real *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, integer *, integer *, integer *), bispev_(real *, integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, integer *, integer *, integer *);
    static real pi2, pi4;
    extern doublereal testsp_(real *, real *);
    static integer ier;
    static real phi[192], scp, eps, sct;
    static integer npp, ntt;
    static real wrk1[12000], wrk2[72];

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___7 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___11 = { 0, 5, 0, fmt_910, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_980, 0 };


/*  set constants */
    pi4 = atan(1.f);
    pi = pi4 * 4;
    pi2 = pi + pi;
    scale = pi4 / 45.f;
/*  we fetch the number of data points. */
    m = 192;
/*  we fetch and print the latitude - longitude coordinates of the data */
/*  points (in degrees). */
    s_wsfe(&io___6);
    e_wsfe();
    s_wsfe(&io___7);
    e_wsfe();
    l2 = 0;
    for (i__ = 1; i__ <= 48; ++i__) {
	l1 = l2 + 1;
	l2 += 4;
	s_rsfe(&io___11);
	i__1 = l2;
	for (l = l1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&teta[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&phi[l - 1], (ftnlen)sizeof(real));
	}
	e_rsfe();
	s_wsfe(&io___15);
	i__1 = l2;
	for (l = l1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&teta[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&phi[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L10: */
    }
/*  we set up the weights, scale into radians the latitude-longitude */
/*  coordinates and calculate the function values. */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 1.f;
	teta[i__ - 1] *= scale;
	phi[i__ - 1] *= scale;
	if (teta[i__ - 1] > pi) {
	    teta[i__ - 1] = pi;
	}
	if (phi[i__ - 1] > pi2) {
	    phi[i__ - 1] = pi2;
	}
	r__[i__ - 1] = testsp_(&teta[i__ - 1], &phi[i__ - 1]);
/* L20: */
    }
/*  we set up the coordinates of the grid points for the evaluation of */
/*  the spline approximations. */
    sct = pi / 8;
    scp = pi2 / 8;
    for (i__ = 1; i__ <= 8; ++i__) {
	ai = (real) (i__ - 1);
	t[i__ - 1] = ai * sct;
	p[i__ - 1] = ai * scp;
/* L30: */
    }
    t[8] = pi;
    p[8] = pi2;
/* we set up the dimension information */
    ntest = 15;
    npest = 19;
    lwrk1 = 12000;
    lwrk2 = 72;
    kwrk = 300;
/*  we choose a value for eps */
    eps = 1e-6f;
/*  main loop for the different spline approximations */
    for (is = 1; is <= 4; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	}
/*  we start computing the least-squares constrained polynomial (large s) */
L110:
	iopt = 0;
	s = 500.f;
	goto L200;
/*  iopt = 1 from the second call on. */
L120:
	iopt = 1;
	s = 135.f;
	goto L200;
L130:
	s = 15.f;
	goto L200;
/*  a least-squares spherical spline with specified knots. */
L140:
	iopt = -1;
/*  we set up the number of knots. */
	nt = 11;
	np = 15;
/*  we set up the position of the interior knots of the spline. */
	ntt = nt - 8;
	i__1 = ntt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ai = (real) i__;
	    j = i__ + 4;
	    tt[j - 1] = ai * pi4;
/* L150: */
	}
	npp = np - 8;
	i__1 = npp;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ai = (real) i__;
	    j = i__ + 4;
	    tp[j - 1] = ai * pi4;
/* L160: */
	}
/*  determination of the spline approximation. */
L200:
	sphere_(&iopt, &m, teta, phi, r__, w, &s, &ntest, &npest, &eps, &nt, tt, &np, tp, c__, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
/*  printing of the fitting results. */
	if (iopt >= 0) {
	    goto L210;
	}
	s_wsfe(&io___45);
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___46);
	e_wsfe();
	s_wsfe(&io___47);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___49);
	do_fio(&c__1, (char *)&nt, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___50);
	e_wsfe();
	s_wsfe(&io___51);
	i__1 = nt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tt[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___52);
	do_fio(&c__1, (char *)&np, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___53);
	e_wsfe();
	s_wsfe(&io___54);
	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tp[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nc = (nt - 4) * (np - 4);
	s_wsfe(&io___56);
	e_wsfe();
	s_wsfe(&io___57);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  evaluation of the spline approximation. */
	bispev_(tt, &nt, tp, &np, c__, &c__3, &c__3, t, &c__9, p, &c__9, f, wrk2, &lwrk2, iwrk, &kwrk, &ier);
	s_wsfe(&io___59);
	for (i__ = 1; i__ <= 9; ++i__) {
	    do_fio(&c__1, (char *)&p[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___60);
	e_wsfe();
	i2 = 0;
	for (i__ = 1; i__ <= 9; ++i__) {
	    i1 = i2 + 1;
	    i2 += 9;
	    s_wsfe(&io___63);
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    i__1 = i2;
	    for (j = i1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&f[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/* L230: */
	}
/* L300: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

doublereal testsp_(real *v, real *u)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    static real cu, cv, su, sv, rad1, rad2, rad3;

/* function program testsp calculates the value of a test function for */
/* the sphere package. */
/* .. */
    cu = cos(*u);
    cv = cos(*v);
    su = sin(*u);
    sv = sin(*v);
/* Computing 2nd power */
    r__1 = cu * sv * .2f;
/* Computing 2nd power */
    r__2 = su * sv;
/* Computing 2nd power */
    r__3 = cv * .5f;
    rad1 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/* Computing 2nd power */
    r__1 = cu * sv;
/* Computing 2nd power */
    r__2 = su * sv * .5f;
/* Computing 2nd power */
    r__3 = cv * .2f;
    rad2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
/* Computing 2nd power */
    r__1 = cu * sv * .5f;
/* Computing 2nd power */
    r__2 = su * sv * .2f;
/* Computing 2nd power */
    r__3 = cv;
    rad3 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
    ret_val = 1.f / sqrt(rad1) + 1.f / sqrt(rad2) + 1.f / sqrt(rad3);
    return ret_val;
} /* testsp_ */

