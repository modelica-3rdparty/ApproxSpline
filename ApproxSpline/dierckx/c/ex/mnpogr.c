/* mnpogr.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static real c_b16 = 0.f;
static integer c__3 = 3;
static integer c__116 = 116;
static integer c__29 = 29;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                  mnpogr : pogrid test program                      cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  ..local scalars.. */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(10f8.3)";
    static char fmt_905[] = "(\0021data value (exact function value) at grid points\002)";
    static char fmt_910[] = "(\002 u(i),i=\002,3x,9(i1,7x))";
    static char fmt_915[] = "(\002 v(j),j=\002)";
    static char fmt_920[] = "(1x,i5,9(2x,f6.3))";
    static char fmt_925[] = "(7x,9(\002 (\002,f5.3,\002)\002))";
    static char fmt_930[] = "(\0020data value at (0,0) = \002,f7.3,5x,\002exact value = \002,f7.3)";
    static char fmt_935[] = "(\0020mean abs. error = \002,f9.3,5x,\002max. abs. error = \002,f9.3)";
    static char fmt_940[] = "(\0020least-squares spline\002)";
    static char fmt_945[] = "(\0020smoothing spline with s=\002,f7.2)";
    static char fmt_950[] = "(1x,\002order of continuity at the origin =\002,i3)";
    static char fmt_955[] = "(1x,\002vanishing partial derivatives at the origin\002)";
    static char fmt_960[] = "(1x,\002vanishing at the boundary of the disc\002)";
    static char fmt_965[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i3)";
    static char fmt_970[] = "(1x,\002total number of knots in the u-direction =\002,i3)";
    static char fmt_975[] = "(1x,\002position of the knots \002)";
    static char fmt_980[] = "(5x,8f9.4)";
    static char fmt_985[] = "(1x,\002total number of knots in the v-direction =\002,i3)";
    static char fmt_990[] = "(\0020b-spline coefficients \002)";
    static char fmt_995[] = "(\0020spline value (approximation error) at grid points\002)";
    static char fmt_1000[] = "(\0020spline value at (0,0) = \002,f7.3,5x,\002error = \002,f7.3)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    double atan2(doublereal, doublereal);
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void), s_wsfe(cilist *), e_wsfe(void);
    double cos(doublereal), sin(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer ider[2], iopt[3], iwrk[100], kwrk, lwrk;
    static real c__[300], f[180];
    static integer i__, j, k, m;
    static real r__, s, u[9], v[20], x, y, z__[180], exact[180], ermax;
    static integer nuest, nvest;
    static real z0, ai;
    static integer nc;
    static real fp, cv, pi;
    static integer is, iw[29], mu, mv, nu, nv;
    static real wk[116], sp[9], sv, tu[50], tv[50];
    extern /* Subroutine */ int pogrid_(integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, real *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, integer *, integer *, integer *, integer *), bispev_(real *, integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, integer *, integer *, integer *);
    extern doublereal tespog_(real *, real *);
    static real er0, del;
    static integer ier;
    static real one, err[9], wrk[1600], sum, exz0;

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 5, 0, fmt_900, 0 };
    static cilist io___15 = { 0, 5, 0, fmt_900, 0 };
    static cilist io___17 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_990, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_995, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___72 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___77 = { 0, 6, 0, fmt_935, 0 };


/*  ..local arrays.. */
/*  ..function references.. */
/*  .. */
/*  set constants */
    one = 1.f;
    pi = atan2(0.f, -one);
/* we set up the radius of the disc */
    r__ = one;
/* we set up the number of u (radius)-values of the grid. */
    mu = 9;
/* we set up the u-coordinates of the grid. */
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) i__;
	u[i__ - 1] = ai * .1f;
/* L10: */
    }
/* we set up the number of v (angle)-values of the grid */
    mv = 20;
/* we set up the v-coordinates of the grid. */
    del = pi * .1f;
    i__1 = mv;
    for (j = 1; j <= i__1; ++j) {
	ai = (real) (j - 1);
	v[j - 1] = ai * del - pi;
/* L20: */
    }
/* we fetch the data values at the grid points. */
    m = mu * mv;
    s_rsfe(&io___13);
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&z__[i__ - 1], (ftnlen)sizeof(real));
    }
    e_rsfe();
/* we fetch the data value at the origin. */
    s_rsfe(&io___15);
    do_fio(&c__1, (char *)&z0, (ftnlen)sizeof(real));
    e_rsfe();
/* we print the data values at the grid points. we also compute and print */
/* the exact value of the test function underlying the data. */
    s_wsfe(&io___17);
    e_wsfe();
    s_wsfe(&io___18);
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
    }
    e_wsfe();
    s_wsfe(&io___19);
    e_wsfe();
    exz0 = tespog_(&c_b16, &c_b16);
    er0 = (r__1 = exz0 - z0, dabs(r__1));
    ermax = er0;
    sum = er0;
    i__1 = mv;
    for (j = 1; j <= i__1; ++j) {
	cv = cos(v[j - 1]);
	sv = sin(v[j - 1]);
	k = j;
	i__2 = mu;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = u[i__ - 1] * cv;
	    y = u[i__ - 1] * sv;
	    exact[k - 1] = tespog_(&x, &y);
	    err[i__ - 1] = (r__1 = exact[k - 1] - z__[k - 1], dabs(r__1));
	    sum += err[i__ - 1];
	    if (err[i__ - 1] > ermax) {
		ermax = err[i__ - 1];
	    }
	    k += mv;
/* L30: */
	}
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	i__2 = m;
	i__3 = mv;
	for (k = j; i__3 < 0 ? k >= i__2 : k <= i__2; k += i__3) {
	    do_fio(&c__1, (char *)&z__[k - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___32);
	i__3 = m;
	i__2 = mv;
	for (k = j; i__2 < 0 ? k >= i__3 : k <= i__3; k += i__2) {
	    do_fio(&c__1, (char *)&exact[k - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L40: */
    }
    ai = (real) (m + 1);
    sum /= ai;
    s_wsfe(&io___33);
    do_fio(&c__1, (char *)&z0, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&exz0, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___34);
    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&ermax, (ftnlen)sizeof(real));
    e_wsfe();
/*  we set up the dimension information */
    nuest = 16;
    nvest = 27;
    kwrk = 100;
    lwrk = 1600;
/* main loop for the different spline approximations */
    for (is = 1; is <= 6; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L150;
	    case 6:  goto L160;
	}
/*  we start computing a set of spline approximations with */
/*  only c0-continuity at the origin, */
L110:
	iopt[1] = 0;
	ider[1] = 0;
/*  non-vanishing at the boundary of the disc, */
	iopt[2] = 0;
/*  with a data value at the origin. */
	ider[0] = 0;
/*  initialisation */
	iopt[0] = 0;
/*  a large value for s for computing the least-squares polynomial */
	s = 5.f;
	goto L200;
/*  iopt(1) = 1 from the second call on */
L120:
	s = .1f;
	iopt[0] = 1;
	goto L200;
/*  an interpolating spline */
L130:
	s = 0.f;
	goto L200;
/*  a second set of approximations with c1-continuity at the origin */
L140:
	iopt[1] = 1;
/*  vanishing at the boundary of the disc. */
	iopt[2] = 1;
/*  exact value at the origin. */
	ider[0] = 1;
	z0 = exz0;
/* reinitialization */
	iopt[0] = 0;
	s = .1f;
	goto L200;
/*  no data value at the origin */
L150:
	ider[0] = -1;
/*  vanishing partial derivatives at the origin */
	ider[1] = 1;
/* reinitialization */
	iopt[0] = 0;
	goto L200;
/* finally we calculate the least-squares spline according to the current */
/*  set of knots */
L160:
	iopt[0] = -1;
L200:
	pogrid_(iopt, ider, &mu, u, &mv, v, z__, &z0, &r__, &s, &nuest, &nvest, &nu, tu, &nv, tv, c__, &fp, wrk, &lwrk, iwrk, &kwrk, &ier);
/* printing of the fitting results. */
	if (iopt[0] >= 0) {
	    goto L210;
	}
	s_wsfe(&io___52);
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___54);
	do_fio(&c__1, (char *)&iopt[1], (ftnlen)sizeof(integer));
	e_wsfe();
	if (ider[1] == 1) {
	    s_wsfe(&io___55);
	    e_wsfe();
	}
	if (iopt[2] == 1) {
	    s_wsfe(&io___56);
	    e_wsfe();
	}
	s_wsfe(&io___57);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___58);
	do_fio(&c__1, (char *)&nu, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___59);
	e_wsfe();
	s_wsfe(&io___60);
	i__1 = nu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tu[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___61);
	do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___62);
	e_wsfe();
	s_wsfe(&io___63);
	i__1 = nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tv[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nc = (nu - 4) * (nv - 4);
	s_wsfe(&io___65);
	e_wsfe();
	s_wsfe(&io___66);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  evaluation of the spline approximation */
	bispev_(tu, &nu, tv, &nv, c__, &c__3, &c__3, u, &mu, v, &mv, f, wk, &c__116, iw, &c__29, &ier);
	s_wsfe(&io___70);
	e_wsfe();
	s_wsfe(&io___71);
	i__1 = mu;
	for (i__ = 1; i__ <= i__1; i__ += 2) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	}
	e_wsfe();
	s_wsfe(&io___72);
	e_wsfe();
	er0 = (r__1 = exz0 - c__[0], dabs(r__1));
	ermax = er0;
	sum = er0;
	i__1 = mv;
	for (j = 1; j <= i__1; ++j) {
	    k = j;
	    i__2 = mu;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sp[i__ - 1] = f[k - 1];
		err[i__ - 1] = (r__1 = exact[k - 1] - f[k - 1], dabs(r__1));
		sum += err[i__ - 1];
		if (err[i__ - 1] > ermax) {
		    ermax = err[i__ - 1];
		}
		k += mv;
/* L230: */
	    }
	    if (j / 3 * 3 != j) {
		goto L240;
	    }
	    s_wsfe(&io___74);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    i__2 = mu;
	    for (i__ = 1; i__ <= i__2; i__ += 2) {
		do_fio(&c__1, (char *)&sp[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___75);
	    i__2 = mu;
	    for (i__ = 1; i__ <= i__2; i__ += 2) {
		do_fio(&c__1, (char *)&err[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
L240:
	    ;
	}
	sum /= ai;
	s_wsfe(&io___76);
	do_fio(&c__1, (char *)&c__[0], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&er0, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___77);
	do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ermax, (ftnlen)sizeof(real));
	e_wsfe();
/* L300: */
    }
    s_stop("", 0L);
    return 0;
} /* MAIN__ */


doublereal tespog_(real *x, real *y)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real f;

/* function program tespog calculates the value of the test function */
/* underlying the data. */
/*  .. */
/*  ..scalar arguments.. */
/*  .. */
/* Computing 2nd power */
    r__1 = *x * 3.f - 1.f;
/* Computing 2nd power */
    r__2 = *y * 3.f - 1.f;
    f = 1.f - (r__1 * r__1 + r__2 * r__2) / (11.f - (*x + *y) * 6.f);
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = *y;
    ret_val = f - (1.f - r__1 * r__1 - r__2 * r__2) * (*x + *y) * 54.f / 121.f;
    return ret_val;
} /* tespog_ */

