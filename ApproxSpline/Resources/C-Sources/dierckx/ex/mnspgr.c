/* mnspgr.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static real c_b21 = 0.f;
static integer c__3 = 3;
static integer c__100 = 100;
static integer c__25 = 25;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                  mnspgr : spgrid test program                      cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*  ..local scalars.. */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_905[] = "(14i3)";
    static char fmt_900[] = "(7f8.3)";
    static char fmt_910[] = "(\0021data value (exact function value) at grid points\002)";
    static char fmt_915[] = "(\0020v(j),j=\002,3x,7(i2,6x))";
    static char fmt_920[] = "(\002 u(i),i=\002)";
    static char fmt_925[] = "(1x,i5,7(2x,f6.3))";
    static char fmt_930[] = "(7x,7(\002 (\002,f5.3,\002)\002))";
    static char fmt_935[] = "(\0020mean abs. error = \002,f9.3,5x,\002max. abs. error = \002,f9.3)";
    static char fmt_940[] = "(\002 function values at the poles \002,f7.3,5x,f7.3)";
    static char fmt_945[] = "(\0020least-squares spline\002)";
    static char fmt_950[] = "(\0020smoothing spline with s=\002,f7.2)";
    static char fmt_955[] = "(1x,\002order of continuity at the poles =\002,2i5)";
    static char fmt_960[] = "(1x,\002vanishing derivatives at the pole u=0\002)";
    static char fmt_965[] = "(1x,\002vanishing derivatives at the pole u=pi\002)";
    static char fmt_970[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i3)";
    static char fmt_975[] = "(1x,\002total number of knots in the u-direction =\002,i3)";
    static char fmt_980[] = "(1x,\002position of the knots \002)";
    static char fmt_985[] = "(5x,8f9.4)";
    static char fmt_990[] = "(1x,\002total number of knots in the v-direction =\002,i3)";
    static char fmt_995[] = "(\0020b-spline coefficients \002)";
    static char fmt_1000[] = "(\0020spline value (approximation error) at grid points\002)";
    static char fmt_1005[] = "(\002 spline values at the poles \002,f7.3,5x,f7.3)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double atan2(doublereal, doublereal);
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void), s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer ider[4], iopt[3], iwrk[70], kwrk, lwrk;
    static real c__[300], f[154];
    static integer i__, j, k, l, m;
    static real r__[154], s, u[11], v[14], exact[154], ermax;
    static integer nuest, nvest;
    static real r0, r1, ai;
    static integer nc;
    static real fp, pi;
    static integer is, iw[25], mu, mv, nu, nv;
    static real wk[100], sp[14], tu[25], tv[25];
    extern /* Subroutine */ int spgrid_(integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, real *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, integer *, integer *, integer *, integer *), bispev_(real *, integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, integer *, integer *, integer *);
    extern doublereal tesspg_(real *, real *);
    static real del, erf;
    static integer ier;
    static real one, err[14], wrk[1500], sum, exr0, exr1;

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 5, 0, fmt_905, 0 };
    static cilist io___11 = { 0, 5, 0, fmt_905, 0 };
    static cilist io___14 = { 0, 5, 0, fmt_900, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___17 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_990, 0 };
    static cilist io___67 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___70 = { 0, 6, 0, fmt_995, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___74 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___77 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___78 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___79 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___80 = { 0, 6, 0, fmt_1005, 0 };


/*  ..local arrays.. */
/*  ..function references.. */
/*  .. */
/*  set constants */
    one = 1.f;
    pi = atan2(0.f, -one);
    del = pi * .05f;
/* we set up the number of u (latitude)-values of the grid. */
    mu = 11;
/* we set up the u-coordinates of the grid. */
    s_rsfe(&io___5);
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&iw[i__ - 1], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) iw[i__ - 1];
	u[i__ - 1] = ai * del;
/* L10: */
    }
/* we set up the number of v (longitude)-values of the grid */
    mv = 14;
/* we set up the v-coordinates of the grid. */
    s_rsfe(&io___11);
    i__1 = mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&iw[i__ - 1], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    i__1 = mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) iw[i__ - 1];
	v[i__ - 1] = ai * del;
/* L20: */
    }
/* we fetch the data values at the grid points. */
    m = mu * mv;
    s_rsfe(&io___14);
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&r__[i__ - 1], (ftnlen)sizeof(real));
    }
    e_rsfe();
/* we print the data values at the grid points. we also compute and print */
/* the exact value of the test function underlying the data. */
    s_wsfe(&io___16);
    e_wsfe();
    s_wsfe(&io___17);
    i__1 = mv;
    for (j = 1; j <= i__1; j += 2) {
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    }
    e_wsfe();
    s_wsfe(&io___19);
    e_wsfe();
    exr0 = tesspg_(&c_b21, &c_b21);
    exr1 = tesspg_(&pi, &c_b21);
    ermax = 0.f;
    sum = 0.f;
    l = 0;
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = (i__ - 1) * mv + 1;
	k = 1;
	for (j = 1; j <= 7; ++j) {
	    exact[l - 1] = tesspg_(&u[i__ - 1], &v[k - 1]);
	    erf = (r__1 = exact[l - 1] - r__[l - 1], dabs(r__1));
	    sum += erf;
	    if (erf > ermax) {
		ermax = erf;
	    }
	    sp[j - 1] = exact[l - 1];
	    err[j - 1] = r__[l - 1];
	    l += 2;
	    k += 2;
/* L30: */
	}
	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	for (j = 1; j <= 7; ++j) {
	    do_fio(&c__1, (char *)&err[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___31);
	for (j = 1; j <= 7; ++j) {
	    do_fio(&c__1, (char *)&sp[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L40: */
    }
    s_wsfe(&io___32);
    i__1 = mv;
    for (j = 2; j <= i__1; j += 2) {
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    }
    e_wsfe();
    s_wsfe(&io___33);
    e_wsfe();
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = (i__ - 1) * mv + 2;
	k = 2;
	for (j = 1; j <= 7; ++j) {
	    exact[l - 1] = tesspg_(&u[i__ - 1], &v[k - 1]);
	    erf = (r__1 = exact[l - 1] - r__[l - 1], dabs(r__1));
	    sum += erf;
	    if (erf > ermax) {
		ermax = erf;
	    }
	    sp[j - 1] = exact[l - 1];
	    err[j - 1] = r__[l - 1];
	    l += 2;
	    k += 2;
/* L50: */
	}
	s_wsfe(&io___34);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	for (j = 1; j <= 7; ++j) {
	    do_fio(&c__1, (char *)&err[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___35);
	for (j = 1; j <= 7; ++j) {
	    do_fio(&c__1, (char *)&sp[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L60: */
    }
    ai = (real) m;
    sum /= ai;
    s_wsfe(&io___36);
    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&ermax, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___37);
    do_fio(&c__1, (char *)&exr0, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&exr1, (ftnlen)sizeof(real));
    e_wsfe();
/*  we set up the dimension information */
    nuest = 19;
    nvest = 21;
    kwrk = 70;
    lwrk = 1500;
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
/*  only c0-continuity at the poles, */
L110:
	iopt[1] = 0;
	ider[1] = 0;
	iopt[2] = 0;
	ider[3] = 0;
/*  with no data values at the poles. */
	ider[0] = -1;
	ider[2] = -1;
/*  initialisation */
	iopt[0] = 0;
/*  a large value for s for computing the least-squares polynomial */
	s = 60.f;
	goto L200;
/*  iopt(1) = 1 from the second call on */
L120:
	s = .05f;
	iopt[0] = 1;
	goto L200;
/*  an interpolating spline */
L130:
	s = 0.f;
	goto L200;
/*  a second set of approximations with c1-continuity at the poles */
L140:
	iopt[1] = 1;
	iopt[2] = 1;
/*  exact values at the poles. */
	ider[0] = 1;
	ider[2] = 1;
	r0 = exr0;
	r1 = exr1;
/* reinitialization */
	iopt[0] = 0;
	s = .05f;
	goto L200;
/*  vanishing derivatives at the poles */
L150:
	ider[1] = 1;
	ider[3] = 1;
/* reinitialization */
	iopt[0] = 0;
	goto L200;
/* finally we calculate the least-squares spline according to the current */
/*  set of knots */
L160:
	iopt[0] = -1;
L200:
	spgrid_(iopt, ider, &mu, u, &mv, v, r__, &r0, &r1, &s, &nuest, &nvest, &nu, tu, &nv, tv, c__, &fp, wrk, &lwrk, iwrk, &kwrk, &ier);
/* printing of the fitting results. */
	if (iopt[0] >= 0) {
	    goto L210;
	}
	s_wsfe(&io___57);
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___58);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___59);
	do_fio(&c__1, (char *)&iopt[1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iopt[2], (ftnlen)sizeof(integer));
	e_wsfe();
	if (ider[1] == 1) {
	    s_wsfe(&io___60);
	    e_wsfe();
	}
	if (ider[3] == 1) {
	    s_wsfe(&io___61);
	    e_wsfe();
	}
	s_wsfe(&io___62);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___63);
	do_fio(&c__1, (char *)&nu, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___64);
	e_wsfe();
	s_wsfe(&io___65);
	i__1 = nu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tu[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___66);
	do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___67);
	e_wsfe();
	s_wsfe(&io___68);
	i__1 = nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tv[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nc = (nu - 4) * (nv - 4);
	s_wsfe(&io___70);
	e_wsfe();
	s_wsfe(&io___71);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  evaluation of the spline approximation */
	bispev_(tu, &nu, tv, &nv, c__, &c__3, &c__3, u, &mu, v, &mv, f, wk, &c__100, iw, &c__25, &ier);
	s_wsfe(&io___74);
	e_wsfe();
	s_wsfe(&io___75);
	i__1 = mv;
	for (j = 1; j <= i__1; j += 2) {
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	}
	e_wsfe();
	s_wsfe(&io___76);
	e_wsfe();
	ermax = 0.f;
	sum = 0.f;
	i__1 = mu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = i__;
	    i__2 = mv;
	    for (j = 1; j <= i__2; ++j) {
		sp[j - 1] = f[k - 1];
		err[j - 1] = (r__1 = exact[k - 1] - f[k - 1], dabs(r__1));
		sum += err[j - 1];
		if (err[j - 1] > ermax) {
		    ermax = err[j - 1];
		}
		++k;
/* L230: */
	    }
	    if (i__ / 2 << 1 != i__) {
		goto L240;
	    }
	    s_wsfe(&io___77);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    i__2 = mv;
	    for (j = 1; j <= i__2; j += 2) {
		do_fio(&c__1, (char *)&sp[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___78);
	    i__2 = mv;
	    for (j = 1; j <= i__2; j += 2) {
		do_fio(&c__1, (char *)&err[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
L240:
	    ;
	}
	sum /= ai;
	s_wsfe(&io___79);
	do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ermax, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___80);
	do_fio(&c__1, (char *)&c__[0], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&c__[nc - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L300: */
    }
    s_stop("", 0L);
    return 0;
} /* MAIN__ */


doublereal tesspg_(real *u, real *v)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

/* function program tesspg calculates the value of the test function */
/* underlying the data. */
/* Computing 2nd power */
    r__1 = sin(*u);
    ret_val = 2.f / (cos(*u * 3.f) + 4.1f + cos(*v + *v + *u * .25f) * 3.f * (r__1 * r__1));
    return ret_val;
} /* tesspg_ */

