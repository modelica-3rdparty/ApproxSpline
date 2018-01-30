/* mnpola.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mnpola : polar test program                        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0021the input data\002)";
    static char fmt_905[] = "(\0020\002,3(3x,\002x\002,6x,\002y\002,6x,\002z\002,5x))";
    static char fmt_910[] = "(12f6.3)";
    static char fmt_915[] = "(\002 \002,3(3f7.3,2x))";
    static char fmt_920[] = "(\0020mean error = \002,f7.4,5x,\002max. error = \002,f7.4)";
    static char fmt_925[] = "(\0020smoothing spline on the disk with s =\002,f5.0)";
    static char fmt_930[] = "(\0020smoothing spline on the ellips with s =\002,f5.0)";
    static char fmt_935[] = "(\0020least-squares spline on the ellips\002)";
    static char fmt_940[] = "(\0020sum of squared residuals =\002,e15.6,5x,\002error flag =\002,i5)";
    static char fmt_945[] = "(1x,\002total number of knots in the u-direction =\002,i3)";
    static char fmt_950[] = "(1x,\002position of the knots \002)";
    static char fmt_955[] = "(5x,8f8.4)";
    static char fmt_960[] = "(1x,\002total number of knots in the v-direction =\002,i3)";
    static char fmt_965[] = "(\0020b-spline coefficients \002)";
    static char fmt_970[] = "(5x,8f9.4)";
    static char fmt_975[] = "(\0020spline values at selected points\002)";
    static char fmt_980[] = "(\0020\002,3(3x,\002x\002,6x,\002y\002,6x,\002f\002,5x))";

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer iopt[3], iwrk[500], kwrk, lwrk1, lwrk2;
    static real c__[300], f[200];
    static integer i__, l, m;
    static real s, u[200], v[200], w[200], x[200], y[200], z__[200], exact[200], ermax;
    extern /* Subroutine */ int polar_(integer *, integer *, real *, real *, real *, real *, E_fp, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, real *, real *, real *, integer *, real *, integer *, integer *, integer *, integer *);
    static real error;
    static integer l1, l2, m1, m2, nuest, nvest;
    static real ai;
    static integer nc;
    static real fp;
    static integer is, nu, nv;
    static real tu[30], tv[30];
    extern doublereal evapol_(real *, integer *, real *, integer *, real *, E_fp, real *, real *), testpo_(real *, real *);
    static integer ier;
    static real eps, sum;
    extern doublereal rad1_(), rad2_();
    static real wrk1[15000], wrk2[5700];

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___4 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___8 = { 0, 5, 0, fmt_910, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_920, 0 };


/*  we fetch the number of data points. */
    m1 = 200;
    m2 = 90;
/*  we fetch and print the coordinates and function values of the data. */
    s_wsfe(&io___3);
    e_wsfe();
    s_wsfe(&io___4);
    e_wsfe();
    l2 = 0;
    for (i__ = 1; i__ <= 50; ++i__) {
	l1 = l2 + 1;
	l2 += 4;
	s_rsfe(&io___8);
	i__1 = l2;
	for (l = l1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&x[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&y[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&z__[l - 1], (ftnlen)sizeof(real));
	}
	e_rsfe();
/* L10: */
    }
    s_wsfe(&io___13);
    i__1 = m1;
    for (l = 1; l <= i__1; ++l) {
	do_fio(&c__1, (char *)&x[l - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[l - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z__[l - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
/*  we calculate the exact function values and set up the weights w(i)= */
/*  (0.01)**(-1) (0.01 is an estimate for the standard deviation of the */
/*  error in z(i)). at the same time we calculate the mean and maximum */
/*  errors for the data values. */
    sum = 0.f;
    ermax = 0.f;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 100.f;
	exact[i__ - 1] = testpo_(&x[i__ - 1], &y[i__ - 1]);
	error = (r__1 = z__[i__ - 1] - exact[i__ - 1], dabs(r__1));
	sum += error;
	if (error > ermax) {
	    ermax = error;
	}
/* L20: */
    }
    ai = (real) m1;
    sum /= ai;
    s_wsfe(&io___20);
    do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&ermax, (ftnlen)sizeof(real));
    e_wsfe();
/*  we set up the dimension information */
    nuest = 15;
    nvest = 19;
    lwrk1 = 15000;
    lwrk2 = 5700;
    kwrk = 500;
/*  we choose a value for eps */
    eps = 1e-6f;
/*  main loop for the different spline approximations */
    for (is = 1; is <= 5; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L160;
	}
/*  we determine a number of smoothing spline approximations on the unit */
/*  disk x**2+y**2 <= 1. */
/*  all the data points are considered. */
L110:
	m = m1;
/*  we set up the smoothing factor. */
	s = 1500.f;
/*  the approximations are not restricted at the boundaries of the disk */
	iopt[2] = 0;
/*  we request c2-continuity at the origin. */
	iopt[1] = 2;
/*  at the first call of polar iopt(1) must be zero. */
	iopt[0] = 0;
	goto L200;
/*  iopt(1) = 1 from the second call on */
L120:
	iopt[0] = 1;
	s = 200.f;
	goto L200;
L130:
	s = 170.f;
	goto L200;
/*  we determine a smoothing spline approximation on the ellips */
/*  3*x**2+3*y**2-4*x*y<=1. */
/*  we only consider the data points inside this domain. */
L140:
	m = m2;
	ai = (real) m;
/*  the given function has a constant value 0.4 at the boundary of the */
/*  ellips. we calculate new data values by substracting this constant */
/*  from the old ones. */
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[i__ - 1] += -.4f;
/* L150: */
	}
/*  given these data we will then determine approximations which are */
/*  identically zero at the boundary of the ellips. */
	iopt[2] = 1;
/*  we still request c2-continuity at the origin. */
	iopt[1] = 2;
/*  reinitialization for the knots. */
	iopt[0] = 0;
/*  we set up the smoothing factor. */
	s = 90.f;
	goto L250;
/*  at the last call we will determine the least-squares spline */
/*  approximation corresponding to the current set of knots */
L160:
	iopt[0] = -1;
	goto L250;
/*  determination of the spline approximation on the disk */
L200:
	polar_(iopt, &m, x, y, z__, w, (E_fp)rad1_, &s, &nuest, &nvest, &eps, &nu, tu, &nv, tv, u, v, c__, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
	nc = (nu - 4) * (nv - 4);
/*  we calculate the function values at the different points. */
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    f[i__ - 1] = evapol_(tu, &nu, tv, &nv, c__, (E_fp)rad1_, &x[i__ - 1], &y[i__ - 1]);
/* L220: */
	}
	s_wsfe(&io___45);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
	goto L300;
/*  determination of the spline approximation on the ellips. */
L250:
	polar_(iopt, &m, x, y, z__, w, (E_fp)rad2_, &s, &nuest, &nvest, &eps, &nu, tu, &nv, tv, u, v, c__, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
/*  we determine the b-spline coefficients for the spline approximations */
/*  of the given function. */
	nc = (nu - 4) * (nv - 4);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__[i__ - 1] += .4f;
/* L260: */
	}
/*  we calculate the function values at the different points. */
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    f[i__ - 1] = evapol_(tu, &nu, tv, &nv, c__, (E_fp)rad2_, &x[i__ - 1], &y[i__ - 1]);
/* L270: */
	}
	if (iopt[0] < 0) {
	    goto L280;
	}
	s_wsfe(&io___46);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
	goto L300;
L280:
	s_wsfe(&io___47);
	e_wsfe();
L300:
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___49);
	do_fio(&c__1, (char *)&nu, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___50);
	e_wsfe();
	s_wsfe(&io___51);
	i__1 = nu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tu[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___52);
	do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___53);
	e_wsfe();
	s_wsfe(&io___54);
	i__1 = nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tv[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___55);
	e_wsfe();
	s_wsfe(&io___56);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  we determine mean and maximum errors. */
	sum = 0.f;
	ermax = 0.f;
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    error = (r__1 = f[i__ - 1] - exact[i__ - 1], dabs(r__1));
	    sum += error;
	    if (error > ermax) {
		ermax = error;
	    }
/* L350: */
	}
	sum /= ai;
	s_wsfe(&io___57);
	e_wsfe();
	s_wsfe(&io___58);
	e_wsfe();
	s_wsfe(&io___59);
	i__1 = m;
	for (l = 2; l <= i__1; l += 3) {
	    do_fio(&c__1, (char *)&x[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&y[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&f[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___60);
	do_fio(&c__1, (char *)&sum, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ermax, (ftnlen)sizeof(real));
	e_wsfe();
/* L400: */
    }
    s_stop("", 0L);
/*  format statements */
    return 0;
} /* MAIN__ */

doublereal rad1_(real *v)
{
    /* System generated locals */
    real ret_val;

/*  function program rad1 defines in polar coordinates, the boundary of */
/*  the approximation domain  x**2+y**2<=1. */
    ret_val = 1.f;
    return ret_val;
} /* rad1_ */

doublereal rad2_(real *v)
{
    /* System generated locals */
    real ret_val;

    /* Builtin functions */
    double sin(doublereal), sqrt(doublereal);

/*  function program rad2 defines in polar coordinates, the boundary of */
/*  the approximation domain  3*x**2+3*y**2-4*x*y<=1. */
    ret_val = 1.f / sqrt(3.f - sin(*v + *v) * 2.f);
    return ret_val;
} /* rad2_ */

doublereal testpo_(real *x, real *y)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

/*  function program testpo evaluates the test function for the polar */
/*  package. */
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = *y;
/* Computing 2nd power */
    r__3 = *x + *y;
    ret_val = (r__1 * r__1 + r__2 * r__2) / (r__3 * r__3 + .5f);
    return ret_val;
} /* testpo_ */

