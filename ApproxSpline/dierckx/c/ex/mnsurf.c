/* mnsurf.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c        mnsurf : surfit test program                                cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(i3)";
    static char fmt_905[] = "(\0021\002,i3,\002 data points\002)";
    static char fmt_910[] = "(\0020\002,2(2x,\002i\002,5x,\002x(i)\002,6x,\002y(i)\002,6x,\002z(i)\002,6x))";
    static char fmt_915[] = "(3f10.4)";
    static char fmt_920[] = "(1x,2(i3,3f10.4,5x))";
    static char fmt_925[] = "(e20.6)";
    static char fmt_930[] = "(1x,\002estimate of standard deviation of z(i) =\002,e15.6)";
    static char fmt_935[] = "(\0020least-squares spline of degrees\002,2i3)";
    static char fmt_940[] = "(\0020smoothing spline of degrees\002,2i3)";
    static char fmt_945[] = "(\002 smoothing factor s=\002,f9.0)";
    static char fmt_950[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i3)";
    static char fmt_955[] = "(1x,\002total number of knots in the x-direction =\002,i3)";
    static char fmt_960[] = "(1x,\002position of the knots \002)";
    static char fmt_965[] = "(5x,10f7.3)";
    static char fmt_970[] = "(1x,\002total number of knots in the y-direction =\002,i3)";
    static char fmt_975[] = "(\0020b-spline coefficients \002)";
    static char fmt_980[] = "(5x,8f9.4)";
    static char fmt_1000[] = "(\0020\002,\002spline evaluation on a given grid\002)";
    static char fmt_985[] = "(\0020\002,\002x\002,2x,11f7.1)";
    static char fmt_990[] = "(3x,\002y\002)";
    static char fmt_995[] = "(1x,f4.1,11f7.3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void), s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer nmax, iopt, iwrk[300], kwrk, lwrk1, lwrk2;
    static real c__[200];
    static integer i__, j, m;
    static real s, w[80], x[80], y[80], z__[80], delta;
    static integer nxest, nyest;
    static real ai;
    static integer nc;
    static real fp, xb, yb;
    static integer is;
    static real xe, ye;
    static integer kx, ky, mx, my, nx, ny;
    static real tx[15], ty[15], ww, xx[11], yy[11], zz[121];
    extern /* Subroutine */ int bispev_(real *, integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, integer *, integer *, integer *), surfit_(integer *, integer *, real *, real *, real *, real *, real *, real *, real *, real *, integer *, integer *, real *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, integer *, integer *, integer *);
    static integer ier;
    static real eps, wrk1[12000], wrk2[6000];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 5, 0, fmt_900, 0 };
    static cilist io___3 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___4 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___6 = { 0, 5, 0, fmt_915, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___12 = { 0, 5, 0, fmt_925, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_990, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_995, 0 };


/*  we fetch the number of data points */
    s_rsfe(&io___1);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_rsfe();
    s_wsfe(&io___3);
    do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
    e_wsfe();
/*  we fetch the co-ordinate and function values of each data point. */
    s_wsfe(&io___4);
    e_wsfe();
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_rsfe(&io___6);
	do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z__[i__ - 1], (ftnlen)sizeof(real));
	e_rsfe();
	if (i__ / 2 << 1 != i__) {
	    goto L10;
	}
	j = i__ - 1;
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z__[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&z__[i__ - 1], (ftnlen)sizeof(real));
	e_wsfe();
L10:
	;
    }
/*  we fetch an estimate of the standard deviation of the data values. */
    s_rsfe(&io___12);
    do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(real));
    e_rsfe();
    s_wsfe(&io___14);
    do_fio(&c__1, (char *)&delta, (ftnlen)sizeof(real));
    e_wsfe();
/*  the weights are set equal to delta**(-1) */
    ww = 1.f / delta;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = ww;
/* L20: */
    }
/*  we set up the boundaries of the approximation domain. */
    xb = -2.f;
    xe = 2.f;
    yb = -2.f;
    ye = 2.f;
/* we generate a rectangular grid for evaluating the splines. */
    mx = 11;
    my = 11;
    for (i__ = 1; i__ <= 11; ++i__) {
	ai = (real) (i__ - 6);
	xx[i__ - 1] = ai * .4f;
	yy[i__ - 1] = xx[i__ - 1];
/* L30: */
    }
/*  we set up the dimension information */
    nxest = 15;
    nyest = 15;
    nmax = 15;
    kwrk = 300;
    lwrk1 = 12000;
    lwrk2 = 6000;
/*  we choose a value for eps */
    eps = 1e-6f;
/*  main loop for the different spline approximations. */
    for (is = 1; is <= 6; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L150;
	    case 6:  goto L160;
	}
/*  we start computing the least-squares bicubic polynomial (large s) */
L110:
	iopt = 0;
	kx = 3;
	ky = 3;
	s = 9e5f;
	goto L200;
/*  iopt=1 from the second call on. */
L120:
	iopt = 1;
	s = 200.f;
	goto L200;
/*  a value for s within its confidence interval */
L130:
	s = (real) m;
	goto L200;
/*  overfitting (s too small) */
L140:
	s = 20.f;
	goto L200;
/*  we change the degrees of the spline */
L150:
	iopt = 0;
	kx = 5;
	ky = 5;
	s = (real) m;
	goto L200;
/*  finally, we also calculate a least-squares spline approximation */
/*  with specified knots. */
L160:
	iopt = -1;
	kx = 3;
	ky = 3;
	nx = 11;
	ny = 11;
	j = kx + 2;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ai = (real) (i__ - 2);
	    tx[j - 1] = ai;
	    ty[j - 1] = ai;
	    ++j;
/* L170: */
	}
/*  determination of the spline approximation. */
L200:
	surfit_(&iopt, &m, x, y, z__, w, &xb, &xe, &yb, &ye, &kx, &ky, &s, &nxest, &nyest, &nmax, &eps, &nx, tx, &ny, ty, c__, &fp, wrk1, &lwrk1, wrk2, &lwrk2, iwrk, &kwrk, &ier);
/*  printing of the fitting results. */
	if (iopt >= 0) {
	    goto L210;
	}
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___49);
	do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___50);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___51);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___52);
	do_fio(&c__1, (char *)&nx, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___53);
	e_wsfe();
	s_wsfe(&io___54);
	i__1 = nx;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tx[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___55);
	do_fio(&c__1, (char *)&ny, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___56);
	e_wsfe();
	s_wsfe(&io___57);
	i__1 = ny;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&ty[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nc = (nx - kx - 1) * (ny - ky - 1);
	s_wsfe(&io___59);
	e_wsfe();
	s_wsfe(&io___60);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  evaluation of the spline approximation. */
	bispev_(tx, &nx, ty, &ny, c__, &kx, &ky, xx, &mx, yy, &my, zz, wrk2, &lwrk2, iwrk, &kwrk, &ier);
	s_wsfe(&io___62);
	e_wsfe();
	s_wsfe(&io___63);
	i__1 = mx;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&xx[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___64);
	e_wsfe();
	i__1 = my;
	for (j = 1; j <= i__1; ++j) {
	    s_wsfe(&io___65);
	    do_fio(&c__1, (char *)&yy[j - 1], (ftnlen)sizeof(real));
	    for (i__ = j; i__ <= 121; i__ += 11) {
		do_fio(&c__1, (char *)&zz[i__ - 1], (ftnlen)sizeof(real));
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

