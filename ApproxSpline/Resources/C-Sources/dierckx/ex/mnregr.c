/* mnregr.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__132 = 132;
static integer c__22 = 22;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c               mnregr : regrid test program                         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(i2)";
    static char fmt_905[] = "(11f5.1)";
    static char fmt_910[] = "(11f7.4)";
    static char fmt_915[] = "(\0021the input data\002)";
    static char fmt_920[] = "(\0020\002,8x,\002y\002,4x,6(4x,f4.1))";
    static char fmt_925[] = "(\002 \002,7x,\002x\002)";
    static char fmt_930[] = "(6x,f4.1,5x,6f8.4)";
    static char fmt_935[] = "(\0020least-squares spline of degrees\002,2i3)";
    static char fmt_940[] = "(\0020smoothing spline of degrees\002,2i3)";
    static char fmt_945[] = "(\002 smoothing factor s=\002,f7.2)";
    static char fmt_950[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i3)";
    static char fmt_955[] = "(1x,\002total number of knots in the x-direction =\002,i3)";
    static char fmt_960[] = "(1x,\002position of the knots \002)";
    static char fmt_965[] = "(5x,10f6.2)";
    static char fmt_970[] = "(1x,\002total number of knots in the y-direction =\002,i3)";
    static char fmt_975[] = "(\0020b-spline coefficients \002)";
    static char fmt_980[] = "(5x,8f9.4)";
    static char fmt_985[] = "(\0020\002,\002spline values at selected grid points\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void), s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer iopt, iwrk[60], kwrk, lwrk;
    static real c__[300], f[121];
    static integer i__, j, m;
    static real s, x[11], y[11], z__[121];
    static integer m1, m2, nxest, nyest;
    static real ai;
    static integer nc;
    static real fp, xb, yb;
    static integer is;
    static real xe, ye;
    static integer iw[22];
    static real wk[132];
    static integer kx, ky, mx, my, nx, ny;
    static real tx[17], ty[17];
    extern /* Subroutine */ int regrid_(integer *, integer *, real *, integer *, real *, real *, real *, real *, real *, real *, integer *, integer *, real *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, integer *, integer *, integer *, integer *), bispev_(real *, integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, integer *, integer *, integer *);
    static integer ier;
    static real wrk[850];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 5, 0, fmt_900, 0 };
    static cilist io___3 = { 0, 5, 0, fmt_905, 0 };
    static cilist io___6 = { 0, 5, 0, fmt_900, 0 };
    static cilist io___8 = { 0, 5, 0, fmt_905, 0 };
    static cilist io___11 = { 0, 5, 0, fmt_910, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___21 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___54 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___63 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___65 = { 0, 6, 0, fmt_930, 0 };


/*  we fetch the number of x-coordinates of the grid. */
    s_rsfe(&io___1);
    do_fio(&c__1, (char *)&mx, (ftnlen)sizeof(integer));
    e_rsfe();
/*  we fetch the x-coordinates of the grid. */
    s_rsfe(&io___3);
    i__1 = mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
    }
    e_rsfe();
/*  we fetch the number of y-coordinates of the grid. */
    s_rsfe(&io___6);
    do_fio(&c__1, (char *)&my, (ftnlen)sizeof(integer));
    e_rsfe();
/*  we fetch the y-coordinates of the grid. */
    s_rsfe(&io___8);
    i__1 = my;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
    }
    e_rsfe();
/*  we fetch the function values at the grid points. */
    m = mx * my;
    s_rsfe(&io___11);
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&z__[i__ - 1], (ftnlen)sizeof(real));
    }
    e_rsfe();
/*  printing of the input data. */
    s_wsfe(&io___13);
    e_wsfe();
    s_wsfe(&io___14);
    for (i__ = 1; i__ <= 6; ++i__) {
	do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___15);
    e_wsfe();
    m1 = 1;
    i__1 = mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m2 = m1 + 5;
	s_wsfe(&io___18);
	do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	i__2 = m2;
	for (j = m1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&z__[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	m1 += my;
/* L10: */
    }
    s_wsfe(&io___20);
    i__1 = my;
    for (i__ = 7; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___21);
    e_wsfe();
    m1 = 7;
    i__1 = mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m2 = m1 + 4;
	s_wsfe(&io___22);
	do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	i__2 = m2;
	for (j = m1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&z__[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	m1 += my;
/* L20: */
    }
/*  we set up the boundaries of the approximation domain. */
    xb = x[0];
    yb = y[0];
    xe = x[mx - 1];
    ye = y[my - 1];
/*  we set up the dimension information */
    nxest = 17;
    nyest = 17;
    lwrk = 850;
    kwrk = 60;
/*  main loop for the different spline approximations */
    for (is = 1; is <= 6; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L150;
	    case 6:  goto L160;
	}
/*  we start computing the least-squares bicubic polynomial */
L110:
	iopt = 0;
	kx = 3;
	ky = 3;
	s = 10.f;
	goto L200;
/*  iopt=1 from the second call on */
L120:
	iopt = 1;
	s = .22f;
	goto L200;
/*  overfitting (s too small) */
L130:
	s = .1f;
	goto L200;
/*  an interpolating spline */
L140:
	s = 0.f;
	goto L200;
/*  we change the degrees of the spline */
L150:
	kx = 5;
	ky = 5;
	s = .2f;
	iopt = 0;
	goto L200;
/*  finally we also calculate a least-squares spline approximation */
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
	    tx[j - 1] = ai * .5f;
	    ty[j - 1] = tx[j - 1];
	    ++j;
/* L170: */
	}
/*  determination of the spline approximation. */
L200:
	regrid_(&iopt, &mx, x, &my, y, z__, &xb, &xe, &yb, &ye, &kx, &ky, &s, &nxest, &nyest, &nx, tx, &ny, ty, c__, &fp, wrk, &lwrk, iwrk, &kwrk, &ier);
/*  printing of the fitting results. */
	if (iopt >= 0) {
	    goto L210;
	}
	s_wsfe(&io___46);
	do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___47);
	do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___48);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___49);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___50);
	do_fio(&c__1, (char *)&nx, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___51);
	e_wsfe();
	s_wsfe(&io___52);
	i__1 = nx;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tx[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___53);
	do_fio(&c__1, (char *)&ny, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___54);
	e_wsfe();
	s_wsfe(&io___55);
	i__1 = ny;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&ty[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nc = (nx - kx - 1) * (ny - ky - 1);
	s_wsfe(&io___57);
	e_wsfe();
	s_wsfe(&io___58);
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  evaluation of the spline approximation. */
	bispev_(tx, &nx, ty, &ny, c__, &kx, &ky, x, &mx, y, &my, f, wk, &c__132, iw, &c__22, &ier);
	s_wsfe(&io___62);
	e_wsfe();
	s_wsfe(&io___63);
	i__1 = my;
	for (i__ = 1; i__ <= i__1; i__ += 2) {
	    do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___64);
	e_wsfe();
	m1 = 1;
	i__1 = mx;
	for (i__ = 1; i__ <= i__1; i__ += 2) {
	    m2 = m1 + my - 1;
	    s_wsfe(&io___65);
	    do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	    i__2 = m2;
	    for (j = m1; j <= i__2; j += 2) {
		do_fio(&c__1, (char *)&f[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    m1 += my << 1;
/* L230: */
	}
/* L300: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

