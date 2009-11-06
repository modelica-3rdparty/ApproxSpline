/* mndbin.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c              mndbin : dblint test program                          cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020tensor product spline of degrees\002,2i3)";
    static char fmt_910[] = "(1x,\002position of the knots in the x-direction\002)";
    static char fmt_920[] = "(1x,15f5.1)";
    static char fmt_930[] = "(1x,\002position of the knots in the y-direction\002)";
    static char fmt_940[] = "(\002 b-spline coefficients \002)";
    static char fmt_950[] = "(1x,8f9.4)";
    static char fmt_960[] = "(\0020\002,5x,\002xb\002,5x,\002xe\002,5x,\002yb\002,5x,\002ye\002,6x,\002dblint\002,7x,\002exint\002)";
    static char fmt_970[] = "(1x,4f7.1,2f12.5)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real facx, aint, c__[100];
    static integer i__, j;
    static real x[6], y[6], exint;
    static integer m0, m1, m2, m3, nc;
    static real xb, yb, xe, ye;
    static integer kx, ky, mx, my, nx, ny;
    extern doublereal dblint_(real *, integer *, real *, integer *, real *, integer *, integer *, real *, real *, real *, real *, real *);
    static real tx[15], ty[15];
    static integer kx1, ky1;
    static real fac, wrk[50];
    static integer nkx1, nky1;

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_970, 0 };


/*  we set up the end points of the integration domains. */
    mx = 6;
    my = 6;
    for (i__ = 1; i__ <= 6; ++i__) {
	x[i__ - 1] = (i__ - 1) * .2f;
	y[i__ - 1] = x[i__ - 1];
/* L10: */
    }
/*  loop for different spline degrees with respect to the x-variable */
    for (kx = 1; kx <= 5; kx += 2) {
/*  the knots in the x-direction */
	tx[kx + 1] = .4f;
	tx[kx + 2] = .7f;
	tx[kx + 3] = .9f;
	kx1 = kx + 1;
	nx = (kx1 << 1) + 3;
	j = nx;
	i__1 = kx1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tx[i__ - 1] = 0.f;
	    tx[j - 1] = 1.f;
	    --j;
/* L20: */
	}
/*  loop for different spline degrees with respect to the y-variable */
	for (ky = 2; ky <= 3; ++ky) {
/*  the knots in the y-direction */
	    ty[ky + 1] = .3f;
	    ty[ky + 2] = .8f;
	    ky1 = ky + 1;
	    ny = (ky1 << 1) + 2;
	    j = ny;
	    i__1 = ky1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ty[i__ - 1] = 0.f;
		ty[j - 1] = 1.f;
		--j;
/* L30: */
	    }
/*  we generate the b-spline coefficients for the test function x*y */
	    nkx1 = nx - kx1;
	    nky1 = ny - ky1;
	    i__1 = nky1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		c__[i__ - 1] = 0.f;
/* L40: */
	    }
	    i__1 = nkx1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		c__[(i__ - 1) * nky1] = 0.f;
/* L50: */
	    }
	    fac = (real) (kx * ky);
	    m0 = 1;
	    i__1 = nkx1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		m1 = m0 + nky1;
		facx = (tx[i__ + kx - 1] - tx[i__ - 1]) / fac;
		i__2 = nky1;
		for (j = 2; j <= i__2; ++j) {
		    m2 = m0 + 1;
		    m3 = m1 + 1;
		    c__[m3 - 1] = c__[m1 - 1] + c__[m2 - 1] - c__[m0 - 1] + facx * (ty[j + ky - 1] - ty[j - 1]);
		    ++m0;
		    ++m1;
/* L60: */
		}
		++m0;
/* L70: */
	    }
/*  printing of the spline information */
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&kx, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ky, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___25);
	    e_wsfe();
	    s_wsfe(&io___26);
	    i__1 = nx;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&tx[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___27);
	    e_wsfe();
	    s_wsfe(&io___28);
	    i__1 = ny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&ty[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    nc = nkx1 * nky1;
	    s_wsfe(&io___30);
	    e_wsfe();
	    s_wsfe(&io___31);
	    i__1 = nc;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/*  integration of the tensor product spline */
	    s_wsfe(&io___32);
	    e_wsfe();
	    j = 6;
	    for (i__ = 1; i__ <= 3; ++i__) {
		xb = x[i__ - 1];
		yb = xb;
		xe = x[j - 1];
		ye = xe;
		--j;
		aint = dblint_(tx, &nx, ty, &ny, c__, &kx, &ky, &xb, &xe, &yb, &ye, wrk);
		exint = (xe - xb) * (xe + xb) * (ye - yb) * (ye + yb) * .25f;
		s_wsfe(&io___40);
		do_fio(&c__1, (char *)&xb, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&xe, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&yb, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ye, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&aint, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&exint, (ftnlen)sizeof(real));
		e_wsfe();
/* L100: */
	    }
/* L200: */
	}
/* L300: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

