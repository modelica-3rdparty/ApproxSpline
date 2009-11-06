/* mnevpo.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c              mnevpo : evapol test program                          cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020polar spline evaluation\002)";
    static char fmt_910[] = "(1x,\002position of the knots in the u-direction\002)";
    static char fmt_920[] = "(1x,15f5.1)";
    static char fmt_930[] = "(1x,\002position of the knots in the v-direction\002)";
    static char fmt_940[] = "(\002 b-spline coefficients \002)";
    static char fmt_950[] = "(5x,8f9.4)";
    static char fmt_960[] = "(\0020\002,\002f(x,y) corresponding to r(v)=1\002)";
    static char fmt_965[] = "(\0020\002,\002f(x,y) corresponding to r(v)=(1+cos(v)**2)/2\002)";
    static char fmt_970[] = "(\0020\002,\002x\002,2x,6f7.1)";
    static char fmt_975[] = "(3x,\002y\002)";
    static char fmt_980[] = "(1x,f4.1,6f7.3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real c__[42], f[36];
    static integer i__, j, m;
    static real x[6], y[6];
    static integer m0, m1, m2;
    extern doublereal r1_(), r2_();
    static integer nc, ir, nu, nv, mx, my;
    static real tu[11], tv[10];
    extern doublereal evapol_(real *, integer *, real *, integer *, real *, E_fp, real *, real *);
    static integer nu4, nv4;
    static real fac;

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___21 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_980, 0 };


/*  we set up the grid points for evaluating the polar spline. */
    mx = 6;
    my = 6;
    for (i__ = 1; i__ <= 6; ++i__) {
	x[i__ - 1] = ((i__ << 1) - 7) * .1f;
	y[i__ - 1] = x[i__ - 1];
/* L10: */
    }
/*  the interior knots with respect to the u-variable. */
    tu[4] = .4f;
    tu[5] = .7f;
    tu[6] = .9f;
    nu = 11;
/*  the interior knots with respect to the v-variable. */
    tv[4] = .3f;
    tv[5] = .8f;
    nv = 10;
/*  the boundary knots */
    for (i__ = 1; i__ <= 4; ++i__) {
	tu[i__ - 1] = 0.f;
	tv[i__ - 1] = 0.f;
	tu[i__ + 6] = 1.f;
	tv[i__ + 5] = 1.f;
/* L20: */
    }
/*  the number of b-spline coefficients */
    nu4 = nu - 4;
    nv4 = nv - 4;
    nc = nu4 * nv4;
/*  we generate the b-spline coefficients for the function s(u,v)=u**2 */
    m0 = 1;
    m1 = m0 + nv4;
    c__[m0 - 1] = 0.f;
    c__[m1 - 1] = 0.f;
    i__1 = nu4;
    for (i__ = 3; i__ <= i__1; ++i__) {
	m2 = m1 + nv4;
	c__[m2 - 1] = c__[m1 - 1] + (tu[i__ + 2] - tu[i__ - 1]) * ((c__[m1 - 1] - c__[m0 - 1]) / (tu[i__ + 1] - tu[i__ - 2]) + (tu[i__ + 1] - tu[i__ - 1]) / 3.f);
	m0 = m1;
	m1 = m2;
/* L70: */
    }
    i__1 = nu4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m0 = (i__ - 1) * nv4 + 1;
	fac = c__[m0 - 1];
	i__2 = nv4;
	for (j = 2; j <= i__2; ++j) {
	    ++m0;
	    c__[m0 - 1] = fac;
/* L80: */
	}
    }
    s_wsfe(&io___19);
    e_wsfe();
    s_wsfe(&io___20);
    e_wsfe();
    s_wsfe(&io___21);
    i__2 = nu;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, (char *)&tu[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___22);
    e_wsfe();
    s_wsfe(&io___23);
    i__2 = nv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_fio(&c__1, (char *)&tv[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___24);
    e_wsfe();
    s_wsfe(&io___25);
    i__2 = nc;
    for (j = 1; j <= i__2; ++j) {
	do_fio(&c__1, (char *)&c__[j - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
/* the spline s(u,v) defines a function f(x,y) through the transformation */
/*    x = r(v)*u*cos(v)   y = r(v)*u*sin(v) */
/* we consider two different functions r(v) */
    for (ir = 1; ir <= 2; ++ir) {
	switch (ir) {
	    case 1:  goto L110;
	    case 2:  goto L130;
	}
/* if r(v) =1 and s(u,v) = u**2 then f(x,y) = x**2+y**2 */
/* evaluation of f(x,y) */
L110:
	m = 0;
	i__2 = my;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = mx;
	    for (j = 1; j <= i__1; ++j) {
		++m;
		f[m - 1] = evapol_(tu, &nu, tv, &nv, c__, (E_fp)r1_, &x[j - 1], &y[i__ - 1]);
/* L120: */
	    }
	}
	s_wsfe(&io___29);
	e_wsfe();
	goto L200;
/* if r(v) = (1+cos(v)**2)/2 and s(u,v) = u**2 then f(x,y) = */
/*    4*(x**2+y**2)**3/(4*x**4 +y**4 +4*x**2*y**2) */
/* evaluation of f(x,y) */
L130:
	m = 0;
	i__1 = my;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = mx;
	    for (j = 1; j <= i__2; ++j) {
		++m;
		f[m - 1] = evapol_(tu, &nu, tv, &nv, c__, (E_fp)r2_, &x[j - 1], &y[i__ - 1]);
/* L140: */
	    }
	}
	s_wsfe(&io___30);
	e_wsfe();
L200:
	s_wsfe(&io___31);
	i__2 = mx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___32);
	e_wsfe();
	m1 = 0;
	i__2 = my;
	for (j = 1; j <= i__2; ++j) {
	    m0 = m1 + 1;
	    m1 += mx;
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&y[j - 1], (ftnlen)sizeof(real));
	    i__1 = m1;
	    for (m = m0; m <= i__1; ++m) {
		do_fio(&c__1, (char *)&f[m - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/* L210: */
	}
/* L300: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

doublereal r1_(real *v)
{
    /* System generated locals */
    real ret_val;

    ret_val = 1.f;
    return ret_val;
} /* r1_ */

doublereal r2_(real *v)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double cos(doublereal);

/* Computing 2nd power */
    r__1 = cos(*v);
    ret_val = (r__1 * r__1 + 1.f) * .5f;
    return ret_val;
} /* r2_ */

