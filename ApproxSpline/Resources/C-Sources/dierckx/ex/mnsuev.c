/* mnsuev.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__108 = 108;
static integer c__48 = 48;
static integer c__12 = 12;
static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c              mnsuev : surev test program                           cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020bicubic spline surface\002)";
    static char fmt_910[] = "(1x,\002position of the knots in the u-direction\002)";
    static char fmt_920[] = "(1x,15f5.1)";
    static char fmt_930[] = "(1x,\002position of the knots in the v-direction\002)";
    static char fmt_940[] = "(\002 b-spline coefficients \002)";
    static char fmt_950[] = "(5x,8f9.4)";
    static char fmt_960[] = "(\0020\002,\002spline values at selected grid points\002)";
    static char fmt_970[] = "(\0020\002,2x,\002v\002,6(3x,f4.1))";
    static char fmt_980[] = "(\002 \002,1x,\002u\002)";
    static char fmt_990[] = "(\002 \002,f4.1)";
    static char fmt_995[] = "(5x,6f7.3)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer idim, iwrk[12];
    static real c__[126], f[108];
    static integer i__, j, l, m;
    static real u[6], v[6];
    static integer m0, m1, m2, m3;
    extern /* Subroutine */ int surev_(integer *, real *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, integer *);
    static integer nc, mu, mv, nu, nv;
    static real tu[11], tv[10];
    static integer nu4, nv4;
    static real fac;
    static integer ier;
    static real wrk[48];

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_990, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_995, 0 };


/*  we set up the grid points for evaluating the spline surface. */
    mu = 6;
    mv = 6;
    for (i__ = 1; i__ <= 6; ++i__) {
	u[i__ - 1] = (i__ - 1) * .2f;
	v[i__ - 1] = u[i__ - 1];
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
/*  we generate the b-spline coefficients for the test surface */
/*        x = u*v    y = v**2    z = u+v     0 <= u,v <= 1 */
/*  the dimension of the surface */
    idim = 3;
/*  the number of b-spline coefficients for each co-ordinate */
    nu4 = nu - 4;
    nv4 = nv - 4;
    nc = nu4 * nv4;
/*  the coefficients for x = u*v */
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__ - 1] = 0.f;
/* L30: */
    }
    i__1 = nu4;
    for (i__ = 2; i__ <= i__1; ++i__) {
	c__[(i__ - 1) * nv4] = 0.f;
/* L40: */
    }
    m0 = 1;
    i__1 = nu4;
    for (i__ = 2; i__ <= i__1; ++i__) {
	m1 = m0 + nv4;
	fac = (tu[i__ + 2] - tu[i__ - 1]) / 9.f;
	i__2 = nv4;
	for (j = 2; j <= i__2; ++j) {
	    m2 = m0 + 1;
	    m3 = m1 + 1;
	    c__[m3 - 1] = c__[m1 - 1] + c__[m2 - 1] - c__[m0 - 1] + fac * (tv[j + 2] - tv[j - 1]);
	    ++m0;
	    ++m1;
/* L50: */
	}
	++m0;
/* L60: */
    }
/*  the coefficients for y = v**2. */
    l = nc;
    m0 = l + 1;
    m1 = m0 + 1;
    c__[m0 - 1] = 0.f;
    c__[m1 - 1] = 0.f;
    i__1 = nv4;
    for (i__ = 3; i__ <= i__1; ++i__) {
	c__[m1] = c__[m1 - 1] + (tv[i__ + 2] - tv[i__ - 1]) * ((c__[m1 - 1] - c__[m0 - 1]) / (tv[i__ + 1] - tv[i__ - 2]) + (tv[i__ + 1] - tv[i__ - 1]) / 3.f);
	m0 = m1;
	m1 = m0 + 1;
/* L70: */
    }
    i__1 = nv4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m0 = l + i__;
	fac = c__[m0 - 1];
	i__2 = nu4;
	for (j = 1; j <= i__2; ++j) {
	    m0 += nv4;
	    c__[m0 - 1] = fac;
/* L80: */
	}
    }
/*  the coefficients for z = u+v */
    l += nc;
    m0 = l + 1;
    c__[m0 - 1] = 0.f;
    i__2 = nv4;
    for (i__ = 2; i__ <= i__2; ++i__) {
	m1 = m0 + 1;
	c__[m1 - 1] = c__[m0 - 1] + (tv[i__ + 2] - tv[i__ - 1]) / 3.f;
	m0 = m1;
/* L90: */
    }
    i__2 = nv4;
    for (i__ = 1; i__ <= i__2; ++i__) {
	m0 = l + i__;
	i__1 = nu4;
	for (j = 2; j <= i__1; ++j) {
	    m1 = m0 + nv4;
	    c__[m1 - 1] = c__[m0 - 1] + (tu[j + 2] - tu[j - 1]) / 3.f;
	    m0 = m1;
/* L100: */
	}
    }
/*  evaluation of the spline surface */
    surev_(&idim, tu, &nu, tv, &nv, c__, u, &mu, v, &mv, f, &c__108, wrk, &c__48, iwrk, &c__12, &ier);
/*  printing of the results */
    s_wsfe(&io___26);
    e_wsfe();
    s_wsfe(&io___27);
    e_wsfe();
    s_wsfe(&io___28);
    i__1 = nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&tu[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___29);
    e_wsfe();
    s_wsfe(&io___30);
    i__1 = nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&tv[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___31);
    e_wsfe();
    m1 = 0;
    i__1 = idim;
    for (l = 1; l <= i__1; ++l) {
	m0 = m1 + 1;
	m1 += nc;
	s_wsfe(&io___32);
	i__2 = m1;
	for (j = m0; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&c__[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L110: */
    }
    s_wsfe(&io___33);
    e_wsfe();
    s_wsfe(&io___34);
    i__1 = mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&v[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___35);
    e_wsfe();
    m = mu * mv;
    m0 = 0;
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&u[i__ - 1], (ftnlen)sizeof(real));
	e_wsfe();
	m1 = m0;
	i__2 = idim;
	for (l = 1; l <= i__2; ++l) {
	    m2 = m1 + 1;
	    m3 = m1 + mv;
	    s_wsfe(&io___38);
	    i__3 = m3;
	    for (j = m2; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&f[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    m1 += m;
/* L120: */
	}
	m0 += mv;
/* L130: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

