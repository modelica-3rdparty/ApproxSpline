/* mncuev.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mncuev : curev test program                        cc */
/* c             evaluation of a closed planar curve                    cc */
/* c                    x = sx(u) , y = sy(u)                           cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020degree of the spline curve k =\002,i2)";
    static char fmt_905[] = "(1x,\002position of the knots\002)";
    static char fmt_910[] = "(5x,12f6.1)";
    static char fmt_915[] = "(1x,\002b-spline coefficients of sx(u)\002)";
    static char fmt_920[] = "(5x,14f5.0)";
    static char fmt_925[] = "(1x,\002b-spline coefficients of sy(u)\002)";
    static char fmt_930[] = "(1x,2(5x,\002i\002,3x,\002u(i)\002,4x,\002sx(u(i))\002,4x,\002sy(u(i))\002))";
    static char fmt_935[] = "(1x,2(i6,f7.2,2f12.5))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer idim;
    static real c__[40];
    static integer i__, j, k, m, n;
    static real t[20], u[20];
    extern /* Subroutine */ int curev_(integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, integer *);
    static integer i1, i2, j1, j2, j3, j4, k1;
    static real ai;
    static integer nc, jn, nk;
    static real sp[40];
    static integer mx, nk1, ier;
    static real per;

    /* Fortran I/O blocks */
    static cilist io___23 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_935, 0 };


/*  we have a planar curve */
    idim = 2;
/*  set up the dimension information */
    nc = 40;
    mx = 40;
/*  set up the points where the curve will be evaluated. */
    m = 20;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 1);
	u[i__ - 1] = ai * .05f;
/* L10: */
    }
/*  main loop for the different spline degrees. */
    for (k = 1; k <= 5; ++k) {
	k1 = k + 1;
/*  n denotes the total number of knots. */
	n = (k1 << 1) + 4;
/*  set up the knots of the spline */
	t[k1 - 1] = 0.f;
	t[k1] = .1f;
	t[k1 + 1] = .3f;
	t[k1 + 2] = .4f;
	t[k1 + 3] = .8f;
	t[k1 + 4] = 1.f;
/*  fetch the b-spline coefficients for sx(u) */
	c__[0] = 1.f;
	c__[1] = 3.f;
	c__[2] = 4.f;
	c__[3] = 5.f;
	c__[4] = -1.f;
/*  fetch the b-spline coefficients for sy(u) */
	c__[n] = 1.f;
	c__[n + 1] = 2.f;
	c__[n + 2] = -3.f;
	c__[n + 3] = 2.f;
	c__[n + 4] = 4.f;
/*  incorporate the boundary conditions for periodic splines */
	nk = n - k;
	per = t[nk - 1] - t[k1 - 1];
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
/*  the boundary knots */
	    i1 = nk + j;
	    i2 = nk - j;
	    j1 = k1 + j;
	    j2 = k1 - j;
	    t[i1 - 1] = t[j1 - 1] + per;
	    t[j2 - 1] = t[i2 - 1] - per;
/*  the boundary coefficients */
	    jn = j + n;
	    c__[j + 4] = c__[j - 1];
	    c__[jn + 4] = c__[jn - 1];
/* L20: */
	}
/*  evaluate the spline curve */
	curev_(&idim, t, &n, c__, &nc, &k, u, &m, sp, &mx, &ier);
/*  print the results. */
	s_wsfe(&io___23);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___24);
	e_wsfe();
	s_wsfe(&io___25);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___26);
	e_wsfe();
	nk1 = n - k1;
	s_wsfe(&io___28);
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___29);
	e_wsfe();
	i1 = n + 1;
	i2 = n + nk1;
	s_wsfe(&io___30);
	i__1 = i2;
	for (i__ = i1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___31);
	e_wsfe();
	i2 = 0;
	j4 = 0;
	for (j = 1; j <= 10; ++j) {
	    i1 = i2 + 1;
	    i2 = i1 + 1;
	    j1 = j4 + 1;
	    j2 = j1 + 1;
	    j3 = j2 + 1;
	    j4 = j3 + 1;
	    s_wsfe(&io___34);
	    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u[i1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j2 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u[i2 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j3 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j4 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L40: */
	}
/* L50: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

