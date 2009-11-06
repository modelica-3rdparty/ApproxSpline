/* mnspal.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mnspal : spalde test program                       cc */
/* c    evaluation of a spline function through its polynomial          cc */
/* c            representation in each knot interval.                   cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020degree of the spline k =\002,i2)";
    static char fmt_905[] = "(1x,\002position of the knots\002)";
    static char fmt_910[] = "(5x,15f5.1)";
    static char fmt_915[] = "(1x,\002b-spline coefficients\002)";
    static char fmt_920[] = "(5x,8f9.5)";
    static char fmt_925[] = "(\0020knot interval (\002,f4.1,\002,\002,f4.1,\002)\002)";
    static char fmt_930[] = "(1x,\002derivative values at the midpoint of the interval\002)";
    static char fmt_935[] = "(2x,6e13.5)";
    static char fmt_940[] = "(1x,\002coefficients in the polynomial representation\002)";
    static char fmt_945[] = "(\0020\002,3(7x,\002i\002,3x,\002x(i)\002,4x,\002s(x(i))\002))";
    static char fmt_950[] = "(1x,3(i8,f7.2,f11.5))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real c__[20], d__[6];
    static integer i__, j, k, l, m, n;
    static real t[20], x[21], y[21];
    static integer i1, i2, k1, l1;
    static real ai, aj;
    static integer jj;
    static real tt;
    extern /* Subroutine */ int spalde_(real *, integer *, real *, integer *, real *, real *, integer *);
    static real xx;
    static integer nk1;
    static real fac, cof[6], arg;
    static integer ier;
    static real pol;

    /* Fortran I/O blocks */
    static cilist io___12 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_950, 0 };


/*  set up the points where the splines will be evaluated. */
    m = 21;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 1);
	x[i__ - 1] = ai * .05f;
/* L10: */
    }
/*  main loop for the different spline degrees. */
    for (k = 3; k <= 5; k += 2) {
	k1 = k + 1;
/*  n denotes the total number of knots. */
	n = (k1 << 1) + 4;
/*  set up the knots of the spline */
	j = n;
/*  the boundary knots */
	i__1 = k1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t[i__ - 1] = 0.f;
	    t[j - 1] = 1.f;
	    --j;
/* L20: */
	}
/*  the interior knots */
	t[k1] = .1f;
	t[k1 + 1] = .3f;
	t[k1 + 2] = .4f;
	t[k1 + 3] = .8f;
/*  generate the b-spline coefficients. */
	nk1 = n - k1;
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ai = (real) i__;
	    c__[i__ - 1] = ai * .01f * (ai - 5.f);
/* L30: */
	}
/*  print the data for the spline. */
	s_wsfe(&io___12);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___13);
	e_wsfe();
	s_wsfe(&io___14);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___15);
	e_wsfe();
	s_wsfe(&io___16);
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	l = k;
	l1 = k1;
/*  main loop for the different points of evaluation. */
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    arg = x[i__ - 1];
/*  search for knot interval t(l)<=x(i)<t(l+1). */
L40:
	    if (arg < t[l1 - 1] || l == nk1) {
		goto L60;
	    }
/*  a new knot interval. */
	    l = l1;
	    l1 = l + 1;
	    if (t[l - 1] == t[l1 - 1]) {
		goto L40;
	    }
	    s_wsfe(&io___20);
	    do_fio(&c__1, (char *)&t[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&t[l1 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/*  calculate the spline derivatives at the midpoint tt of the interval */
	    tt = (t[l - 1] + t[l1 - 1]) * .5f;
	    spalde_(t, &n, c__, &k1, &tt, d__, &ier);
	    s_wsfe(&io___24);
	    e_wsfe();
	    s_wsfe(&io___25);
	    i__2 = k1;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&d__[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/*  calculate the coefficients cof in the polynomial representation of */
/*  the spline in the current knot interval,i.e. */
/*    s(x) = cof(1)+cof(2)*(x-tt)+...+cof(k1)*(x-tt)**k */
	    fac = 1.f;
	    i__2 = k1;
	    for (j = 1; j <= i__2; ++j) {
		cof[j - 1] = d__[j - 1] / fac;
		aj = (real) j;
		fac *= aj;
/* L50: */
	    }
	    s_wsfe(&io___29);
	    e_wsfe();
	    s_wsfe(&io___30);
	    i__2 = k1;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&cof[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    goto L40;
/*  evaluate the polynomial */
L60:
	    xx = arg - tt;
	    pol = cof[k1 - 1];
	    jj = k1;
	    i__2 = k;
	    for (j = 1; j <= i__2; ++j) {
		--jj;
		pol = pol * xx + cof[jj - 1];
/* L70: */
	    }
	    y[i__ - 1] = pol;
/* L80: */
	}
	s_wsfe(&io___35);
	e_wsfe();
	i2 = 0;
	for (j = 1; j <= 7; ++j) {
	    i1 = i2 + 1;
	    i2 = i1 + 2;
	    s_wsfe(&io___38);
	    i__1 = i2;
	    for (i__ = i1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/* L90: */
	}
/* L100: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

