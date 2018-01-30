/* mnspev.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mnspev : splev test program                        cc */
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
    static char fmt_925[] = "(1x,\002error flag ier =\002,i3)";
    static char fmt_930[] = "(1x,3(7x,\002i\002,3x,\002x(i)\002,4x,\002s(x(i))\002))";
    static char fmt_935[] = "(1x,3(i8,f7.2,f11.5))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real c__[20];
    static integer i__, j, k, m, n;
    static real t[20], x[21], y[21];
    extern /* Subroutine */ int splev_(real *, integer *, real *, integer *, real *, real *, integer *, integer *);
    static integer i1, i2, k1;
    static real ai;
    static integer nk1, ier;

    /* Fortran I/O blocks */
    static cilist io___14 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___17 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_935, 0 };


/*  set up the points where the splines will be evaluated. */
    m = 21;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 1);
	x[i__ - 1] = ai * .05f;
/* L10: */
    }
/*  main loop for the different spline degrees. */
    for (k = 1; k <= 5; ++k) {
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
/*  evaluate the spline. */
	splev_(t, &n, c__, &k, x, y, &m, &ier);
/*  print the results. */
	s_wsfe(&io___14);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___15);
	e_wsfe();
	s_wsfe(&io___16);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___17);
	e_wsfe();
	s_wsfe(&io___18);
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___19);
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___20);
	e_wsfe();
	i2 = 0;
	for (j = 1; j <= 7; ++j) {
	    i1 = i2 + 1;
	    i2 = i1 + 2;
	    s_wsfe(&io___23);
	    i__1 = i2;
	    for (i__ = i1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/* L40: */
	}
/* L50: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

