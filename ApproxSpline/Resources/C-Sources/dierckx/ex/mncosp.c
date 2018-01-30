/* mncosp.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mncosp : cocosp test program                       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real x[10] = { .25f,.5f,.75f,1.25f,1.75f,2.25f,2.75f,3.25f,6.25f,12.25f };
    static real y[10] = { 17.f,15.2f,13.8f,12.2f,11.f,10.1f,9.4f,8.6f,6.1f,3.5f };

    /* Format strings */
    static char fmt_900[] = "(\0020convex spline approximation\002)";
    static char fmt_905[] = "(\0020concave spline approximation\002)";
    static char fmt_910[] = "(\0020unconstrained spline approximation\002)";
    static char fmt_915[] = "(\002 error flag ier=\002,i2)";
    static char fmt_920[] = "(1x,\002sum of squared residuals sq=\002,e10.3)";
    static char fmt_925[] = "(1x,\002total number of knots n=\002,i2)";
    static char fmt_930[] = "(1x,\002position of the knots\002)";
    static char fmt_935[] = "(5x,8f7.2)";
    static char fmt_940[] = "(1x,\002the knots where s''(x)=0\002)";
    static char fmt_945[] = "(5x,f7.2)";
    static char fmt_950[] = "(1x,\002b-spline coefficients\002)";
    static char fmt_955[] = "(5x,4f12.4)";
    static char fmt_960[] = "(\0020 i\002,6x,\002x(i)\002,5x,\002y(i)\002,4x,\002s(x(i))\002,3x,\002s''(x(i))\002)";
    static char fmt_965[] = "(1x,i2,5x,f5.2,5x,f4.1,5x,f5.2,5x,f5.2)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static logical bind[20];
    static integer iwrk[450], kwrk, lwrk;
    static real c__[20], e[20];
    static integer i__, j, m, n;
    static real t[20], w[10];
    static integer maxtr, j3, n4, n6;
    static real s2[10];
    static integer is;
    static real sq, sx[10];
    static integer maxbin;
    extern /* Subroutine */ int cocosp_(integer *, real *, real *, real *, integer *, real *, real *, integer *, integer *, real *, real *, real *, logical *, real *, integer *, integer *, integer *, integer *), splder_(real *, integer *, real *, integer *, integer *, real *, real *, integer *, real *, integer *);
    static integer ier;
    static real wrk[550];

    /* Fortran I/O blocks */
    static cilist io___15 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___19 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_965, 0 };


/*  the absciss values of the data points. */
/*  the ordinate values of the data points. */
/*  m denotes the number of data points. */
    m = 10;
/*  we set up the weights of the data points. */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 1.f;
/* L10: */
    }
/*  we set up the dimension information. */
    maxtr = 100;
    maxbin = 10;
    lwrk = 550;
    kwrk = 450;
/*  we fetch the knots of the cubic spline */
    n = 11;
    n4 = n - 4;
    n6 = n - 6;
/*  the interior knots */
    t[4] = 1.6f;
    t[5] = 2.5f;
    t[6] = 6.f;
/*  the boundary knots */
    for (i__ = 1; i__ <= 4; ++i__) {
	t[i__ - 1] = x[0];
	t[i__ + 6] = x[m - 1];
/* L20: */
    }
/*  loop for the different spline approximations */
    for (is = 1; is <= 3; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L130;
	    case 3:  goto L150;
	}
/*  a convex spline approximation */
L110:
	s_wsfe(&io___15);
	e_wsfe();
	i__1 = n6;
	for (j = 1; j <= i__1; ++j) {
	    e[j - 1] = -1.f;
/* L120: */
	}
	goto L200;
/*  a concave spline approximation (a straight line) */
L130:
	s_wsfe(&io___18);
	e_wsfe();
	i__1 = n6;
	for (j = 1; j <= i__1; ++j) {
	    e[j - 1] = 1.f;
/* L140: */
	}
	goto L200;
/*  no convexity/concavity constraints */
L150:
	s_wsfe(&io___19);
	e_wsfe();
	i__1 = n6;
	for (j = 1; j <= i__1; ++j) {
	    e[j - 1] = 0.f;
/* L160: */
	}
/*  we determine the spline approximation. */
L200:
	cocosp_(&m, x, y, w, &n, t, e, &maxtr, &maxbin, c__, &sq, sx, bind, wrk, &lwrk, iwrk, &kwrk, &ier);
/*  printing of the results. */
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___28);
	do_fio(&c__1, (char *)&sq, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___29);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___30);
	e_wsfe();
	s_wsfe(&io___31);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___32);
	e_wsfe();
	i__1 = n6;
	for (j = 1; j <= i__1; ++j) {
	    j3 = j + 3;
	    if (bind[j - 1]) {
		s_wsfe(&io___34);
		do_fio(&c__1, (char *)&t[j3 - 1], (ftnlen)sizeof(real));
		e_wsfe();
	    }
/* L300: */
	}
	s_wsfe(&io___35);
	e_wsfe();
	s_wsfe(&io___36);
	i__1 = n4;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  we evaluate the second order derivative of the spline. */
	splder_(t, &n, c__, &c__3, &c__2, x, s2, &m, wrk, &ier);
	s_wsfe(&io___38);
	e_wsfe();
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s_wsfe(&io___39);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sx[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&s2[i__ - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L400: */
	}
/* L500: */
    }
    s_stop("", 0L);
/*  format statements */
    return 0;
} /* MAIN__ */

