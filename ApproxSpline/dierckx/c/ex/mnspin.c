/* mnspin.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static real c_b4 = .8f;
static integer c__20 = 20;
static real c_b7 = .4f;
static real c_b10 = .3f;
static real c_b13 = .1f;
static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mnspin : splint test program                       cc */
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
    static char fmt_925[] = "(\0020\002,5x,\002a\002,6x,\002b\002,6x,\002splint\002,7x,\002exint\002)";
    static char fmt_930[] = "(1x,2f7.1,2f12.5)";

    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double pow_ri(real *, integer *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real aint, a, b, c__[20];
    static integer i__, j, k, n;
    static real t[20], exint;
    static integer k1;
    static real ak;
    extern /* Subroutine */ int insert_(integer *, real *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
    extern doublereal splint_(real *, integer *, real *, integer *, real *, real *, real *);
    static integer nk1, ier;
    static real wrk[20];

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___18 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_930, 0 };


/*  as an example we calculate some integrals of the form */
/*          / b */
/*         !     (1-x)**k  dx */
/*      a / */

/*  main loop for the different spline degrees. */
    for (k = 1; k <= 5; ++k) {
	k1 = k + 1;
	ak = (real) k1;
/*  find the b-spline representation of the polynomial (1-x)**k. */
	n = k1 << 1;
	j = n;
	i__1 = k1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__[i__ - 1] = 0.f;
	    t[i__ - 1] = 0.f;
	    t[j - 1] = 1.f;
	    --j;
/* L10: */
	}
	c__[0] = 1.f;
/*  insert a number of knots */
	insert_(&c__0, t, &n, c__, &k, &c_b4, t, &n, c__, &c__20, &ier);
	insert_(&c__0, t, &n, c__, &k, &c_b7, t, &n, c__, &c__20, &ier);
	insert_(&c__0, t, &n, c__, &k, &c_b10, t, &n, c__, &c__20, &ier);
	insert_(&c__0, t, &n, c__, &k, &c_b13, t, &n, c__, &c__20, &ier);
/*  print the data for the spline. */
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___11);
	e_wsfe();
	s_wsfe(&io___12);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nk1 = n - k1;
	s_wsfe(&io___14);
	e_wsfe();
	s_wsfe(&io___15);
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  loop for the different integration limits a and b. */
	a = 0.f;
	b = 1.f;
	s_wsfe(&io___18);
	e_wsfe();
	for (j = 1; j <= 4; ++j) {
/*  calculate the value of the spline integral */
	    aint = splint_(t, &n, c__, &k, &a, &b, wrk);
/*  calculate the exact value of the integral */
	    r__1 = 1.f - a;
	    r__2 = 1.f - b;
	    exint = (pow_ri(&r__1, &k1) - pow_ri(&r__2, &k1)) / ak;
	    s_wsfe(&io___22);
	    do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&aint, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&exint, (ftnlen)sizeof(real));
	    e_wsfe();
	    a += .1f;
	    b += -.3f;
/* L20: */
	}
/* L30: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

