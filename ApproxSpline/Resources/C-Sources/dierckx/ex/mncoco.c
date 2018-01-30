/* mncoco.f -- translated by f2c (version 19961017).
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
/* c                 mncoco : concon test program                       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real x[16] = { .1f,.3f,.5f,.7f,.9f,1.25f,1.75f,2.25f,2.75f,3.5f,4.5f,5.5f,6.5f,7.5f,8.5f,9.5f };
    static real y[16] = { .124f,.234f,.256f,.277f,.278f,.291f,.308f,.311f,.315f,.322f,.317f,.326f,.323f,.321f,.322f,.328f };

    /* Format strings */
    static char fmt_900[] = "(\0020upper limit for the sum of squared residuals s=\002,e8.1,5x,\002error flag ier=\002,i2)";
    static char fmt_905[] = "(1x,\002sum of squared residuals sq=\002,e10.3)";
    static char fmt_910[] = "(1x,\002total number of knots n=\002,i2)";
    static char fmt_915[] = "(1x,\002position of the knots\002)";
    static char fmt_920[] = "(5x,8f7.2)";
    static char fmt_925[] = "(1x,\002the knots where s''(x)=0\002)";
    static char fmt_930[] = "(5x,f7.2)";
    static char fmt_935[] = "(1x,\002b-spline coefficients\002)";
    static char fmt_940[] = "(5x,4f12.6)";
    static char fmt_945[] = "(\0020 i\002,5x,\002x(i)\002,6x,\002y(i)\002,4x,\002s(x(i))\002,4x,\002s''(x(i))\002)";
    static char fmt_950[] = "(1x,i2,5x,f4.2,5x,f5.3,5x,f5.3,5x,f8.4)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static logical bind[20];
    static integer nest, iopt, iwrk[450], kwrk, lwrk;
    static real c__[20];
    static integer i__, j, m, n;
    static real s[3], t[20], v[16], w[16];
    static integer maxtr, j3, n4, n6;
    static real s2[16];
    static integer is;
    static real sq, sx[16];
    static integer maxbin;
    extern /* Subroutine */ int concon_(integer *, integer *, real *, real *, real *, real *, real *, integer *, integer *, integer *, integer *, real *, real *, real *, real *, logical *, real *, integer *, integer *, integer *, integer *), splder_(real *, integer *, real *, integer *, integer *, real *, real *, integer *, real *, integer *);
    static integer ier;
    static real wrk[550];

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_950, 0 };


/*  the absciss values of the data points. */
/*  the ordinate values of the data points. */
/*  m denotes the number of data points. */
    m = 16;
/*  we set up the weights of the data points. */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 1.f;
/* L10: */
    }
    w[0] = 10.f;
    w[1] = 3.f;
    w[15] = 10.f;
/*  we will determine concave approximations */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = 1.f;
/* L20: */
    }
/*  we set up the dimension information. */
    nest = 20;
    maxtr = 100;
    maxbin = 10;
    lwrk = 550;
    kwrk = 450;
/*  we set up the different s-values. */
    s[0] = .2f;
    s[1] = .04f;
    s[2] = 2e-4f;
/*  initialization. */
    iopt = 0;
/*  loop for the different spline approximations */
    for (is = 1; is <= 3; ++is) {
/*  we determine the concave spline approximation. */
	concon_(&iopt, &m, x, y, w, v, &s[is - 1], &nest, &maxtr, &maxbin, &n, t, c__, &sq, sx, bind, wrk, &lwrk, iwrk, &kwrk, &ier);
/*  printing of the results. */
	s_wsfe(&io___24);
	do_fio(&c__1, (char *)&s[is - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___25);
	do_fio(&c__1, (char *)&sq, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___26);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___27);
	e_wsfe();
	s_wsfe(&io___28);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___29);
	e_wsfe();
	n6 = n - 6;
	i__1 = n6;
	for (j = 1; j <= i__1; ++j) {
	    j3 = j + 3;
	    if (bind[j - 1]) {
		s_wsfe(&io___33);
		do_fio(&c__1, (char *)&t[j3 - 1], (ftnlen)sizeof(real));
		e_wsfe();
	    }
/* L100: */
	}
	s_wsfe(&io___34);
	e_wsfe();
	n4 = n - 4;
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
/* L200: */
	}
/*  iopt=1 from the second call on. */
	iopt = 1;
/* L300: */
    }
    s_stop("", 0L);
/*  format statements */
    return 0;
} /* MAIN__ */

