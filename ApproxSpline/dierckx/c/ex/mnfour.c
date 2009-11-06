/* mnfour.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mnfour : fourco test program                       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020degree of the spline k =\002,i2)";
    static char fmt_905[] = "(1x,\002position of the knots\002)";
    static char fmt_910[] = "(5x,12f5.1)";
    static char fmt_915[] = "(1x,\002b-spline coefficients\002)";
    static char fmt_920[] = "(5x,8f9.5)";
    static char fmt_925[] = "(\0020\002,2x,\002alfa\002,9x,\002ress\002,9x,\002exas\002,9x,\002resc\002,9x,\002exac\002)";
    static char fmt_930[] = "(1x,e8.1,4f13.5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real alfa[10], resc[10], ress[10], c__[20];
    static integer i__, j, k, m, n;
    static real t[20];
    static integer k1;
    static real ak, rc, rs;
    extern /* Subroutine */ int fourco_(real *, integer *, real *, real *, integer *, real *, real *, real *, real *, integer *), exfour_(real *, real *, real *);
    static integer nk1, ier;
    static real wrk1[20], wrk2[20];

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_930, 0 };


/*  as an example we calculate some integrals of the form */
/*          / 1                               / 1 */
/*         !    x * sin(alfa*x) dx  and      !   x * cos(alfa*x) dx */
/*      0 /                               0 / */

/*  we will represent y = x as a cubic spline. */
    k = 3;
    k1 = k + 1;
/*  we fetch the knots of the cubic spline */
    n = (k1 << 1) + 4;
/*  the boundary knots */
    j = n;
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t[i__ - 1] = 0.f;
	t[j - 1] = 1.f;
	--j;
/* L10: */
    }
/*  the interior knots */
    t[4] = .1f;
    t[5] = .3f;
    t[6] = .4f;
    t[7] = .8f;
/*  find the b-spline representation of y=x */
    nk1 = n - k1;
    ak = (real) k;
    c__[0] = 0.f;
    i__1 = nk1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = i__ + k;
	c__[i__ - 1] = c__[i__ - 2] + (t[j - 1] - t[i__ - 1]) / ak;
/* L20: */
    }
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
    s_wsfe(&io___13);
    e_wsfe();
    s_wsfe(&io___14);
    i__1 = nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
/*  fetch the different values for alfa */
    m = 8;
    alfa[0] = 0.f;
    alfa[1] = .001f;
    i__1 = m;
    for (i__ = 3; i__ <= i__1; ++i__) {
	alfa[i__ - 1] = -alfa[i__ - 2] * 10.f;
/* L30: */
    }
/*  calculate the fourier integrals of the cubic spline */
    fourco_(t, &n, c__, alfa, &m, ress, resc, wrk1, wrk2, &ier);
/*  print the results */
    s_wsfe(&io___22);
    e_wsfe();
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch the exact values of the integrals */
	exfour_(&alfa[i__ - 1], &rs, &rc);
	s_wsfe(&io___25);
	do_fio(&c__1, (char *)&alfa[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ress[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rs, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&resc[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&rc, (ftnlen)sizeof(real));
	e_wsfe();
/* L40: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

/* Subroutine */ int exfour_(real *alfa, real *rs, real *rc)
{
    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real half;
    static integer k;
    static real three, c1;
    static integer k2;
    static real s1, aa, cc, ak, ss, one;

/*  subroutine exfour calculates the integrals */
/*                 / 1 */
/*      rs =      !    x*sin(alfa*x) dx    and */
/*             0 / */
/*                 / 1 */
/*      rc =      !    x*cos(alfa*x) dx */
/*             0 / */
/*  ..function references.. */
/*  set constants */
    one = 1.f;
    three = 3.f;
    half = .5f;
    if (dabs(*alfa) < one) {
	goto L10;
    }
/*  integration by parts */
    aa = one / *alfa;
    cc = cos(*alfa);
    ss = sin(*alfa);
    *rs = (ss * aa - cc) * aa;
    *rc = ((cc - one) * aa + ss) * aa;
    goto L50;
/*  using the series expansions of sin(alfa*x) and cos(alfa*x) */
L10:
    *rs = 0.f;
    *rc = half;
    if (*alfa == 0.f) {
	goto L50;
    }
    *rs = *alfa / three;
    ss = *rs;
    cc = *rc;
    aa = -(*alfa) * *alfa;
    for (k = 1; k <= 21; ++k) {
	k2 = k - 1 << 1;
	ak = (real) ((k2 + 2) * (k2 + 5));
	ss = ss * aa / ak;
	s1 = *rs + ss;
	if (s1 == *rs) {
	    goto L30;
	}
	*rs = s1;
/* L20: */
    }
L30:
    for (k = 1; k <= 21; ++k) {
	k2 = k - 1 << 1;
	ak = (real) ((k2 + 1) * (k2 + 4));
	cc = cc * aa / ak;
	c1 = *rc + cc;
	if (c1 == *rc) {
	    goto L50;
	}
	*rc = c1;
/* L40: */
    }
L50:
    return 0;
} /* exfour_ */

