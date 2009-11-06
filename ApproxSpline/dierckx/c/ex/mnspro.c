/* mnspro.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mnspro : sproot test program                       cc */
/* c      application : to find the intersection of a planar            cc */
/* c      cubic spline curve   x = sx(u)   y = sy(u)   with             cc */
/* c      a straight line   alfa*x + beta*y = gamma                     cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020degree of the spline curve k =\002,i2)";
    static char fmt_905[] = "(1x,\002position of the knots\002)";
    static char fmt_910[] = "(5x,7f6.1)";
    static char fmt_915[] = "(1x,\002b-spline coefficients of sx(u)\002)";
    static char fmt_920[] = "(5x,14f5.0)";
    static char fmt_925[] = "(1x,\002b-spline coefficients of sy(u)\002)";
    static char fmt_930[] = "(\0020intersection with\002,f6.1,\002 *x +\002,f5.1,\002 *y =\002,f5.1)";
    static char fmt_935[] = "(1x,\002number of intersection points m =\002,i3)";
    static char fmt_940[] = "(6x,\002i\002,7x,\002u(i)\002,5x,\002sx(u(i))\002,4x,\002sy(u(i))\002)";
    static char fmt_945[] = "(1x,i6,3f12.5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static real alfa, beta;
    static integer idim, mest;
    static real zero[20], c__[26];
    static integer i__, j, k, m, n;
    static real gamma, t[13];
    extern /* Subroutine */ int curev_(integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, integer *);
    static integer i1, i2, k1, l1, l2;
    static real cc[13];
    static integer nc, is;
    static real sp[40];
    static integer nk1;
    extern /* Subroutine */ int sproot_(real *, integer *, real *, real *, integer *, integer *, integer *);
    static integer ier;
    static real per;

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___14 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___16 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___17 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_945, 0 };


/*  we have a planar curve */
    idim = 2;
/*  we have a cubic spline curve. */
    k = 3;
    k1 = k + 1;
/*  set up the dimension information */
    nc = 26;
    mest = 20;
/*  n denotes the total number of knots. */
    n = 13;
/*  set up the knots of the spline curve */
    t[3] = 0.f;
    t[4] = .2f;
    t[5] = .3f;
    t[6] = .5f;
    t[7] = .6f;
    t[8] = .7f;
    t[9] = 1.f;
/*  fetch the b-spline coefficients for sx(u) */
    c__[0] = 1.f;
    c__[1] = 3.f;
    c__[2] = 4.f;
    c__[3] = 5.f;
    c__[4] = 3.f;
    c__[5] = -1.f;
/*  fetch the b-spline coefficients for sy(u) */
    c__[13] = 1.f;
    c__[14] = 2.f;
    c__[15] = -3.f;
    c__[16] = 2.f;
    c__[17] = 1.f;
    c__[18] = 4.f;
/*  we have a closed curve. */
/*  incorporate the boundary conditions for periodic splines */
    per = t[9] - t[3];
    for (i__ = 1; i__ <= 3; ++i__) {
/*  the boundary knots */
	t[i__ - 1] = t[i__ + 5] - per;
	t[i__ + 9] = t[i__ + 3] + per;
/*  the boundary coefficients */
	c__[i__ + 5] = c__[i__ - 1];
	c__[i__ + 18] = c__[i__ + 12];
/* L10: */
    }
/*  print the data of the spline curve. */
    s_wsfe(&io___11);
    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___12);
    e_wsfe();
    s_wsfe(&io___13);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___14);
    e_wsfe();
    nk1 = n - k1;
    s_wsfe(&io___16);
    i__1 = nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___17);
    e_wsfe();
    i1 = n + 1;
    i2 = n + nk1;
    s_wsfe(&io___20);
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
/*  loop for the different lines. */
    for (is = 1; is <= 5; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L150;
	}
/*  fetch the parameters of the straight line. */
L110:
	alfa = 0.f;
	beta = 1.f;
	gamma = 0.f;
	goto L160;
L120:
	alfa = 1.f;
	beta = 0.f;
	goto L160;
L130:
	beta = -1.f;
	goto L160;
L140:
	alfa = .4f;
	beta = .3f;
	gamma = 1.2f;
	goto L160;
L150:
	beta = .4f;
	gamma = 0.f;
/*  print the parameters of the straight line. */
L160:
	s_wsfe(&io___25);
	do_fio(&c__1, (char *)&alfa, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&beta, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&gamma, (ftnlen)sizeof(real));
	e_wsfe();
/*  calculate the coefficients of s(u) = sx(u)*alfa + sy(u)*beta - gamma */
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    j = i__ + n;
	    cc[i__ - 1] = alfa * c__[i__ - 1] + beta * c__[j - 1] - gamma;
/* L170: */
	}
/*  find the zeros of s(u) */
	sproot_(t, &n, cc, zero, &mest, &m, &ier);
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&m, (ftnlen)sizeof(integer));
	e_wsfe();
	if (m == 0) {
	    goto L200;
	}
/*  find the intersection points */
	curev_(&idim, t, &n, c__, &nc, &k, zero, &m, sp, &nc, &ier);
/*  print the intersection points */
	s_wsfe(&io___33);
	e_wsfe();
	l2 = 0;
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    l1 = l2 + 1;
	    l2 = l1 + 1;
	    s_wsfe(&io___36);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&zero[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[l2 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L180: */
	}
L200:
	;
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

