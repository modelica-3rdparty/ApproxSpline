/* mncual.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mncual : cualde test program                       cc */
/* c         evaluation of a closed planar spline curve                 cc */
/* c                    x = sx(u) , y = sy(u)                           cc */
/* c            through its polynomial representation                   cc */
/* c                    in each knot interval.                          cc */
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
    static char fmt_930[] = "(\0020knot interval (\002,f4.1,\002,\002,f4.1,\002)\002)";
    static char fmt_935[] = "(1x,\002curve derivatives at the midpoint of the interval\002)";
    static char fmt_940[] = "(1x,3(1x,2e12.4))";
    static char fmt_945[] = "(1x,\002coefficients in the polynomial represent. of sx(u)\002)";
    static char fmt_950[] = "(2x,6e13.5)";
    static char fmt_955[] = "(1x,\002coefficients in the polynomial represent. of sy(u)\002)";
    static char fmt_960[] = "(1x,2(5x,\002i\002,3x,\002u(i)\002,4x,\002sx(u(i))\002,4x,\002sy(u(i))\002))";
    static char fmt_965[] = "(1x,2(i6,f7.2,2f12.5))";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer idim;
    static real c__[40], d__[12];
    static integer i__, j, k, l, m, n;
    static real t[20], u[20];
    static integer i1, i2, j1, j2, j3, j4, k1, l1;
    static real ai, aj;
    static integer nc, ii, nd, jj, kk, jn, ip, nk;
    extern /* Subroutine */ int cualde_(integer *, real *, integer *, real *, integer *, integer *, real *, real *, integer *, integer *);
    static real sp[40], tt, uu;
    static integer nk1;
    static real fac, cof[12]	/* was [2][6] */, arg;
    static integer ier;
    static real per, pol;

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___23 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_965, 0 };



#define cof_ref(a_1,a_2) cof[(a_2)*2 + a_1 - 3]

/*  we have a planar curve */
    idim = 2;
/*  set up the dimension information */
    nc = 40;
    nd = 12;
/*  set up the points where the curve will be evaluated. */
    m = 20;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 1);
	u[i__ - 1] = ai * .05f;
/* L10: */
    }
/*  main loop for the different spline degrees. */
    for (k = 3; k <= 5; k += 2) {
/*  the order of the spline. */
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
/*  print the data for the spline. */
	s_wsfe(&io___21);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___22);
	e_wsfe();
	s_wsfe(&io___23);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___24);
	e_wsfe();
	nk1 = n - k1;
	s_wsfe(&io___26);
	i__1 = nk1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___27);
	e_wsfe();
	i1 = n + 1;
	i2 = n + nk1;
	s_wsfe(&io___28);
	i__1 = i2;
	for (i__ = i1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	l = k;
	l1 = k1;
	kk = k1 * idim;
/*  main loop for the different points of evaluation. */
	ip = 0;
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    arg = u[i__ - 1];
/*  search for knot interval t(l)<=u(i)<t(l+1). */
L40:
	    if (arg < t[l1 - 1] || l == nk1) {
		goto L70;
	    }
/*  a new knot interval. */
	    l = l1;
	    l1 = l + 1;
	    if (t[l - 1] == t[l1 - 1]) {
		goto L40;
	    }
	    s_wsfe(&io___34);
	    do_fio(&c__1, (char *)&t[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&t[l1 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/*  calculate the spline derivatives at the midpoint tt of the interval */
	    tt = (t[l - 1] + t[l1 - 1]) * .5f;
	    cualde_(&idim, t, &n, c__, &nc, &k1, &tt, d__, &nd, &ier);
	    s_wsfe(&io___38);
	    e_wsfe();
	    s_wsfe(&io___39);
	    i__2 = kk;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&d__[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/*  calculate the coefficients cof in the polynomial representation of */
/*  the spline curve in the current knot interval,i.e. */
/*    sx(u) = cof(1,1)+cof(1,2)*(u-tt)+...+cof(1,k1)*(u-tt)**k */
/*    sy(u) = cof(2,1)+cof(2,2)*(u-tt)+...+cof(2,k1)*(u-tt)**k */
	    fac = 1.f;
	    jj = 0;
	    i__2 = k1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = idim;
		for (ii = 1; ii <= i__3; ++ii) {
		    ++jj;
		    cof_ref(ii, j) = d__[jj - 1] / fac;
/* L50: */
		}
		aj = (real) j;
		fac *= aj;
/* L60: */
	    }
	    s_wsfe(&io___45);
	    e_wsfe();
	    s_wsfe(&io___46);
	    i__2 = k1;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&cof_ref(1, j), (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___47);
	    e_wsfe();
	    s_wsfe(&io___48);
	    i__2 = k1;
	    for (j = 1; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&cof_ref(2, j), (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    goto L40;
/*  evaluate the polynomial curve */
L70:
	    uu = arg - tt;
	    i__2 = idim;
	    for (ii = 1; ii <= i__2; ++ii) {
		pol = cof_ref(ii, k1);
		jj = k1;
		i__3 = k;
		for (j = 1; j <= i__3; ++j) {
		    --jj;
		    pol = pol * uu + cof_ref(ii, jj);
/* L80: */
		}
		++ip;
		sp[ip - 1] = pol;
/* L90: */
	    }
/* L100: */
	}
	s_wsfe(&io___52);
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
	    s_wsfe(&io___55);
	    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u[i1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j2 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&u[i2 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j3 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j4 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L110: */
	}
/* L120: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

#undef cof_ref


