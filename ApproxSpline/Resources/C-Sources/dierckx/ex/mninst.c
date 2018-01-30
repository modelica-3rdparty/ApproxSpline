/* mninst.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static real c_b14 = 1.f;
static real c_b15 = 0.f;
static real c_b41 = .1f;
static real c_b42 = .8f;
static real c_b43 = .9f;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                 mninst : insert test program                       cc */
/* c      application : to find the sum of two periodic splines         cc */
/* c                 with different sets of knots                       cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0020insertion algorithm for ordinary splines\002)";
    static char fmt_905[] = "(\0020insertion algorithm for periodic splines\002)";
    static char fmt_910[] = "(1x,\002degree of the splines k =\002,i2)";
    static char fmt_915[] = "(1x,\002position of the knots of s1(x)\002)";
    static char fmt_920[] = "(5x,15f5.1)";
    static char fmt_925[] = "(1x,\002b-spline coefficients of s1(x)\002)";
    static char fmt_930[] = "(5x,8f9.5)";
    static char fmt_935[] = "(1x,\002position of the knots of s2(x)\002)";
    static char fmt_940[] = "(1x,\002b-spline coefficients of s2(x)\002)";
    static char fmt_945[] = "(1x,\002position of the knots after insertion\002)";
    static char fmt_950[] = "(1x,\002b-spline coefficients of s1+s2(x)\002)";
    static char fmt_955[] = "(\0020 i\002,6x,\002x\002,7x,\002s1(x)\002,7x,\002s2(x)\002,6x,\002s1+s2(x)\002)";
    static char fmt_960[] = "(1x,i2,f8.2,3f12.5)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer nest, iopt, i__, j, k, m;
    static real x[21], y[21], c1[30], c2[30];
    extern /* Subroutine */ int splev_(real *, integer *, real *, integer *, real *, real *, integer *, integer *);
    static integer i1, i2, j1, j2, k1, n1, n2;
    static real t1[30], t2[30], y1[21], y2[21], ai;
    static integer ip, nk;
    extern /* Subroutine */ int insert_(integer *, real *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
    static integer n1k1, n2k1, ier;
    static real per;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___11 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_960, 0 };


/*  set up the points where the splines will be evaluated. */
    m = 21;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 1);
	x[i__ - 1] = ai * .05f;
/* L10: */
    }
/*  set up the dimension information */
    nest = 30;
/*  main loop for the different spline degrees. */
    for (k = 3; k <= 5; k += 2) {
	k1 = k + 1;
	for (ip = 1; ip <= 2; ++ip) {
/*  if iopt = 1 the splines will be considered as periodic splines. */
/*  if iopt = 0 they will be considered as ordinary splines. */
	    iopt = 2 - ip;
	    if (iopt == 0) {
		s_wsfe(&io___10);
		e_wsfe();
	    }
	    if (iopt != 0) {
		s_wsfe(&io___11);
		e_wsfe();
	    }
	    s_wsfe(&io___12);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfe();
/*  fetch the knots and b-spline coefficients of the first spline s1(x). */
	    n1 = (k1 << 1) + 5;
	    n1k1 = n1 - k1;
	    t1[k1 - 1] = 0.f;
	    t1[k1] = .2f;
	    t1[k1 + 1] = .3f;
	    t1[k1 + 2] = .4f;
	    t1[k1 + 3] = .7f;
	    t1[k1 + 4] = .9f;
	    t1[k1 + 5] = 1.f;
	    c1[0] = 1.f;
	    c1[1] = 2.f;
	    c1[2] = -1.f;
	    c1[3] = 3.f;
	    c1[4] = 3.f;
	    c1[5] = -3.f;
/*  fetch the knots and b-spline coefficients of the second spline s2(x). */
	    n2 = (k1 << 1) + 6;
	    n2k1 = n2 - k1;
	    t2[k1 - 1] = 0.f;
	    t2[k1] = .1f;
	    t2[k1 + 1] = .2f;
	    t2[k1 + 2] = .3f;
	    t2[k1 + 3] = .4f;
	    t2[k1 + 4] = .7f;
	    t2[k1 + 5] = .8f;
	    t2[k1 + 6] = 1.f;
	    c2[0] = 2.f;
	    c2[1] = -2.f;
	    c2[2] = 1.f;
	    c2[3] = -3.f;
	    c2[4] = 4.f;
	    c2[5] = 4.f;
	    c2[6] = 4.f;
/*  incorporate the boundary conditions for periodic splines. */
	    per = 1.f;
	    nk = n1 - k;
	    i__1 = k;
	    for (j = 1; j <= i__1; ++j) {
		i1 = nk + j;
		i2 = nk - j;
		j1 = k1 + j;
		j2 = k1 - j;
/*  the boundary knots */
		t1[i1 - 1] = t1[j1 - 1] + per;
		t1[j2 - 1] = t1[i2 - 1] - per;
		t2[i1] = t2[j1 - 1] + per;
		t2[j2 - 1] = t2[i2] - per;
/*  the boundary coefficients */
		c1[j + 5] = c1[j - 1];
		c2[j + 6] = c2[j - 1];
/* L20: */
	    }
	    if (iopt != 0) {
		goto L100;
	    }
/*  if iopt=0 we insert k knots at the boundaries of the interval to */
/*  find the representation with coincident boundary knots */
	    i__1 = k;
	    for (j = 1; j <= i__1; ++j) {
		insert_(&iopt, t1, &n1, c1, &k, &c_b14, t1, &n1, c1, &nest, &ier);
		--n1;
		insert_(&iopt, t1, &n1, c1, &k, &c_b15, t1, &n1, c1, &nest, &ier);
		--n1;
		i__2 = n1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    t1[i__ - 1] = t1[i__];
		    c1[i__ - 1] = c1[i__];
/* L30: */
		}
/* L40: */
	    }
	    i__1 = k;
	    for (j = 1; j <= i__1; ++j) {
		insert_(&iopt, t2, &n2, c2, &k, &c_b14, t2, &n2, c2, &nest, &ier);
		--n2;
		insert_(&iopt, t2, &n2, c2, &k, &c_b15, t2, &n2, c2, &nest, &ier);
		--n2;
		i__2 = n2;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    t2[i__ - 1] = t2[i__];
		    c2[i__ - 1] = c2[i__];
/* L50: */
		}
/* L60: */
	    }
/*  print knots and b-spline coefficients of the two splines. */
L100:
	    s_wsfe(&io___29);
	    e_wsfe();
	    s_wsfe(&io___30);
	    i__1 = n1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t1[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___31);
	    e_wsfe();
	    s_wsfe(&io___32);
	    i__1 = n1k1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c1[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___33);
	    e_wsfe();
	    s_wsfe(&io___34);
	    i__1 = n2;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t2[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___35);
	    e_wsfe();
	    s_wsfe(&io___36);
	    i__1 = n2k1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c2[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/*  evaluate the two splines */
	    splev_(t1, &n1, c1, &k, x, y1, &m, &ier);
	    splev_(t2, &n2, c2, &k, x, y2, &m, &ier);
/*  insert the knots of the second spline into those of the first one */
	    insert_(&iopt, t1, &n1, c1, &k, &c_b41, t1, &n1, c1, &nest, &ier);
	    insert_(&iopt, t1, &n1, c1, &k, &c_b42, t1, &n1, c1, &nest, &ier);
/*  insert the knots of the first spline into those of the second one */
	    insert_(&iopt, t2, &n2, c2, &k, &c_b43, t2, &n2, c2, &nest, &ier);
/*  print the knots and coefficients of the splines in their new */
/*  representation */
	    n1k1 = n1 - k1;
	    s_wsfe(&io___39);
	    e_wsfe();
	    s_wsfe(&io___40);
	    i__1 = n1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t1[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___41);
	    e_wsfe();
	    s_wsfe(&io___42);
	    i__1 = n1k1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c1[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___43);
	    e_wsfe();
	    s_wsfe(&io___44);
	    i__1 = n1k1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c2[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/*  find the coefficients of the sum of the two splines. */
	    i__1 = n1k1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		c1[i__ - 1] += c2[i__ - 1];
/* L110: */
	    }
	    s_wsfe(&io___45);
	    e_wsfe();
	    s_wsfe(&io___46);
	    i__1 = n1k1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c1[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/*  evaluate this new spline and compare results */
	    splev_(t1, &n1, c1, &k, x, y, &m, &ier);
	    s_wsfe(&io___48);
	    e_wsfe();
	    i__1 = m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		s_wsfe(&io___49);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&x[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&y1[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&y2[i__ - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&y[i__ - 1], (ftnlen)sizeof(real));
		e_wsfe();
/* L200: */
	    }
/* L700: */
	}
/* L800: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

