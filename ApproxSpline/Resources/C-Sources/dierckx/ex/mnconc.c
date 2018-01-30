/* mnconc.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                mnconc : concur test program                        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real x[62] = { -3.109f,3.04f,-2.188f,2.876f,-1.351f,2.634f,-.605f,2.183f,.093f,1.586f,.451f,1.01f,.652f,.382f,.701f,-.218f,.518f,-.632f,.277f,-.879f,.008f,-.981f,-.291f,-.886f,-.562f,-.642f,-.679f,-.195f,-.637f,.373f,-.425f,1.07f,-.049f,1.607f,.575f,2.165f,1.334f,2.618f,2.167f,2.905f,3.206f,2.991f,4.099f,2.897f,4.872f,2.615f,5.71f,2.164f,6.33f,1.617f,6.741f,.977f,6.928f,.383f,6.965f,-.194f,6.842f,-.665f,6.593f,-.901f,6.269f,-1.01f };

    /* Format strings */
    static char fmt_910[] = "(\0020least-squares curve of degree \002,i1)";
    static char fmt_915[] = "(\0020smoothing curve of degree \002,i1)";
    static char fmt_920[] = "(\002 smoothing factor s=\002,f5.0)";
    static char fmt_925[] = "(\002 number of derivative constraints ib=\002,i2,5x,\002ie=\002,i2)";
    static char fmt_930[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i2)";
    static char fmt_935[] = "(1x,\002total number of knots n=\002,i3)";
    static char fmt_940[] = "(1x,\002position of the knots \002)";
    static char fmt_945[] = "(5x,8f8.4)";
    static char fmt_950[] = "(1x,\002b-spline coefficients of sx(u)\002)";
    static char fmt_955[] = "(1x,\002b-spline coefficients of sy(u)\002)";
    static char fmt_960[] = "(1x,\002derivatives at the begin point\002)";
    static char fmt_970[] = "(5x,\002order=\002,i2,2f9.4)";
    static char fmt_965[] = "(1x,\002derivatives at the end point\002)";
    static char fmt_975[] = "(\0020\002,2(4x,\002xi\002,7x,\002yi\002,6x,\002sx(ui)\002,3x,\002sy(ui)\002))";
    static char fmt_980[] = "(\002 \002,8f9.4)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer idim, nest, iopt, iwrk[50], lwrk;
    static real c__[100];
    static integer i__, j, k, l, m, n;
    static real s, t[50], u[31], w[31], sigma;
    extern /* Subroutine */ int curev_(integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, integer *);
    static integer i1, i2, j1, k1, l1, l2;
    static real db[6], dd[12], de[6], ai;
    static integer ib, ie, nb, nc;
    static real cp[24];
    static integer ne;
    static real fp;
    static integer kk;
    static real pi;
    static integer is;
    extern /* Subroutine */ int cualde_(integer *, real *, integer *, real *, integer *, integer *, real *, real *, integer *, integer *);
    static integer np;
    static real sp[62];
    static integer mx;
    static real ww, xx[62];
    extern /* Subroutine */ int concur_(integer *, integer *, integer *, real *, integer *, real *, real *, real *, integer *, real *, integer *, integer *, real *, integer *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, integer *, integer *, integer *);
    static integer nk1;
    static real del;
    static integer ndd, ier;
    static real wrk[1400];

    /* Fortran I/O blocks */
    static cilist io___39 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___60 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___61 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___62 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_980, 0 };


/*  the data absciss values */
/*  the data ordinate values */
/*  m denotes the number of data points */
    m = 31;
/*  set up the parameter values for the data points */
    pi = atan(1.f) * 4.f;
    del = pi * .1f;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 11);
	u[i__ - 1] = ai * del;
/* L10: */
    }
/*  the weights are taken as 1/sigma with sigma an estimate of the */
/*  standard deviation of the data points. */
    sigma = .04f;
    ww = 1.f / sigma;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = ww;
/* L20: */
    }
/*  accordingly, the smoothing factor is chosen as s = m */
    s = (real) m;
/*  we have a planar curve  x = sx(u) , y = sy(u) */
    idim = 2;
/*  begin point derivatives of the curve */
    db[0] = -pi;
    db[1] = 3.f;
    db[2] = 3.f;
    db[3] = 0.f;
    db[4] = 0.f;
    db[5] = -2.f;
/*  end point derivatives of the curve */
    de[0] = pi * 2.f;
    de[1] = -1.f;
    de[2] = -1.f;
    de[3] = 0.f;
    de[4] = 0.f;
    de[5] = 2.f;
/*  we set up the dimension information. */
    ndd = 12;
    np = 24;
    nb = 6;
    ne = 6;
    nest = 50;
    lwrk = 1400;
    nc = 100;
    mx = 62;
/*  for the first approximations we will use cubic splines. */
    k = 3;
/*  loop for the different spline curves */
    for (is = 1; is <= 7; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L150;
	    case 6:  goto L160;
	    case 7:  goto L170;
	}
/*  no derivative constraints */
L110:
	iopt = 0;
	ib = 0;
	ie = 0;
	goto L200;
/*  fixed end points. */
L120:
	iopt = 0;
	ib = 1;
	ie = 1;
	goto L200;
/*  first derivative constraint at the end point */
L130:
	iopt = 0;
	ib = 2;
	ie = 1;
	goto L200;
/*  first derivative constraints at begin and end point. */
L140:
	iopt = 0;
	ib = 2;
	ie = 2;
	goto L200;
/*  we choose quintic splines with second derivative constraints. */
L150:
	iopt = 0;
	k = 5;
	ib = 3;
	ie = 3;
	goto L200;
/*  we choose another s-value and continue with the set of knots found at */
/*  the last call of concur. */
L160:
	iopt = 1;
	s = 26.f;
	goto L200;
/*  finally we also calculate a least-squares curve with specified knots */
L170:
	iopt = -1;
	j = k + 2;
	for (l = 1; l <= 5; ++l) {
	    ai = (real) (l - 2);
	    t[j - 1] = ai * pi * .5f;
	    ++j;
/* L180: */
	}
	n = (k << 1) + 7;
/*  determination of the spline curve. */
L200:
	concur_(&iopt, &idim, &m, u, &mx, x, xx, w, &ib, db, &nb, &ie, de, &ne, &k, &s, &nest, &n, t, &nc, c__, &np, cp, &fp, wrk, &lwrk, iwrk, &ier);
/*  printing of the results. */
	if (iopt >= 0) {
	    goto L210;
	}
	s_wsfe(&io___39);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___40);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___42);
	do_fio(&c__1, (char *)&ib, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___43);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___44);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___45);
	e_wsfe();
	s_wsfe(&io___46);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nk1 = n - k - 1;
	s_wsfe(&io___48);
	e_wsfe();
	s_wsfe(&io___49);
	i__1 = nk1;
	for (l = 1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&c__[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___50);
	e_wsfe();
	i1 = n + 1;
	i2 = n + nk1;
	s_wsfe(&io___53);
	i__1 = i2;
	for (l = i1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&c__[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/*  calculate derivatives at the begin point. */
	k1 = k + 1;
	kk = k1 / 2;
	cualde_(&idim, t, &n, c__, &nc, &k1, u, dd, &ndd, &ier);
	s_wsfe(&io___57);
	e_wsfe();
	i__1 = kk;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    l = i__ - 1;
	    l1 = l * idim + 1;
	    l2 = l1 + 1;
	    s_wsfe(&io___60);
	    do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dd[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dd[l2 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L300: */
	}
/*  calculate derivatives at the end point. */
	cualde_(&idim, t, &n, c__, &nc, &k1, &u[m - 1], dd, &ndd, &ier);
	s_wsfe(&io___61);
	e_wsfe();
	i__1 = kk;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    l = i__ - 1;
	    l1 = l * idim + 1;
	    l2 = l1 + 1;
	    s_wsfe(&io___62);
	    do_fio(&c__1, (char *)&l, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&dd[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&dd[l2 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L350: */
	}
/*  we evaluate the spline curve */
	curev_(&idim, t, &n, c__, &nc, &k, u, &m, sp, &mx, &ier);
	s_wsfe(&io___64);
	e_wsfe();
	for (i__ = 1; i__ <= 5; ++i__) {
	    l = (i__ - 1) * 12 + 3;
	    l1 = l + 1;
	    j = l + 6;
	    j1 = j + 1;
	    s_wsfe(&io___66);
	    do_fio(&c__1, (char *)&x[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&x[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&x[j - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&x[j1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j1 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L400: */
	}
/* L500: */
    }
    s_stop("", 0L);
    return 0;
} /* MAIN__ */

