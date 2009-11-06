/* mnpasu.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__693 = 693;
static integer c__128 = 128;
static integer c__32 = 32;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c               mnpasu : parsur test program                         cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Format strings */
    static char fmt_900[] = "(\0021the input data\002)";
    static char fmt_905[] = "(\0020\002,2x,\002v\002,11(3x,f4.1))";
    static char fmt_910[] = "(\002 \002,1x,\002u\002)";
    static char fmt_915[] = "(\002 \002,f4.1)";
    static char fmt_920[] = "(11f7.3)";
    static char fmt_925[] = "(5x,11f7.3)";
    static char fmt_935[] = "(\0020least-squares surface of periodicity\002,2i3)";
    static char fmt_940[] = "(\0020smoothing surface of periodicity\002,2i3)";
    static char fmt_945[] = "(\002 smoothing factor s=\002,f8.2)";
    static char fmt_950[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i3)";
    static char fmt_955[] = "(1x,\002total number of knots in the u-direction =\002,i3)";
    static char fmt_960[] = "(1x,\002position of the knots \002)";
    static char fmt_965[] = "(5x,10f6.2)";
    static char fmt_970[] = "(1x,\002total number of knots in the v-direction =\002,i3)";
    static char fmt_975[] = "(\0020b-spline coefficients \002)";
    static char fmt_980[] = "(5x,8f9.4)";
    static char fmt_985[] = "(\0020\002,\002spline values at selected grid points\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen), s_rsfe(cilist *), e_rsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer idim, ipar[2], iopt, iwrk[80], kwrk, lwrk;
    static real c__[900], f[693];
    static integer i__, j, l, m;
    static real s, u[21], v[11], z__[693];
    static integer j0, j1, j2, j3, nuest, nvest;
    extern /* Subroutine */ int surev_(integer *, real *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, integer *);
    static real ai;
    static integer nc;
    static real fp;
    static integer is, iw[32];
    static real wk[128];
    static integer mu, mv, nu, nv;
    static real tu[27], tv[17];
    extern /* Subroutine */ int parsur_(integer *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, integer *, integer *, integer *, integer *);
    static integer ier;
    static real wrk[2000];

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___8 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___9 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___17 = { 0, 5, 0, fmt_920, 0 };
    static cilist io___20 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_970, 0 };
    static cilist io___47 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_965, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_975, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_980, 0 };
    static cilist io___55 = { 0, 6, 0, fmt_985, 0 };
    static cilist io___56 = { 0, 6, 0, fmt_905, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___58 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___59 = { 0, 6, 0, fmt_925, 0 };


/*  we generate the u-coordinates of the grid. */
    mu = 21;
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u[i__ - 1] = (real) (i__ - 1);
/* L10: */
    }
/*  we generate the v-coordinates of the grid. */
    mv = 11;
    i__1 = mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ - 1] = u[(i__ << 1) - 2];
/* L20: */
    }
/*  the dimension of the surface */
    idim = 3;
/*  we fetch and print the surface co-ordinates at the grid points */
    s_wsfe(&io___7);
    e_wsfe();
    s_wsfe(&io___8);
    i__1 = mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&v[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
    s_wsfe(&io___9);
    e_wsfe();
    m = mu * mv;
    j0 = 0;
    i__1 = mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___12);
	do_fio(&c__1, (char *)&u[i__ - 1], (ftnlen)sizeof(real));
	e_wsfe();
	j1 = j0;
	i__2 = idim;
	for (l = 1; l <= i__2; ++l) {
	    j2 = j1 + 1;
	    j3 = j1 + mv;
	    s_rsfe(&io___17);
	    i__3 = j3;
	    for (j = j2; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&f[j - 1], (ftnlen)sizeof(real));
	    }
	    e_rsfe();
	    s_wsfe(&io___20);
	    i__3 = j3;
	    for (j = j2; j <= i__3; ++j) {
		do_fio(&c__1, (char *)&f[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    j1 += m;
/* L30: */
	}
	j0 += mv;
/* L40: */
    }
/*  we set up the dimension information */
    nuest = 27;
    nvest = 17;
    lwrk = 2000;
    kwrk = 80;
/*  main loop for the different spline approximations */
    for (is = 1; is <= 4; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	}
/*  a smoothing surface with no periodicity conditions */
L110:
	iopt = 0;
	s = .07f;
	ipar[0] = 0;
	ipar[1] = 0;
	goto L200;
/*  a smoothing surface periodic in the v-variable */
L120:
	ipar[1] = 1;
	goto L200;
/*  a smoothing surface periodic in both variables */
L130:
	ipar[0] = 1;
	goto L200;
/*  finally we also calculate a least-squares spline surface */
/*  with specified knots. */
L140:
	iopt = -1;
	nu = 11;
	nv = 11;
	j = 5;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ai = (real) (i__ * 5);
	    tu[j - 1] = ai;
	    tv[j - 1] = tu[j - 1];
	    ++j;
/* L150: */
	}
/*  determination of the spline surface. */
L200:
	parsur_(&iopt, ipar, &idim, &mu, u, &mv, v, f, &s, &nuest, &nvest, &nu, tu, &nv, tv, c__, &fp, wrk, &lwrk, iwrk, &kwrk, &ier);
/*  printing of the fitting results. */
	if (iopt >= 0) {
	    goto L210;
	}
	s_wsfe(&io___39);
	do_fio(&c__1, (char *)&ipar[0], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar[1], (ftnlen)sizeof(integer));
	e_wsfe();
	goto L220;
L210:
	s_wsfe(&io___40);
	do_fio(&c__1, (char *)&ipar[0], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar[1], (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L220:
	s_wsfe(&io___42);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___43);
	do_fio(&c__1, (char *)&nu, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___44);
	e_wsfe();
	s_wsfe(&io___45);
	i__1 = nu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tu[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___46);
	do_fio(&c__1, (char *)&nv, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___47);
	e_wsfe();
	s_wsfe(&io___48);
	i__1 = nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&tv[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nc = (nu - 4) * (nv - 4);
	s_wsfe(&io___50);
	e_wsfe();
	j1 = 0;
	i__1 = idim;
	for (l = 1; l <= i__1; ++l) {
	    j0 = j1 + 1;
	    j1 += nc;
	    s_wsfe(&io___51);
	    i__2 = j1;
	    for (j = j0; j <= i__2; ++j) {
		do_fio(&c__1, (char *)&c__[j - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
/* L230: */
	}
/*  evaluation of the spline surface. */
	surev_(&idim, tu, &nu, tv, &nv, c__, u, &mu, v, &mv, z__, &c__693, wk, &c__128, iw, &c__32, &ier);
	s_wsfe(&io___55);
	e_wsfe();
	s_wsfe(&io___56);
	i__1 = mv;
	for (i__ = 1; i__ <= i__1; i__ += 2) {
	    do_fio(&c__1, (char *)&v[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___57);
	e_wsfe();
	j0 = 0;
	i__1 = mu;
	for (i__ = 1; i__ <= i__1; i__ += 4) {
	    s_wsfe(&io___58);
	    do_fio(&c__1, (char *)&u[i__ - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	    j1 = j0;
	    i__2 = idim;
	    for (l = 1; l <= i__2; ++l) {
		j2 = j1 + 1;
		j3 = j1 + mv;
		s_wsfe(&io___59);
		i__3 = j3;
		for (j = j2; j <= i__3; j += 2) {
		    do_fio(&c__1, (char *)&z__[j - 1], (ftnlen)sizeof(real));
		}
		e_wsfe();
		j1 += m;
/* L240: */
	    }
	    j0 += mv << 2;
/* L250: */
	}
/* L300: */
    }
    s_stop("", 0L);
/*  format statements. */
    return 0;
} /* MAIN__ */

