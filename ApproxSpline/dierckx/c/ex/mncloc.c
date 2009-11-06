/* mncloc.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                mncloc : clocur test program                        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real x[38] = { -4.7f,0.f,-7.048f,2.565f,-6.894f,5.785f,-3.75f,6.495f,-1.042f,5.909f,.938f,5.318f,2.5f,4.33f,3.524f,2.957f,4.511f,1.642f,5.f,0.f,4.886f,-1.779f,3.524f,-2.957f,3.2f,-5.543f,1.302f,-7.386f,-1.424f,-8.075f,-3.f,-5.196f,-3.064f,-2.571f,-3.665f,-1.334f };

    /* Format strings */
    static char fmt_910[] = "(\0020least-squares closed curve of degree \002,i1,\002  ipar=\002,i1)";
    static char fmt_915[] = "(\0020smoothing closed curve of degree \002,i1,\002  ipar=\002,i1)";
    static char fmt_920[] = "(\002 smoothing factor s=\002,f7.1)";
    static char fmt_925[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i2)";
    static char fmt_930[] = "(1x,\002total number of knots n=\002,i3)";
    static char fmt_935[] = "(1x,\002position of the knots \002)";
    static char fmt_940[] = "(5x,10f6.0)";
    static char fmt_950[] = "(5x,8f9.4)";
    static char fmt_945[] = "(1x,\002b-spline coefficients of sx(u)\002)";
    static char fmt_955[] = "(1x,\002b-spline coefficients of sy(u)\002)";
    static char fmt_960[] = "(\0020\002,2(4x,\002xi\002,7x,\002yi\002,6x,\002sx(ui)\002,3x,\002sy(ui)\002))";
    static char fmt_965[] = "(\002 \002,8f9.4)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer idim, ipar, nest, iopt, iwrk[40], lwrk;
    static real c__[80];
    static integer i__, j, k, l, m, n;
    static real s, t[40], u[19], w[19];
    extern /* Subroutine */ int curev_(integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, integer *);
    static integer i1, i2, j1, l1;
    static real al;
    static integer nc;
    static real fp;
    static integer is;
    static real sp[40];
    static integer mx;
    extern /* Subroutine */ int clocur_(integer *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, real *, real *, integer *, integer *, integer *);
    static integer nk1;
    static real del;
    static integer ier;
    static real wrk[1500];

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___46 = { 0, 6, 0, fmt_965, 0 };


/*  the data absciss values */
/*  the data ordinate values */
/*  m denotes the number of data points */
    m = 19;
/*  the first and last data point coincide */
    x[(m << 1) - 2] = x[0];
    x[(m << 1) - 1] = x[1];
/*  we set up the weights and parameter values of the data points */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 1.f;
	al = (real) ((i__ - 1) * 20);
	u[i__ - 1] = al;
/* L10: */
    }
/*  we set up the dimension information. */
    nest = 40;
    lwrk = 1500;
    nc = 80;
    mx = 38;
/*  we will determine a planar closed curve   x=sx(u) , y=sy(u) */
    idim = 2;
/*  for the first approximations we will use cubic splines */
    k = 3;
/*  we will also supply the parameter values u(i) */
    ipar = 1;
/*  loop for the different approximating spline curves */
    for (is = 1; is <= 9; ++is) {
	switch (is) {
	    case 1:  goto L110;
	    case 2:  goto L120;
	    case 3:  goto L130;
	    case 4:  goto L140;
	    case 5:  goto L150;
	    case 6:  goto L160;
	    case 7:  goto L170;
	    case 8:  goto L180;
	    case 9:  goto L190;
	}
/*  we start computing the least-squares point ( s very large) */
L110:
	iopt = 0;
	s = 900.f;
	goto L300;
/*  iopt =  1 from the second call on */
L120:
	iopt = 1;
	s = 10.f;
	goto L300;
/*  a smaller value for s to get a closer approximation */
L130:
	s = .1f;
	goto L300;
/*  a larger value for s to get a smoother approximation */
L140:
	s = .5f;
	goto L300;
/*  if a satisfactory fit is obtained we can calculate a curve of equal */
/*  quality of fit (same value for s) but possibly with fewer knots by */
/*  specifying iopt=0 */
L150:
	iopt = 0;
	s = .5f;
	goto L300;
/*  we determine a spline curve with respect to the same smoothing */
/*  factor s,  but now we let the program determine parameter values u(i) */
L160:
	ipar = 0;
	iopt = 0;
	s = .5f;
	goto L300;
/*  we choose a different degree of spline approximation */
L170:
	k = 5;
	iopt = 0;
	s = .5f;
	goto L300;
/*  we determine an interpolating curve */
L180:
	s = 0.f;
	goto L300;
/*  finally we calculate a least-squares spline curve with specified */
/*  knots */
L190:
	iopt = -1;
	n = (k << 1) + 9;
	j = k + 2;
	del = (u[m - 1] - u[0]) * .125f;
	for (l = 1; l <= 7; ++l) {
	    al = (real) l;
	    t[j - 1] = u[0] + al * del;
	    ++j;
/* L200: */
	}
/*  determine the approximating closed curve */
L300:
	clocur_(&iopt, &ipar, &idim, &m, u, &mx, x, w, &k, &s, &nest, &n, t, &nc, c__, &fp, wrk, &lwrk, iwrk, &ier);
/*  printing of the results. */
	if (iopt >= 0) {
	    goto L310;
	}
	s_wsfe(&io___27);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L320;
L310:
	s_wsfe(&io___28);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___29);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L320:
	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___32);
	e_wsfe();
	if (ipar == 1) {
	    s_wsfe(&io___33);
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	if (ipar == 0) {
	    s_wsfe(&io___34);
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	nk1 = n - k - 1;
	s_wsfe(&io___36);
	e_wsfe();
	s_wsfe(&io___37);
	i__1 = nk1;
	for (l = 1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&c__[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___38);
	e_wsfe();
	i1 = n + 1;
	i2 = n + nk1;
	s_wsfe(&io___41);
	i__1 = i2;
	for (l = i1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&c__[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___42);
	e_wsfe();
/*  we evaluate the spline curve */
	curev_(&idim, t, &n, c__, &nc, &k, u, &m, sp, &mx, &ier);
	for (i__ = 1; i__ <= 9; ++i__) {
	    l = (i__ - 1 << 2) + 1;
	    l1 = l + 1;
	    j = l + 2;
	    j1 = j + 1;
	    s_wsfe(&io___46);
	    do_fio(&c__1, (char *)&x[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&x[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[l - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[l1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&x[j - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&x[j1 - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&sp[j1 - 1], (ftnlen)sizeof(real));
	    e_wsfe();
/* L330: */
	}
/* L400: */
    }
    s_stop("", 0L);
    return 0;
} /* MAIN__ */

