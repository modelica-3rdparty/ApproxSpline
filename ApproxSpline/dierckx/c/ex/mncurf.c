/* mncurf.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                mncurf : curfit test program                        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real y[25] = { 1.f,1.f,1.4f,1.1f,1.f,1.f,4.f,9.f,13.f,13.4f,12.8f,13.1f,13.f,14.f,13.f,13.5f,10.f,2.f,3.f,2.5f,2.5f,2.5f,3.f,4.f,3.5f };

    /* Format strings */
    static char fmt_910[] = "(\0020least-squares spline of degree \002,i1)";
    static char fmt_915[] = "(\0020smoothing spline of degree \002,i1)";
    static char fmt_920[] = "(\002 smoothing factor s=\002,f5.0)";
    static char fmt_925[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i2)";
    static char fmt_930[] = "(1x,\002total number of knots n=\002,i3)";
    static char fmt_935[] = "(1x,\002position of the knots \002)";
    static char fmt_940[] = "(5x,12f6.1)";
    static char fmt_945[] = "(\0020b-spline coefficients \002)";
    static char fmt_950[] = "(5x,8f9.4)";
    static char fmt_955[] = "(\0020\002,5(1x,\002xi\002,3x,\002yi\002,2x,\002s(xi)\002,1x))";
    static char fmt_960[] = "(\002 \002,5(f4.1,1x,f4.1,1x,f4.1,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer nest, iopt, iwrk[35], lwrk;
    static real c__[35];
    static integer i__, j, k, l, m, n;
    static real s, t[35], w[25], x[25];
    extern /* Subroutine */ int splev_(real *, integer *, real *, integer *, real *, real *, integer *, integer *);
    static integer l1, l2;
    static real ai, fp, xb;
    static integer is;
    static real xe, sp[25];
    extern /* Subroutine */ int curfit_(integer *, integer *, real *, real *, real *, real *, real *, integer *, real *, integer *, integer *, real *, real *, real *, real *, integer *, integer *, integer *);
    static integer nk1, ier;
    static real wrk[1000];

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_960, 0 };


/*  the ordinate values of the data points */
/*  m denotes the number of data points */
    m = 25;
/*  we set up the abscissae and weights of the data points */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ai = (real) (i__ - 1);
	x[i__ - 1] = ai;
	w[i__ - 1] = 1.f;
/* L10: */
    }
/*  we set up the boundaries of the approximation interval */
    xb = x[0];
    xe = x[m - 1];
/*  we set up the dimension information. */
    nest = 35;
    lwrk = 1000;
/*  loop for the different spline degrees. */
    for (k = 3; k <= 5; k += 2) {
/*  loop for the different spline approximations of degree k */
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
/*  we start computing the least-squares polynomial (large value for s). */
L110:
	    iopt = 0;
	    s = 1e3f;
	    goto L200;
/*  iopt=1 from the second call on */
L120:
	    iopt = 1;
	    s = 60.f;
	    goto L200;
/*  a smaller value for s to get a closer approximation */
L130:
	    s = 10.f;
	    goto L200;
/*  a larger value for s to get a smoother approximation */
L140:
	    s = 30.f;
	    goto L200;
/*  if a satisfactory fit is obtained  we can calculate a spline of equal */
/*  quality of fit ( same value for s ) but possibly with fewer knots by */
/*  specifying iopt=0 */
L150:
	    s = 30.f;
	    iopt = 0;
	    goto L200;
/*  we calculate an interpolating spline */
L160:
	    s = 0.f;
	    goto L200;
/*  finally, we also calculate a least-squares spline function with */
/*  specified knots */
L170:
	    iopt = -1;
	    j = k + 2;
	    for (l = 1; l <= 7; ++l) {
		ai = (real) (l * 3);
		t[j - 1] = ai;
		++j;
/* L180: */
	    }
	    n = (k << 1) + 9;
L200:
	    curfit_(&iopt, &m, x, y, w, &xb, &xe, &k, &s, &nest, &n, t, c__, &fp, wrk, &lwrk, iwrk, &ier);
/*  printing of the results. */
	    if (iopt >= 0) {
		goto L210;
	    }
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfe();
	    goto L220;
L210:
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	    e_wsfe();
L220:
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___28);
	    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___29);
	    e_wsfe();
	    s_wsfe(&io___30);
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    nk1 = n - k - 1;
	    s_wsfe(&io___32);
	    e_wsfe();
	    s_wsfe(&io___33);
	    i__1 = nk1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___34);
	    e_wsfe();
/*  evaluation of the spline approximation */
	    splev_(t, &n, c__, &k, x, sp, &m, &ier);
	    for (i__ = 1; i__ <= 5; ++i__) {
		l1 = (i__ - 1) * 5 + 1;
		l2 = l1 + 4;
		s_wsfe(&io___38);
		i__1 = l2;
		for (l = l1; l <= i__1; ++l) {
		    do_fio(&c__1, (char *)&x[l - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&y[l - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&sp[l - 1], (ftnlen)sizeof(real));
		}
		e_wsfe();
/* L230: */
	    }
/* L300: */
	}
/* L400: */
    }
    s_stop("", 0L);
    return 0;
} /* MAIN__ */

