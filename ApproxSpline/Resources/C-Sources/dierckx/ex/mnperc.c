/* mnperc.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                mnperc : percur test program                        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real x[27] = { 0.f,3.922f,7.843f,11.765f,15.686f,19.608f,23.509f,27.451f,31.373f,35.294f,39.216f,43.137f,47.059f,50.98f,54.902f,58.824f,62.745f,66.667f,70.588f,74.51f,78.431f,82.353f,86.275f,90.196f,94.118f,98.039f };
    static real y[27] = { 10.099f,14.835f,21.453f,25.022f,22.427f,22.315f,22.07f,19.673f,16.754f,13.983f,11.973f,12.286f,16.129f,21.56f,28.041f,39.205f,59.489f,72.559f,75.96f,79.137f,75.925f,68.809f,55.758f,39.915f,22.006f,12.076f };

    /* Format strings */
    static char fmt_910[] = "(\0020least-squares periodic spline of degree \002,i1)";
    static char fmt_915[] = "(\0020smoothing periodic spline of degree \002,i1)";
    static char fmt_920[] = "(\002 smoothing factor s=\002,f7.0)";
    static char fmt_925[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i2)";
    static char fmt_930[] = "(1x,\002total number of knots n=\002,i3)";
    static char fmt_935[] = "(1x,\002position of the knots \002)";
    static char fmt_940[] = "(5x,8f8.3)";
    static char fmt_945[] = "(\0020b-spline coefficients \002)";
    static char fmt_950[] = "(5x,8f8.4)";
    static char fmt_955[] = "(\0020\002,3(3x,\002xi\002,6x,\002yi\002,4x,\002s(xi)\002,3x))";
    static char fmt_960[] = "(\002 \002,3(f7.3,1x,f7.3,1x,f7.3,2x))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer nest, iopt, iwrk[37], lwrk;
    static real c__[37];
    static integer i__, j, k, l, m, n;
    static real s, t[37], w[27];
    extern /* Subroutine */ int splev_(real *, integer *, real *, integer *, real *, real *, integer *, integer *);
    static integer l1, l2, m1;
    static real al, fp;
    static integer is;
    static real sp[27];
    extern /* Subroutine */ int percur_(integer *, integer *, real *, real *, real *, integer *, real *, integer *, integer *, real *, real *, real *, real *, integer *, integer *, integer *);
    static integer nk1, ier;
    static real wrk[1400];

    /* Fortran I/O blocks */
    static cilist io___23 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___24 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___25 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___26 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___27 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_960, 0 };


/*  the data absciss values */
/*  the data ordinate values */
/*  m denotes the number of data points */
    m = 27;
/*  the period of the spline is determined by x(m) */
    x[m - 1] = 100.f;
    y[m - 1] = y[0];
/*  we set up the weights of the data points */
    m1 = m - 1;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 1.f;
/* L10: */
    }
/*  we set up the dimension information. */
    nest = 37;
    lwrk = 1400;
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
/*  we start computing the least-squares constant (large value for s). */
L110:
	    iopt = 0;
	    s = 6.5e4f;
	    goto L200;
/*  iopt=1 from the second call on */
L120:
	    iopt = 1;
	    s = 500.f;
	    goto L200;
/*  a smaller value for s to get a closer approximation */
L130:
	    s = 5.f;
	    goto L200;
/*  a larger value for s to get a smoother approximation */
L140:
	    s = 20.f;
	    goto L200;
/*  if a satisfactory fit is obtained  we can calculate a spline of equal */
/*  quality of fit ( same value for s ) but possibly with fewer knots by */
/*  specifying iopt=0 */
L150:
	    s = 20.f;
	    iopt = 0;
	    goto L200;
/*  we calculate an interpolating periodic spline. */
L160:
	    s = 0.f;
	    goto L200;
/*  finally, we also calculate a least-squares periodic spline function */
/*  with specified knots. */
L170:
	    iopt = -1;
	    n = (k << 1) + 11;
	    j = k + 2;
	    for (l = 1; l <= 9; ++l) {
		al = (real) (l * 10);
		t[j - 1] = al;
		++j;
/* L180: */
	    }
/*  determine the periodic spline approximation */
L200:
	    percur_(&iopt, &m, x, y, w, &k, &s, &nest, &n, t, c__, &fp, wrk, &lwrk, iwrk, &ier);
/*  printing of the results. */
	    if (iopt >= 0) {
		goto L210;
	    }
	    s_wsfe(&io___23);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfe();
	    goto L220;
L210:
	    s_wsfe(&io___24);
	    do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___25);
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	    e_wsfe();
L220:
	    s_wsfe(&io___26);
	    do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___27);
	    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___28);
	    e_wsfe();
	    s_wsfe(&io___29);
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    nk1 = n - k - 1;
	    s_wsfe(&io___31);
	    e_wsfe();
	    s_wsfe(&io___32);
	    i__1 = nk1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&c__[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	    s_wsfe(&io___33);
	    e_wsfe();
/*  evaluation of the spline approximation */
	    splev_(t, &n, c__, &k, x, sp, &m, &ier);
	    for (i__ = 1; i__ <= 9; ++i__) {
		l1 = (i__ - 1) * 3 + 1;
		l2 = l1 + 2;
		s_wsfe(&io___37);
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

