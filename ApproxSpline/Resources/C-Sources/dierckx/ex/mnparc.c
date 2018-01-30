/* mnparc.f -- translated by f2c (version 19961017).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* c                                                                    cc */
/* c                mnparc : parcur test program                        cc */
/* c                                                                    cc */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Main program */ MAIN__(void)
{
    /* Initialized data */

    static real u[32] = { 120.f,128.f,133.f,136.f,138.f,141.f,144.f,146.f,149.f,151.f,154.f,161.f,170.f,180.f,190.f,200.f,210.f,220.f,230.f,240.f,250.f,262.f,269.f,273.f,278.f,282.f,287.f,291.f,295.f,299.f,305.f,315.f };
    static real x[64] = { -1.5141f,.515f,-2.0906f,1.3412f,-1.9253f,2.6094f,-.8724f,3.2358f,-.3074f,2.7401f,-.5534f,2.7823f,.0192f,3.5932f,1.2298f,3.8353f,2.5479f,2.5863f,2.471f,1.3105f,1.7063f,.6841f,1.1183f,.2575f,.5534f,.246f,.4727f,.3689f,.3574f,.246f,.1998f,.2998f,.2882f,.3651f,.2613f,.3343f,.2652f,.3881f,.2805f,.4573f,.4112f,.5918f,.9377f,.711f,1.3527f,.4035f,1.5564f,.0769f,1.6141f,-.392f,1.6333f,-.857f,1.1567f,-1.3412f,.8109f,-1.5641f,.2498f,-1.7409f,-.2306f,-1.7178f,-.7571f,-1.2989f,
	    -1.1222f,-.5572f };

    /* Format strings */
    static char fmt_910[] = "(\0020least-squares curve of degree \002,i1,\002  ipar=\002,i1)";
    static char fmt_915[] = "(\0020smoothing curve of degree \002,i1,\002  ipar=\002,i1)";
    static char fmt_920[] = "(\002 smoothing factor s=\002,f7.2)";
    static char fmt_925[] = "(1x,\002sum squared residuals =\002,e15.6,5x,\002error flag=\002,i2)";
    static char fmt_930[] = "(1x,\002total number of knots n=\002,i3)";
    static char fmt_935[] = "(1x,\002position of the knots \002)";
    static char fmt_940[] = "(5x,10f6.0)";
    static char fmt_950[] = "(5x,8f8.4)";
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
    static real s, t[40], w[32];
    extern /* Subroutine */ int curev_(integer *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer *, integer *);
    static integer i1, i2, j1, l1;
    static real al;
    static integer nc;
    static real fp, ub, ue;
    static integer is;
    static real sp[64];
    static integer mx;
    extern /* Subroutine */ int parcur_(integer *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, real *, real *, integer *, integer *, integer *);
    static integer nk1;
    static real del;
    static integer ier;
    static real wrk[1200];

    /* Fortran I/O blocks */
    static cilist io___29 = { 0, 6, 0, fmt_910, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_915, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_920, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_925, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_930, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_935, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_940, 0 };
    static cilist io___36 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_945, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_955, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_950, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_960, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_965, 0 };


/*  the data parameter values */
/*  the data absciss values */
/*  the data ordinate values */
/*  m denotes the number of data points */
    m = 32;
/*  we set up the weights of the data points */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__ - 1] = 1.f;
/* L10: */
    }
/*  we set up the dimension information. */
    nest = 40;
    lwrk = 1200;
    nc = 80;
    mx = 64;
/*  we will determine a planar curve   x=sx(u) , y=sy(u) */
    idim = 2;
/*  for the first approximations we will use cubic splines */
    k = 3;
/*  we will also supply the parameter values u(i) */
    ipar = 1;
    ub = 120.f;
    ue = 320.f;
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
/*  we start computing a polynomial curve ( s very large) */
L110:
	iopt = 0;
	s = 100.f;
	goto L300;
/*  iopt =  1 from the second call on */
L120:
	iopt = 1;
	s = 1.f;
	goto L300;
/*  a smaller value for s to get a closer approximation */
L130:
	s = .05f;
	goto L300;
/*  a larger value for s to get a smoother approximation */
L140:
	s = .25f;
	goto L300;
/*  if a satisfactory fit is obtained we can calculate a curve of equal */
/*  quality of fit (same value for s) but possibly with fewer knots by */
/*  specifying iopt=0 */
L150:
	iopt = 0;
	s = .25f;
	goto L300;
/*  we determine a spline curve with respect to the same smoothing */
/*  factor s,  but now we let the program determine parameter values u(i) */
L160:
	ipar = 0;
	iopt = 0;
	s = .25f;
	goto L300;
/*  we choose a different degree of spline approximation */
L170:
	k = 5;
	iopt = 0;
	s = .25f;
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
	del = (ue - ub) * .125f;
	for (l = 1; l <= 7; ++l) {
	    al = (real) l;
	    t[j - 1] = ub + al * del;
	    ++j;
/* L200: */
	}
/*  determine the approximating curve */
L300:
	parcur_(&iopt, &ipar, &idim, &m, u, &mx, x, w, &ub, &ue, &k, &s, &nest, &n, t, &nc, c__, &fp, wrk, &lwrk, iwrk, &ier);
/*  printing of the results. */
	if (iopt >= 0) {
	    goto L310;
	}
	s_wsfe(&io___29);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L320;
L310:
	s_wsfe(&io___30);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipar, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&s, (ftnlen)sizeof(real));
	e_wsfe();
L320:
	s_wsfe(&io___32);
	do_fio(&c__1, (char *)&fp, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ier, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___33);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___34);
	e_wsfe();
	if (ipar == 1) {
	    s_wsfe(&io___35);
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	if (ipar == 0) {
	    s_wsfe(&io___36);
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&t[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	nk1 = n - k - 1;
	s_wsfe(&io___38);
	e_wsfe();
	s_wsfe(&io___39);
	i__1 = nk1;
	for (l = 1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&c__[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___40);
	e_wsfe();
	i1 = n + 1;
	i2 = n + nk1;
	s_wsfe(&io___43);
	i__1 = i2;
	for (l = i1; l <= i__1; ++l) {
	    do_fio(&c__1, (char *)&c__[l - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___44);
	e_wsfe();
/*  we evaluate the spline curve */
	curev_(&idim, t, &n, c__, &nc, &k, u, &m, sp, &mx, &ier);
	for (i__ = 1; i__ <= 8; ++i__) {
	    l = (i__ - 1 << 3) + 3;
	    l1 = l + 1;
	    j = l + 4;
	    j1 = j + 1;
	    s_wsfe(&io___48);
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

