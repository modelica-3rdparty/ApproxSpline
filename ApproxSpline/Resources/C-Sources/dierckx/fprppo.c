/* fprppo.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/*<       subroutine fprppo(nu,nv,if1,if2,cosi,ratio,c,f,ncoff) >*/
/* Subroutine */ int fprppo_(integer *nu, integer *nv, integer *if1, integer *
	if2, doublereal *cosi, doublereal *ratio, doublereal *c__, doublereal 
	*f, integer *ncoff)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, ii, nu4, nvv, iopt;

/*  given the coefficients of a constrained bicubic spline, as determined */
/*  in subroutine fppola, subroutine fprppo calculates the coefficients */
/*  in the standard b-spline representation of bicubic splines. */
/*  .. */
/*  ..scalar arguments.. */
/*<       real ratio >*/
/*<       integer nu,nv,if1,if2,ncoff >*/
/*  ..array arguments */
/*<       real c(ncoff),f(ncoff),cosi(5,nv) >*/
/*  ..local scalars.. */
/*<       integer i,iopt,ii,j,k,l,nu4,nvv >*/
/*  .. */
/*<       nu4 = nu-4 >*/
    /* Parameter adjustments */
    cosi -= 6;
    --f;
    --c__;

    /* Function Body */
    nu4 = *nu - 4;
/*<       nvv = nv-7 >*/
    nvv = *nv - 7;
/*<       iopt = if1+1 >*/
    iopt = *if1 + 1;
/*<       do 10 i=1,ncoff >*/
    i__1 = *ncoff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          f(i) = 0. >*/
	f[i__] = 0.;
/*<   10  continue >*/
/* L10: */
    }
/*<       i = 0 >*/
    i__ = 0;
/*<       do 120 l=1,nu4 >*/
    i__1 = nu4;
    for (l = 1; l <= i__1; ++l) {
/*<          ii = i >*/
	ii = i__;
/*<          if(l.gt.iopt) go to 80 >*/
	if (l > iopt) {
	    goto L80;
	}
/*<          go to (20,40,60),l >*/
	switch (l) {
	    case 1:  goto L20;
	    case 2:  goto L40;
	    case 3:  goto L60;
	}
/*<   20     do 30 k=1,nvv >*/
L20:
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             f(i) = c(1) >*/
	    f[i__] = c__[1];
/*<   30     continue >*/
/* L30: */
	}
/*<          j = 1 >*/
	j = 1;
/*<          go to 100 >*/
	goto L100;
/*<   40     do 50 k=1,nvv >*/
L40:
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             f(i) = c(1)+c(2)*cosi(1,k)+c(3)*cosi(2,k) >*/
	    f[i__] = c__[1] + c__[2] * cosi[k * 5 + 1] + c__[3] * cosi[k * 5 
		    + 2];
/*<   50     continue >*/
/* L50: */
	}
/*<          j = 3 >*/
	j = 3;
/*<          go to 100 >*/
	goto L100;
/*<   60     do 70 k=1,nvv >*/
L60:
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<    >*/
	    f[i__] = c__[1] + *ratio * (c__[2] * cosi[k * 5 + 1] + c__[3] * 
		    cosi[k * 5 + 2]) + c__[4] * cosi[k * 5 + 3] + c__[5] * 
		    cosi[k * 5 + 4] + c__[6] * cosi[k * 5 + 5];
/*<   70     continue >*/
/* L70: */
	}
/*<          j = 6 >*/
	j = 6;
/*<          go to 100 >*/
	goto L100;
/*<   80     if(l.eq.nu4 .and. if2.ne.0) go to 120 >*/
L80:
	if (l == nu4 && *if2 != 0) {
	    goto L120;
	}
/*<          do 90 k=1,nvv >*/
	i__2 = nvv;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             j = j+1 >*/
	    ++j;
/*<             f(i) = c(j) >*/
	    f[i__] = c__[j];
/*<   90     continue >*/
/* L90: */
	}
/*<  100     do 110 k=1,3 >*/
L100:
	for (k = 1; k <= 3; ++k) {
/*<             ii = ii+1 >*/
	    ++ii;
/*<             i = i+1 >*/
	    ++i__;
/*<             f(i) = f(ii) >*/
	    f[i__] = f[ii];
/*<  110     continue >*/
/* L110: */
	}
/*<  120  continue >*/
L120:
	;
    }
/*<       do 130 i=1,ncoff >*/
    i__1 = *ncoff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          c(i) = f(i) >*/
	c__[i__] = f[i__];
/*<  130  continue >*/
/* L130: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fprppo_ */

