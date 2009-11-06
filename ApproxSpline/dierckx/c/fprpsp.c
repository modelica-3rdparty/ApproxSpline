/* fprpsp.f -- translated by f2c (version 20061008).
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

/*<       subroutine fprpsp(nt,np,co,si,c,f,ncoff) >*/
/* Subroutine */ int fprpsp_(integer *nt, integer *np, doublereal *co, 
	doublereal *si, doublereal *c__, doublereal *f, integer *ncoff)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal c1, c2, c3, cn;
    static integer ii, np4, nt4, npp, ncof;

/*  given the coefficients of a spherical spline function, subroutine */
/*  fprpsp calculates the coefficients in the standard b-spline re- */
/*  presentation of this bicubic spline. */
/*  .. */
/*  ..scalar arguments */
/*<       integer nt,np,ncoff >*/
/*  ..array arguments */
/*<       real co(np),si(np),c(ncoff),f(ncoff) >*/
/*  ..local scalars */
/*<       real cn,c1,c2,c3 >*/
/*<       integer i,ii,j,k,l,ncof,npp,np4,nt4 >*/
/*  .. */
/*<       nt4 = nt-4 >*/
    /* Parameter adjustments */
    --si;
    --co;
    --f;
    --c__;

    /* Function Body */
    nt4 = *nt - 4;
/*<       np4 = np-4 >*/
    np4 = *np - 4;
/*<       npp = np4-3 >*/
    npp = np4 - 3;
/*<       ncof = 6+npp*(nt4-4) >*/
    ncof = npp * (nt4 - 4) + 6;
/*<       c1 = c(1) >*/
    c1 = c__[1];
/*<       cn = c(ncof) >*/
    cn = c__[ncof];
/*<       j = ncoff >*/
    j = *ncoff;
/*<       do 10 i=1,np4 >*/
    i__1 = np4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          f(i) = c1 >*/
	f[i__] = c1;
/*<          f(j) = cn >*/
	f[j] = cn;
/*<          j = j-1 >*/
	--j;
/*<   10  continue >*/
/* L10: */
    }
/*<       i = np4 >*/
    i__ = np4;
/*<       j=1 >*/
    j = 1;
/*<       do 70 l=3,nt4 >*/
    i__1 = nt4;
    for (l = 3; l <= i__1; ++l) {
/*<          ii = i >*/
	ii = i__;
/*<          if(l.eq.3 .or. l.eq.nt4) go to 30 >*/
	if (l == 3 || l == nt4) {
	    goto L30;
	}
/*<          do 20 k=1,npp >*/
	i__2 = npp;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             j = j+1 >*/
	    ++j;
/*<             f(i) = c(j) >*/
	    f[i__] = c__[j];
/*<   20     continue >*/
/* L20: */
	}
/*<          go to 50 >*/
	goto L50;
/*<   30     if(l.eq.nt4) c1 = cn >*/
L30:
	if (l == nt4) {
	    c1 = cn;
	}
/*<          c2 = c(j+1) >*/
	c2 = c__[j + 1];
/*<          c3 = c(j+2) >*/
	c3 = c__[j + 2];
/*<          j = j+2 >*/
	j += 2;
/*<          do 40 k=1,npp >*/
	i__2 = npp;
	for (k = 1; k <= i__2; ++k) {
/*<             i = i+1 >*/
	    ++i__;
/*<             f(i) = c1+c2*co(k)+c3*si(k) >*/
	    f[i__] = c1 + c2 * co[k] + c3 * si[k];
/*<   40     continue >*/
/* L40: */
	}
/*<   50     do 60 k=1,3 >*/
L50:
	for (k = 1; k <= 3; ++k) {
/*<             ii = ii+1 >*/
	    ++ii;
/*<             i = i+1 >*/
	    ++i__;
/*<             f(i) = f(ii) >*/
	    f[i__] = f[ii];
/*<   60     continue >*/
/* L60: */
	}
/*<   70  continue >*/
/* L70: */
    }
/*<       do 80 i=1,ncoff >*/
    i__1 = *ncoff;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          c(i) = f(i) >*/
	c__[i__] = f[i__];
/*<   80  continue >*/
/* L80: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fprpsp_ */

