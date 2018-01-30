/* fpcsin.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpcsin(a,b,par,sia,coa,sib,cob,ress,resc) >*/
/* Subroutine */ int fpcsin_(doublereal *a, doublereal *b, doublereal *par, 
	doublereal *sia, doublereal *coa, doublereal *sib, doublereal *cob, 
	doublereal *ress, doublereal *resc)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal b2, b4, f1, f2, ab, ai, ab4, fac, one, eps, six, two, 
	    alfa, beta, three, quart;

/*  fpcsin calculates the integrals ress=integral((b-x)**3*sin(par*x)) */
/*  and resc=integral((b-x)**3*cos(par*x)) over the interval (a,b), */
/*  given sia=sin(par*a),coa=cos(par*a),sib=sin(par*b) and cob=cos(par*b) */
/*  .. */
/*  ..scalar arguments.. */
/*<       real a,b,par,sia,coa,sib,cob,ress,resc >*/
/*  ..local scalars.. */
/*<       integer i,j >*/
/*<    >*/
/*  ..function references.. */
/*<       real abs >*/
/*  .. */
/*<       one = 0.1e+01 >*/
    one = 1.;
/*<       two = 0.2e+01 >*/
    two = 2.;
/*<       three = 0.3e+01 >*/
    three = 3.;
/*<       six = 0.6e+01 >*/
    six = 6.;
/*<       quart = 0.25e+0 >*/
    quart = .25;
/*<       eps = 0.1e-09 >*/
    eps = 1e-10;
/*<       ab = b-a >*/
    ab = *b - *a;
/*<       ab4 = ab**4 >*/
/* Computing 4th power */
    d__1 = ab, d__1 *= d__1;
    ab4 = d__1 * d__1;
/*<       alfa = ab*par >*/
    alfa = ab * *par;
/* the way of calculating the integrals ress and resc depends on */
/* the value of alfa = (b-a)*par. */
/*<       if(abs(alfa).le.one) go to 100 >*/
    if (abs(alfa) <= one) {
	goto L100;
    }
/* integration by parts. */
/*<       beta = one/alfa >*/
    beta = one / alfa;
/*<       b2 = beta**2 >*/
/* Computing 2nd power */
    d__1 = beta;
    b2 = d__1 * d__1;
/*<       b4 = six*b2**2 >*/
/* Computing 2nd power */
    d__1 = b2;
    b4 = six * (d__1 * d__1);
/*<       f1 = three*b2*(one-two*b2) >*/
    f1 = three * b2 * (one - two * b2);
/*<       f2 = beta*(one-six*b2) >*/
    f2 = beta * (one - six * b2);
/*<       ress = ab4*(coa*f2+sia*f1+sib*b4) >*/
    *ress = ab4 * (*coa * f2 + *sia * f1 + *sib * b4);
/*<       resc = ab4*(coa*f1-sia*f2+cob*b4) >*/
    *resc = ab4 * (*coa * f1 - *sia * f2 + *cob * b4);
/*<       go to 400 >*/
    goto L400;
/* ress and resc are found by evaluating a series expansion. */
/*<  100  fac = quart >*/
L100:
    fac = quart;
/*<       f1 = fac >*/
    f1 = fac;
/*<       f2 = 0. >*/
    f2 = 0.;
/*<       i = 4 >*/
    i__ = 4;
/*<       do 200 j=1,5 >*/
    for (j = 1; j <= 5; ++j) {
/*<         i = i+1 >*/
	++i__;
/*<         ai = i >*/
	ai = (doublereal) i__;
/*<         fac = fac*alfa/ai >*/
	fac = fac * alfa / ai;
/*<         f2 = f2+fac >*/
	f2 += fac;
/*<         if(abs(fac).le.eps) go to 300 >*/
	if (abs(fac) <= eps) {
	    goto L300;
	}
/*<         i = i+1 >*/
	++i__;
/*<         ai = i >*/
	ai = (doublereal) i__;
/*<         fac = -fac*alfa/ai >*/
	fac = -fac * alfa / ai;
/*<         f1 = f1+fac >*/
	f1 += fac;
/*<         if(abs(fac).le.eps) go to 300 >*/
	if (abs(fac) <= eps) {
	    goto L300;
	}
/*<  200  continue >*/
/* L200: */
    }
/*<  300  ress = ab4*(coa*f2+sia*f1) >*/
L300:
    *ress = ab4 * (*coa * f2 + *sia * f1);
/*<       resc = ab4*(coa*f1-sia*f2) >*/
    *resc = ab4 * (*coa * f1 - *sia * f2);
/*<  400  return >*/
L400:
    return 0;
/*<       end >*/
} /* fpcsin_ */

