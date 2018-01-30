/* fpcyt1.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpcyt1(a,n,nn) >*/
/* Subroutine */ int fpcyt1_(doublereal *a, integer *n, integer *nn)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__;
    static doublereal v;
    static integer n1, n2;
    static doublereal aa, one, sum, beta, teta, gamma;

/* (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero */
/* elements stored as follows */

/*    | a(1,2) a(1,3)                                    a(1,1)  | */
/*    | a(2,1) a(2,2) a(2,3)                                     | */
/*    |        a(3,1) a(3,2) a(3,3)                              | */
/*    |               ...............                            | */
/*    |                               a(n-1,1) a(n-1,2) a(n-1,3) | */
/*    | a(n,3)                                  a(n,1)   a(n,2)  | */

/*  .. */
/*  ..scalar arguments.. */
/*<       integer n,nn >*/
/*  ..array arguments.. */
/*<       real a(nn,6) >*/
/*  ..local scalars.. */
/*<       real aa,beta,gamma,sum,teta,v,one >*/
/*<       integer i,n1,n2 >*/
/*  .. */
/*  set constant */
/*<       one = 1 >*/
    /* Parameter adjustments */
    a_dim1 = *nn;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    one = 1.;
/*<       n2 = n-2 >*/
    n2 = *n - 2;
/*<       beta = one/a(1,2) >*/
    beta = one / a[(a_dim1 << 1) + 1];
/*<       gamma = a(n,3) >*/
    gamma = a[*n + a_dim1 * 3];
/*<       teta = a(1,1)*beta >*/
    teta = a[a_dim1 + 1] * beta;
/*<       a(1,4) = beta >*/
    a[(a_dim1 << 2) + 1] = beta;
/*<       a(1,5) = gamma >*/
    a[a_dim1 * 5 + 1] = gamma;
/*<       a(1,6) = teta >*/
    a[a_dim1 * 6 + 1] = teta;
/*<       sum = gamma*teta >*/
    sum = gamma * teta;
/*<       do 10 i=2,n2 >*/
    i__1 = n2;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          v = a(i-1,3)*beta >*/
	v = a[i__ - 1 + a_dim1 * 3] * beta;
/*<          aa = a(i,1) >*/
	aa = a[i__ + a_dim1];
/*<          beta = one/(a(i,2)-aa*v) >*/
	beta = one / (a[i__ + (a_dim1 << 1)] - aa * v);
/*<          gamma = -gamma*v >*/
	gamma = -gamma * v;
/*<          teta = -teta*aa*beta >*/
	teta = -teta * aa * beta;
/*<          a(i,4) = beta >*/
	a[i__ + (a_dim1 << 2)] = beta;
/*<          a(i,5) = gamma >*/
	a[i__ + a_dim1 * 5] = gamma;
/*<          a(i,6) = teta >*/
	a[i__ + a_dim1 * 6] = teta;
/*<          sum = sum+gamma*teta >*/
	sum += gamma * teta;
/*<   10  continue >*/
/* L10: */
    }
/*<       n1 = n-1 >*/
    n1 = *n - 1;
/*<       v = a(n2,3)*beta >*/
    v = a[n2 + a_dim1 * 3] * beta;
/*<       aa = a(n1,1) >*/
    aa = a[n1 + a_dim1];
/*<       beta = one/(a(n1,2)-aa*v) >*/
    beta = one / (a[n1 + (a_dim1 << 1)] - aa * v);
/*<       gamma = a(n,1)-gamma*v >*/
    gamma = a[*n + a_dim1] - gamma * v;
/*<       teta = (a(n1,3)-teta*aa)*beta >*/
    teta = (a[n1 + a_dim1 * 3] - teta * aa) * beta;
/*<       a(n1,4) = beta >*/
    a[n1 + (a_dim1 << 2)] = beta;
/*<       a(n1,5) = gamma >*/
    a[n1 + a_dim1 * 5] = gamma;
/*<       a(n1,6) = teta >*/
    a[n1 + a_dim1 * 6] = teta;
/*<       a(n,4) = one/(a(n,2)-(sum+gamma*teta)) >*/
    a[*n + (a_dim1 << 2)] = one / (a[*n + (a_dim1 << 1)] - (sum + gamma * 
	    teta));
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpcyt1_ */

