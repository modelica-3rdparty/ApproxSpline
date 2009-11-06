/* fppocu.f -- translated by f2c (version 20061008).
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

/*<       subroutine fppocu(idim,k,a,b,ib,db,nb,ie,de,ne,cp,np) >*/
/* Subroutine */ int fppocu_(integer *idim, integer *k, doublereal *a, 
	doublereal *b, integer *ib, doublereal *db, integer *nb, integer *ie, 
	doublereal *de, integer *ne, doublereal *cp, integer *np)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, k1, k2;
    static doublereal ab;
    static integer id, jj, ll;
    static doublereal aki, work[36]	/* was [6][6] */;

/*  subroutine fppocu finds a idim-dimensional polynomial curve p(u) = */
/*  (p1(u),p2(u),...,pidim(u)) of degree k, satisfying certain derivative */
/*  constraints at the end points a and b, i.e. */
/*                  (l) */
/*    if ib > 0 : pj   (a) = db(idim*l+j), l=0,1,...,ib-1 */
/*                  (l) */
/*    if ie > 0 : pj   (b) = de(idim*l+j), l=0,1,...,ie-1 */

/*  the polynomial curve is returned in its b-spline representation */
/*  ( coefficients cp(j), j=1,2,...,np ) */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer idim,k,ib,nb,ie,ne,np >*/
/*<       real a,b >*/
/*  ..array arguments.. */
/*<       real db(nb),de(ne),cp(np) >*/
/*  ..local scalars.. */
/*<       real ab,aki >*/
/*<       integer i,id,j,jj,l,ll,k1,k2 >*/
/*  ..local array.. */
/*<       real work(6,6) >*/
/*  .. */
/*<       k1 = k+1 >*/
    /* Parameter adjustments */
    --db;
    --de;
    --cp;

    /* Function Body */
    k1 = *k + 1;
/*<       k2 = 2*k1 >*/
    k2 = k1 << 1;
/*<       ab = b-a >*/
    ab = *b - *a;
/*<       do 110 id=1,idim >*/
    i__1 = *idim;
    for (id = 1; id <= i__1; ++id) {
/*<         do 10 j=1,k1 >*/
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
/*<           work(j,1) = 0. >*/
	    work[j - 1] = 0.;
/*<   10    continue >*/
/* L10: */
	}
/*<         if(ib.eq.0) go to 50 >*/
	if (*ib == 0) {
	    goto L50;
	}
/*<         l = id >*/
	l = id;
/*<         do 20 i=1,ib >*/
	i__2 = *ib;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           work(1,i) = db(l) >*/
	    work[i__ * 6 - 6] = db[l];
/*<           l = l+idim >*/
	    l += *idim;
/*<   20    continue >*/
/* L20: */
	}
/*<         if(ib.eq.1) go to 50 >*/
	if (*ib == 1) {
	    goto L50;
	}
/*<         ll = ib >*/
	ll = *ib;
/*<         do 40 j=2,ib >*/
	i__2 = *ib;
	for (j = 2; j <= i__2; ++j) {
/*<           ll =  ll-1 >*/
	    --ll;
/*<           do 30 i=1,ll >*/
	    i__3 = ll;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             aki = k1-i >*/
		aki = (doublereal) (k1 - i__);
/*<             work(j,i) = ab*work(j-1,i+1)/aki + work(j-1,i) >*/
		work[j + i__ * 6 - 7] = ab * work[j - 1 + (i__ + 1) * 6 - 7] /
			 aki + work[j - 1 + i__ * 6 - 7];
/*<   30      continue >*/
/* L30: */
	    }
/*<   40    continue >*/
/* L40: */
	}
/*<   50    if(ie.eq.0) go to 90 >*/
L50:
	if (*ie == 0) {
	    goto L90;
	}
/*<         l = id >*/
	l = id;
/*<         j = k1 >*/
	j = k1;
/*<         do 60 i=1,ie >*/
	i__2 = *ie;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<           work(j,i) = de(l) >*/
	    work[j + i__ * 6 - 7] = de[l];
/*<           l = l+idim >*/
	    l += *idim;
/*<           j = j-1 >*/
	    --j;
/*<   60    continue >*/
/* L60: */
	}
/*<         if(ie.eq.1) go to 90 >*/
	if (*ie == 1) {
	    goto L90;
	}
/*<         ll = ie >*/
	ll = *ie;
/*<         do 80 jj=2,ie >*/
	i__2 = *ie;
	for (jj = 2; jj <= i__2; ++jj) {
/*<           ll =  ll-1 >*/
	    --ll;
/*<           j = k1+1-jj >*/
	    j = k1 + 1 - jj;
/*<           do 70 i=1,ll >*/
	    i__3 = ll;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*<             aki = k1-i >*/
		aki = (doublereal) (k1 - i__);
/*<             work(j,i) = work(j+1,i) - ab*work(j,i+1)/aki >*/
		work[j + i__ * 6 - 7] = work[j + 1 + i__ * 6 - 7] - ab * work[
			j + (i__ + 1) * 6 - 7] / aki;
/*<             j = j-1 >*/
		--j;
/*<   70      continue >*/
/* L70: */
	    }
/*<   80    continue >*/
/* L80: */
	}
/*<   90    l = (id-1)*k2 >*/
L90:
	l = (id - 1) * k2;
/*<         do 100 j=1,k1 >*/
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
/*<           l = l+1 >*/
	    ++l;
/*<           cp(l) = work(j,1) >*/
	    cp[l] = work[j - 1];
/*<  100    continue >*/
/* L100: */
	}
/*<  110  continue >*/
/* L110: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fppocu_ */

