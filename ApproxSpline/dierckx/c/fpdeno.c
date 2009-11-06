/* fpdeno.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpdeno(maxtr,up,left,right,nbind,merk) >*/
/* Subroutine */ int fpdeno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *nbind, integer *merk)
{
    static integer i__, j, k, l, point, niveau;

/*  subroutine fpdeno frees the nodes of all branches of a triply linked */
/*  tree with length < nbind by putting to zero their up field. */
/*  on exit the parameter merk points to the terminal node of the */
/*  most left branch of length nbind or takes the value 1 if there */
/*  is no such branch. */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer maxtr,nbind,merk >*/
/*  ..array arguments.. */
/*<       integer up(maxtr),left(maxtr),right(maxtr) >*/
/*  ..local scalars .. */
/*<       integer i,j,k,l,niveau,point >*/
/*  .. */
/*<       i = 1 >*/
    /* Parameter adjustments */
    --right;
    --left;
    --up;

    /* Function Body */
    i__ = 1;
/*<       niveau = 0 >*/
    niveau = 0;
/*<   10  point = i >*/
L10:
    point = i__;
/*<       i = left(point) >*/
    i__ = left[point];
/*<       if(i.eq.0) go to 20 >*/
    if (i__ == 0) {
	goto L20;
    }
/*<       niveau = niveau+1 >*/
    ++niveau;
/*<       go to 10 >*/
    goto L10;
/*<   20  if(niveau.eq.nbind) go to 70 >*/
L20:
    if (niveau == *nbind) {
	goto L70;
    }
/*<   30  i = right(point) >*/
L30:
    i__ = right[point];
/*<       j = up(point) >*/
    j = up[point];
/*<       up(point) = 0 >*/
    up[point] = 0;
/*<       k = left(j) >*/
    k = left[j];
/*<       if(point.ne.k) go to 50 >*/
    if (point != k) {
	goto L50;
    }
/*<       if(i.ne.0) go to 40 >*/
    if (i__ != 0) {
	goto L40;
    }
/*<       niveau = niveau-1 >*/
    --niveau;
/*<       if(niveau.eq.0) go to 80 >*/
    if (niveau == 0) {
	goto L80;
    }
/*<       point = j >*/
    point = j;
/*<       go to 30 >*/
    goto L30;
/*<   40  left(j) = i >*/
L40:
    left[j] = i__;
/*<       go to 10 >*/
    goto L10;
/*<   50  l = right(k) >*/
L50:
    l = right[k];
/*<       if(point.eq.l) go to 60 >*/
    if (point == l) {
	goto L60;
    }
/*<       k = l >*/
    k = l;
/*<       go to 50 >*/
    goto L50;
/*<   60  right(k) = i >*/
L60:
    right[k] = i__;
/*<       point = k >*/
    point = k;
/*<   70  i = right(point) >*/
L70:
    i__ = right[point];
/*<       if(i.ne.0) go to 10 >*/
    if (i__ != 0) {
	goto L10;
    }
/*<       i = up(point) >*/
    i__ = up[point];
/*<       niveau = niveau-1 >*/
    --niveau;
/*<       if(niveau.eq.0) go to 80 >*/
    if (niveau == 0) {
	goto L80;
    }
/*<       point = i >*/
    point = i__;
/*<       go to 70 >*/
    goto L70;
/*<   80  k = 1 >*/
L80:
    k = 1;
/*<       l = left(k) >*/
    l = left[k];
/*<       if(up(l).eq.0) return >*/
    if (up[l] == 0) {
	return 0;
    }
/*<   90  merk = k >*/
L90:
    *merk = k;
/*<       k = left(k) >*/
    k = left[k];
/*<       if(k.ne.0) go to 90 >*/
    if (k != 0) {
	goto L90;
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpdeno_ */

