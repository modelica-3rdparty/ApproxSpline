/* fpfrno.f -- translated by f2c (version 20061008).
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

/*<    >*/
/* Subroutine */ int fpfrno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *info, integer *point, integer *merk, integer 
	*n1, integer *count, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, n, niveau;

/*  subroutine fpfrno collects the free nodes (up field zero) of the */
/*  triply linked tree the information of which is kept in the arrays */
/*  up,left,right and info. the maximal length of the branches of the */
/*  tree is given by n1. if no free nodes are found, the error flag */
/*  ier is set to 1. */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer maxtr,point,merk,n1,count,ier >*/
/*  ..array arguments.. */
/*<       integer up(maxtr),left(maxtr),right(maxtr),info(maxtr) >*/
/*  ..local scalars */
/*<       integer i,j,k,l,n,niveau >*/
/*  .. */
/*<       ier = 1 >*/
    /* Parameter adjustments */
    --info;
    --right;
    --left;
    --up;

    /* Function Body */
    *ier = 1;
/*<       if(n1.eq.2) go to 140 >*/
    if (*n1 == 2) {
	goto L140;
    }
/*<       niveau = 1 >*/
    niveau = 1;
/*<       count = 2 >*/
    *count = 2;
/*<   10  j = 0 >*/
L10:
    j = 0;
/*<       i = 1 >*/
    i__ = 1;
/*<   20  if(j.eq.niveau) go to 30 >*/
L20:
    if (j == niveau) {
	goto L30;
    }
/*<       k = 0 >*/
    k = 0;
/*<       l = left(i) >*/
    l = left[i__];
/*<       if(l.eq.0) go to 110 >*/
    if (l == 0) {
	goto L110;
    }
/*<       i = l >*/
    i__ = l;
/*<       j = j+1 >*/
    ++j;
/*<       go to 20 >*/
    goto L20;
/*<   30  if(i-count) 110,100,40 >*/
L30:
    if ((i__1 = i__ - *count) < 0) {
	goto L110;
    } else if (i__1 == 0) {
	goto L100;
    } else {
	goto L40;
    }
/*<   40  if(up(count).eq.0) go to 50 >*/
L40:
    if (up[*count] == 0) {
	goto L50;
    }
/*<       count = count+1 >*/
    ++(*count);
/*<       go to 30 >*/
    goto L30;
/*<   50  up(count) = up(i) >*/
L50:
    up[*count] = up[i__];
/*<       left(count) = left(i) >*/
    left[*count] = left[i__];
/*<       right(count) = right(i) >*/
    right[*count] = right[i__];
/*<       info(count) = info(i) >*/
    info[*count] = info[i__];
/*<       if(merk.eq.i) merk = count >*/
    if (*merk == i__) {
	*merk = *count;
    }
/*<       if(point.eq.i) point = count >*/
    if (*point == i__) {
	*point = *count;
    }
/*<       if(k.eq.0) go to 60 >*/
    if (k == 0) {
	goto L60;
    }
/*<       right(k) = count >*/
    right[k] = *count;
/*<       go to 70 >*/
    goto L70;
/*<   60  n = up(i) >*/
L60:
    n = up[i__];
/*<       left(n) = count >*/
    left[n] = *count;
/*<   70  l = left(i) >*/
L70:
    l = left[i__];
/*<   80  if(l.eq.0) go to 90 >*/
L80:
    if (l == 0) {
	goto L90;
    }
/*<       up(l) = count >*/
    up[l] = *count;
/*<       l = right(l) >*/
    l = right[l];
/*<       go to 80 >*/
    goto L80;
/*<   90  up(i) = 0 >*/
L90:
    up[i__] = 0;
/*<       i = count >*/
    i__ = *count;
/*<  100  count = count+1 >*/
L100:
    ++(*count);
/*<  110  l = right(i) >*/
L110:
    l = right[i__];
/*<       k = i >*/
    k = i__;
/*<       if(l.eq.0) go to 120 >*/
    if (l == 0) {
	goto L120;
    }
/*<       i = l >*/
    i__ = l;
/*<       go to 20 >*/
    goto L20;
/*<  120  l = up(i) >*/
L120:
    l = up[i__];
/*<       j = j-1 >*/
    --j;
/*<       if(j.eq.0) go to 130 >*/
    if (j == 0) {
	goto L130;
    }
/*<       i = l >*/
    i__ = l;
/*<       go to 110 >*/
    goto L110;
/*<  130  niveau = niveau+1 >*/
L130:
    ++niveau;
/*<       if(niveau.le.n1) go to 10 >*/
    if (niveau <= *n1) {
	goto L10;
    }
/*<       if(count.gt.maxtr) go to 140 >*/
    if (*count > *maxtr) {
	goto L140;
    }
/*<       ier = 0 >*/
    *ier = 0;
/*<  140  return >*/
L140:
    return 0;
/*<       end >*/
} /* fpfrno_ */

