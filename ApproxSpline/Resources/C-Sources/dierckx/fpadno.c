/* fpadno.f -- translated by f2c (version 20061008).
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
/* Subroutine */ int fpadno_(integer *maxtr, integer *up, integer *left, 
	integer *right, integer *info, integer *count, integer *merk, integer 
	*jbind, integer *n1, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static logical bool;
    static integer point, niveau;
    extern /* Subroutine */ int fpfrno_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);

/*  subroutine fpadno adds a branch of length n1 to the triply linked */
/*  tree,the information of which is kept in the arrays up,left,right */
/*  and info. the information field of the nodes of this new branch is */
/*  given in the array jbind. in linking the new branch fpadno takes */
/*  account of the property of the tree that */
/*    info(k) < info(right(k)) ; info(k) < info(left(k)) */
/*  if necessary the subroutine calls subroutine fpfrno to collect the */
/*  free nodes of the tree. if no computer words are available at that */
/*  moment, the error parameter ier is set to 1. */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer maxtr,count,merk,n1,ier >*/
/*  ..array arguments.. */
/*<       integer up(maxtr),left(maxtr),right(maxtr),info(maxtr),jbind(n1) >*/
/*  ..local scalars.. */
/*<       integer k,niveau,point >*/
/*<       logical bool >*/
/*  ..subroutine references.. */
/*    fpfrno */
/*  .. */
/*<       point = 1 >*/
    /* Parameter adjustments */
    --info;
    --right;
    --left;
    --up;
    --jbind;

    /* Function Body */
    point = 1;
/*<       niveau = 1 >*/
    niveau = 1;
/*<   10  k = left(point) >*/
L10:
    k = left[point];
/*<       bool = .true. >*/
    bool = TRUE_;
/*<   20  if(k.eq.0) go to 50 >*/
L20:
    if (k == 0) {
	goto L50;
    }
/*<       if(info(k)-jbind(niveau)) 30,40,50 >*/
    if ((i__1 = info[k] - jbind[niveau]) < 0) {
	goto L30;
    } else if (i__1 == 0) {
	goto L40;
    } else {
	goto L50;
    }
/*<   30  point = k >*/
L30:
    point = k;
/*<       k = right(point) >*/
    k = right[point];
/*<       bool = .false. >*/
    bool = FALSE_;
/*<       go to 20 >*/
    goto L20;
/*<   40  point = k >*/
L40:
    point = k;
/*<       niveau = niveau+1 >*/
    ++niveau;
/*<       go to 10 >*/
    goto L10;
/*<   50  if(niveau.gt.n1) go to 90 >*/
L50:
    if (niveau > *n1) {
	goto L90;
    }
/*<       count = count+1 >*/
    ++(*count);
/*<       if(count.le.maxtr) go to 60 >*/
    if (*count <= *maxtr) {
	goto L60;
    }
/*<       call fpfrno(maxtr,up,left,right,info,point,merk,n1,count,ier) >*/
    fpfrno_(maxtr, &up[1], &left[1], &right[1], &info[1], &point, merk, n1, 
	    count, ier);
/*<       if(ier.ne.0) go to 100 >*/
    if (*ier != 0) {
	goto L100;
    }
/*<   60  info(count) = jbind(niveau) >*/
L60:
    info[*count] = jbind[niveau];
/*<       left(count) = 0 >*/
    left[*count] = 0;
/*<       right(count) = k >*/
    right[*count] = k;
/*<       if(bool) go to 70 >*/
    if (bool) {
	goto L70;
    }
/*<       bool = .true. >*/
    bool = TRUE_;
/*<       right(point) = count >*/
    right[point] = *count;
/*<       up(count) = up(point) >*/
    up[*count] = up[point];
/*<       go to 80 >*/
    goto L80;
/*<   70  up(count) = point >*/
L70:
    up[*count] = point;
/*<       left(point) = count >*/
    left[point] = *count;
/*<   80  point = count >*/
L80:
    point = *count;
/*<       niveau = niveau+1 >*/
    ++niveau;
/*<       k = 0 >*/
    k = 0;
/*<       go to 50 >*/
    goto L50;
/*<   90  ier = 0 >*/
L90:
    *ier = 0;
/*<  100  return >*/
L100:
    return 0;
/*<       end >*/
} /* fpadno_ */

