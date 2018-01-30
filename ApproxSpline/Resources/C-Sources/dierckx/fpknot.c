/* fpknot.f -- translated by f2c (version 20061008).
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

/*<       subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart) >*/
/* Subroutine */ int fpknot_(doublereal *x, integer *m, doublereal *t, 
	integer *n, doublereal *fpint, integer *nrdata, integer *nrint, 
	integer *nest, integer *istart)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static doublereal am, an;
    static integer jj, jk, nrx, next, ihalf;
    static doublereal fpmax;
    static integer maxpt, jbegin, maxbeg, number, jpoint;

/*  subroutine fpknot locates an additional knot for a spline of degree */
/*  k and adjusts the corresponding parameters,i.e. */
/*    t     : the position of the knots. */
/*    n     : the number of knots. */
/*    nrint : the number of knotintervals. */
/*    fpint : the sum of squares of residual right hand sides */
/*            for each knot interval. */
/*    nrdata: the number of data points inside each knot interval. */
/*  istart indicates that the smallest data point at which the new knot */
/*  may be added is x(istart+1) */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer m,n,nrint,nest,istart >*/
/*  ..array arguments.. */
/*<       real x(m),t(nest),fpint(nest) >*/
/*<       integer nrdata(nest) >*/
/*  ..local scalars.. */
/*<       real an,am,fpmax >*/
/*<    >*/
/*  .. */
/*<       k = (n-nrint-1)/2 >*/
    /* Parameter adjustments */
    --x;
    --nrdata;
    --fpint;
    --t;

    /* Function Body */
    k = (*n - *nrint - 1) / 2;
/*  search for knot interval t(number+k) <= x <= t(number+k+1) where */
/*  fpint(number) is maximal on the condition that nrdata(number) */
/*  not equals zero. */
/*<       fpmax = 0. >*/
    fpmax = 0.;
/*<       jbegin = istart >*/
    jbegin = *istart;
/*<       do 20 j=1,nrint >*/
    i__1 = *nrint;
    for (j = 1; j <= i__1; ++j) {
/*<         jpoint = nrdata(j) >*/
	jpoint = nrdata[j];
/*<         if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10 >*/
	if (fpmax >= fpint[j] || jpoint == 0) {
	    goto L10;
	}
/*<         fpmax = fpint(j) >*/
	fpmax = fpint[j];
/*<         number = j >*/
	number = j;
/*<         maxpt = jpoint >*/
	maxpt = jpoint;
/*<         maxbeg = jbegin >*/
	maxbeg = jbegin;
/*<   10    jbegin = jbegin+jpoint+1 >*/
L10:
	jbegin = jbegin + jpoint + 1;
/*<   20  continue >*/
/* L20: */
    }
/*  let coincide the new knot t(number+k+1) with a data point x(nrx) */
/*  inside the old knot interval t(number+k) <= x <= t(number+k+1). */
/*<       ihalf = maxpt/2+1 >*/
    ihalf = maxpt / 2 + 1;
/*<       nrx = maxbeg+ihalf >*/
    nrx = maxbeg + ihalf;
/*<       next = number+1 >*/
    next = number + 1;
/*<       if(next.gt.nrint) go to 40 >*/
    if (next > *nrint) {
	goto L40;
    }
/*  adjust the different parameters. */
/*<       do 30 j=next,nrint >*/
    i__1 = *nrint;
    for (j = next; j <= i__1; ++j) {
/*<         jj = next+nrint-j >*/
	jj = next + *nrint - j;
/*<         fpint(jj+1) = fpint(jj) >*/
	fpint[jj + 1] = fpint[jj];
/*<         nrdata(jj+1) = nrdata(jj) >*/
	nrdata[jj + 1] = nrdata[jj];
/*<         jk = jj+k >*/
	jk = jj + k;
/*<         t(jk+1) = t(jk) >*/
	t[jk + 1] = t[jk];
/*<   30  continue >*/
/* L30: */
    }
/*<   40  nrdata(number) = ihalf-1 >*/
L40:
    nrdata[number] = ihalf - 1;
/*<       nrdata(next) = maxpt-ihalf >*/
    nrdata[next] = maxpt - ihalf;
/*<       am = maxpt >*/
    am = (doublereal) maxpt;
/*<       an = nrdata(number) >*/
    an = (doublereal) nrdata[number];
/*<       fpint(number) = fpmax*an/am >*/
    fpint[number] = fpmax * an / am;
/*<       an = nrdata(next) >*/
    an = (doublereal) nrdata[next];
/*<       fpint(next) = fpmax*an/am >*/
    fpint[next] = fpmax * an / am;
/*<       jk = next+k >*/
    jk = next + k;
/*<       t(jk) = x(nrx) >*/
    t[jk] = x[nrx];
/*<       n = n+1 >*/
    ++(*n);
/*<       nrint = nrint+1 >*/
    ++(*nrint);
/*<       return >*/
    return 0;
/*<       end >*/
} /* fpknot_ */

