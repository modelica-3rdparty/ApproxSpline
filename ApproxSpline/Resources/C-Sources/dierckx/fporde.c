/* fporde.f -- translated by f2c (version 20061008).
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

/*<       subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg) >*/
/* Subroutine */ int fporde_(doublereal *x, doublereal *y, integer *m, 
	integer *kx, integer *ky, doublereal *tx, integer *nx, doublereal *ty,
	 integer *ny, integer *nummer, integer *index, integer *nreg)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, l, k1, l1, im;
    static doublereal xi, yi;
    static integer kx1, ky1, num, nyy, nk1x, nk1y;

/*  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m */
/*  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong */
/*  to. for each panel a stack is constructed  containing the numbers */
/*  of data points lying inside; index(j),j=1,2,...,nreg points to the */
/*  first data point in the jth panel while nummer(i),i=1,2,...,m gives */
/*  the number of the next data point in the panel. */
/*  .. */
/*  ..scalar arguments.. */
/*<       integer m,kx,ky,nx,ny,nreg >*/
/*  ..array arguments.. */
/*<       real x(m),y(m),tx(nx),ty(ny) >*/
/*<       integer nummer(m),index(nreg) >*/
/*  ..local scalars.. */
/*<       real xi,yi >*/
/*<       integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy >*/
/*  .. */
/*<       kx1 = kx+1 >*/
    /* Parameter adjustments */
    --nummer;
    --y;
    --x;
    --tx;
    --ty;
    --index;

    /* Function Body */
    kx1 = *kx + 1;
/*<       ky1 = ky+1 >*/
    ky1 = *ky + 1;
/*<       nk1x = nx-kx1 >*/
    nk1x = *nx - kx1;
/*<       nk1y = ny-ky1 >*/
    nk1y = *ny - ky1;
/*<       nyy = nk1y-ky >*/
    nyy = nk1y - *ky;
/*<       do 10 i=1,nreg >*/
    i__1 = *nreg;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         index(i) = 0 >*/
	index[i__] = 0;
/*<   10  continue >*/
/* L10: */
    }
/*<       do 60 im=1,m >*/
    i__1 = *m;
    for (im = 1; im <= i__1; ++im) {
/*<         xi = x(im) >*/
	xi = x[im];
/*<         yi = y(im) >*/
	yi = y[im];
/*<         l = kx1 >*/
	l = kx1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<   20    if(xi.lt.tx(l1) .or. l.eq.nk1x) go to 30 >*/
L20:
	if (xi < tx[l1] || l == nk1x) {
	    goto L30;
	}
/*<         l = l1 >*/
	l = l1;
/*<         l1 = l+1 >*/
	l1 = l + 1;
/*<         go to 20 >*/
	goto L20;
/*<   30    k = ky1 >*/
L30:
	k = ky1;
/*<         k1 = k+1 >*/
	k1 = k + 1;
/*<   40    if(yi.lt.ty(k1) .or. k.eq.nk1y) go to 50 >*/
L40:
	if (yi < ty[k1] || k == nk1y) {
	    goto L50;
	}
/*<         k = k1 >*/
	k = k1;
/*<         k1 = k+1 >*/
	k1 = k + 1;
/*<         go to 40 >*/
	goto L40;
/*<   50    num = (l-kx1)*nyy+k-ky >*/
L50:
	num = (l - kx1) * nyy + k - *ky;
/*<         nummer(im) = index(num) >*/
	nummer[im] = index[num];
/*<         index(num) = im >*/
	index[num] = im;
/*<   60  continue >*/
/* L60: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* fporde_ */

