
/*
	wrapper functions needed in the Modelica interface
	
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "ModelicaUtilities.h"
#include "f2c.h"
#include "curfit.P"
#include "splev.P"
#include "percur.P"
#include "splder.P"
#include "surfit.P"
#include "bispev.P"
#include "parder.P"
#include "regrid.P"

/* define maximum number of derivative grid points to evaluated at once */
#define MX_MAX 10
#define MY_MAX 10

/* index of row-wise stored 2-D zero-based array */ 
#define IDXR(ir,ic,nr,nc) ((ir)*(nc)+(ic))

/* index of column-wise stored 2-D zero-based array */ 
#define IDXC(ir,ic,nr,nc) ((ic)*(nr)+(ir))

/* Modelica style index of 2-D array (one-based) */ 
#define IDXM(ir,ic,nc) (((ir)-1)*(nc)+((ic)-1))


double mapPeriodic(double *x_lim, double x);
int surf2dLengthOfWrk1(int nxest, int nyest, int kx, int ky, int m) ;
int surf2dLengthOfWrk2(int nxest, int nyest, int kx, int ky, int m) ;


typedef struct { /* 1-D spline structure */ 
	int k;		/* order of spline */
	int n;		/* total number of knots of the spline approximation returned */
	int nest;	/* over-estimate of the total number of knots  (length of arrays t,c and iwrk) */ 
	double fp;	/* weighted sum of squared residuals of the spline approximation */
	double *t;	/* this array will contain the knots of the spline */
	double *c;	/* this array will contain the coefficients c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x) */
	double x_lim[2]; /* boundaries of interpolation interval */
	double *wrk;/* double working array of length */
	int isPeriodic; /* set true for periodic splines */
} Curve1d;


typedef struct { /* 2-D spline surface */ 
	int kx,ky;		/* order of spline */
	int nx,ny;		/* total number of knots of the spline approximation returned */
	int nxest,nyest;/* over-estimate of the total number of knots  (length of arrays t,c and iwrk) */ 
	double fp;	/* weighted sum of squared residuals of the spline approximation */
	double *tx,*ty;	/* this arrays will contain the knots of the spline */
	double *c;	/* this array will contain the coefficients c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x) */
	double x_lim[2]; /* boundaries of interpolation interval */
	double y_lim[2]; /* boundaries of interpolation interval */
	int lwrk;		/* length of double working array wrk */
	int kwrk;		/* length of int working array iwrk */
	double *wrk;	/* double work array (used in derivative calculation) */
	int *iwrk;		/* int work array (used in derivative calculation) */
} Surf2d;

void indexx(int n, double arrin[], int indx[])
{
        int l,j,ir,indxt,i;
        double q;

		arrin--;
		indx--;

        for (j=1;j<=n;j++) indx[j]=j;
        if (n == 1) return;
        l=(n >> 1) + 1;
        ir=n;
        for (;;) {
                if (l > 1)
                        q=arrin[(indxt=indx[--l])];
                else {
                        q=arrin[(indxt=indx[ir])];
                        indx[ir]=indx[1];
                        if (--ir == 1) {
                                indx[1]=indxt;
                                return;
                        }
                }
                i=l;
                j=l << 1;
                while (j <= ir) {
                        if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
                        if (q < arrin[indx[j]]) {
                                indx[i]=indx[j];
                                j += (i=j);
                        }
                        else j=ir+1;
                }
                indx[i]=indxt;
        }
}
void swapDouble(double *x, double *y) {
	double s1;
	s1 = *x;
	*x = *y;
	*y = s1;
}

double meanValue(double *x, int n) {
	double y = 0;
	int i;
	for(i=0;i<n;i++) {
		y += x[i];
	}
	return y/n;
}

double standardDeviation(double *x, int n) {
	int i;
	double dx,sx2 = 0;
	double xbar = meanValue(x,n);
	for (i=0;i<n;i++) {
		dx = (x[i] - xbar);
		sx2 += dx*dx;
	}
	return sqrt(sx2/(n - 1));
}
void transposeDouble(double *A, int nr, int nc) { 
	/* transpose array A of dimension [n x m] */
	int ir,ic;
	for (ir=0;ir<nr;ir++) {
		for (ic=0;ic<nc;ic++) {
			if(ic==ir) break;
			swapDouble(A+IDXC(ir,ic,nr,nc), A+IDXR(ir,ic,nr,nc));
		}
	}
}

void printDoubleArrayC(char *name, double *A, int nr, int nc) {
	int ir,ic;
	ModelicaFormatMessage("\n\ncolumn-wise: %s = {\n", name);
	for (ir=0;ir<nr;ir++) {
		ModelicaFormatMessage("{");
		for (ic=0;ic<nc;ic++) {
			ModelicaFormatMessage("%f", A[IDXC(ir,ic,nr,nc)]);
			ModelicaFormatMessage("%s", ic<(nc-1) ? " " : "}\n");
		}
	}
	ModelicaFormatMessage("}\n\n");
}
void printDoubleArrayR(char *name, double *A, int nr, int nc) {
	int ir,ic;
	ModelicaFormatMessage("\n\nrow-wise %s = {\n", name);
	for (ir=0;ir<nr;ir++) {
		ModelicaFormatMessage("{");
		for (ic=0;ic<nc;ic++) {
			ModelicaFormatMessage("%f", A[IDXR(ir,ic,nr,nc)]);
			ModelicaFormatMessage("%s", ic<(nc-1) ? " " : "}\n");
		}
	}
	ModelicaFormatMessage("}\n\n");
}

void printDoubleVector(char *name, double *v, int n) {
	int i;
	ModelicaFormatMessage("\n\n%s = {", name);
	for (i=0;i<n;i++) {
		ModelicaFormatMessage("%f", v[i]);
		ModelicaFormatMessage("%s", i<(n-1) ? " " : "}\n\n");
	}	
}

void separateMdcArray(double *data, int nr, int nc, double *x, double *y, double *z) {
	/* separate vectors x, y  and array z  from Modelica-Style data table */
	/* Modelica data table is row-major, result data table z is column-major */
	/* x, y and z must point to suitable memory */
	int ic,ir;

	for (ir=0;ir<nr;ir++) { x[ir] = data[IDXR(ir+1,0,nr+1,nc+1)]; }
	for (ic=0;ic<nc;ic++) { y[ic] = data[IDXR(0,ic+1,nr+1,nc+1)]; }

	for (ir=0;ir<nr;ir++) {
		for (ic=0;ic<nc;ic++) {
			z[IDXR(ir,ic,nr,nc)] = data[IDXR(ir+1,ic+1,nr+1,nc+1)];
		}
	}
}


int dblcmp(const void *p1, const void *p2) {
	double *x1, *x2;
	int y;
	x1 = (double *) p1;
	x2 = (double *) p2;
	if(*x1 < *x2) {
		y = -1;
	} else if(*x1 > *x2) {
		y = 1;
	} else {
		y = 0;
	}
	return y;
}
void *curve1dNonPeriodicNew(double *data, int m, int nc, double *x_lim, int k, double s, double *t, int n) {

	Curve1d *spl;
	int i,iopt,ierr,lwrk,kwrk,*iwrk;
	double *x,*y,*w,*wrk;

	/* allocate data arrays */
	x = (double *) calloc(m,sizeof(double));
	y = (double *) calloc(m,sizeof(double));
	w = (double *) calloc(m,sizeof(double));
	if ((!x) || (!y) || (!w)) ModelicaError("Out of memory while allocating data arrays");

	/* as data is ordered row-major and read-only, we have to copy the whole data into arrays */ 
	switch (nc) {
		case 2: 
			qsort(data, m, 2*sizeof(double), dblcmp);
			for(i=0;i<m;i++) {
				x[i] = data[2*i];
				y[i] = data[2*i + 1];
				w[i] = 1.0;
			}
			break;
		case 3:
			qsort(data, m, 3*sizeof(double), dblcmp);
			for(i=0;i<m;i++) {
				x[i] = data[3*i];
				y[i] = data[3*i + 1];
				w[i] = data[3*i + 2];
			}
			break;
		default:
			ModelicaFormatError("data array must have 2 or 3 columns but it has %d", nc);
	}

	/* scale s with standard deviation of x values */
	//s = s*standardDeviation(x,m);

	/* allocate curve1d structure */
	spl = (Curve1d *) calloc(1,sizeof(Curve1d));
	if (!spl) ModelicaError("Out of memory while allocating spline struct");

	/* set degree of spline */
	spl->k = k;

	/* mark spline as non-periodic */
	spl->isPeriodic = 0;

	/* calculate array sizes */
	spl->nest = m + k + 1 ;
	lwrk = (m*(k+1)+spl->nest*(7+3*k));
	kwrk = spl->nest;

	/* allocate result arrays */
	spl->c = (double *) calloc(spl->nest,sizeof(double));
	spl->t = (double *) calloc(spl->nest,sizeof(double));
	if ((!spl->c) || (!spl->t)) ModelicaError("Out of memory while allocating result arrays");

	if(n>0){
		/* if knots are given, we generate a least square spline */
		iopt = -1;
		spl->n = n;
		for (i=0;i<n;i++) 
			spl->t[i] = t[i];
	} else {
		/* knots are not given, we generate a smoothing splines */
		iopt = 0;
		spl->n = 0;
	}

	/* copy interpolation interval boundaries */
	spl->x_lim[0] = x_lim[0];
	spl->x_lim[1] = x_lim[1];

	/* allocate working arrays (local) */
	wrk = (double *) calloc(lwrk,sizeof(double));
	iwrk = (int *)   calloc(kwrk,sizeof(int));
	if ((!wrk) || (!iwrk)) ModelicaError("Out of memory while allocating local work arrays");

	/* allocate working array needed in splder */
	spl->wrk = (double *) calloc(spl->nest,sizeof(double));
	if (!spl->wrk) ModelicaError("Out of memory while allocating work array for splder");

	/* the real thing: call the FORTRAN subroutine */
core:
	curfit_(&iopt, &m, x, y, w, x_lim, x_lim+1, &(spl->k), &s, &(spl->nest), &(spl->n),
		spl->t, spl->c, &(spl->fp), wrk, &lwrk, iwrk, &ierr); 
	
	if (ierr>0){
		int i;	

		if ((ierr<=3) && (s>0)) {
			s = 2*s;
			ModelicaFormatMessage("\n\ncurfit: s probably too small, trying s=%f\n", s);
			goto core;
		}


		ModelicaFormatMessage("\n\ncurfit: iopt=%d\n", iopt);
		if ((iopt==-1) && (n<2*k+2) || (n>min(spl->nest,m+k+1))) {
			ModelicaFormatMessage("curfit: n=%d (range: %d ... %d)\n", n, 2*k+2, min(spl->nest,m+k+1));
		} else {
			ModelicaFormatMessage("curfit: n=%d\n", spl->n);
		}
		ModelicaFormatMessage("curfit: m=%d\n",  m);
		ModelicaFormatMessage("curfit: k=%d\n", k);
		ModelicaFormatMessage("curfit: nest=%d\n", spl->nest);
		ModelicaFormatMessage("curfit: fp=%f\n", spl->fp);
		for(i=0;i<m;i++) {
			ModelicaFormatMessage("curfit: x[%d]=%f  y[%d]=%f  w[%d]=%f\n", i, x[i], i, y[i], i, w[i]);
		}
		ModelicaFormatMessage("curfit: x_b=%f  x_e=%f\n", x_lim[0], x_lim[1]);

		if (ierr==10) {
			int s;
			ModelicaFormatMessage("%s\n", "check input paramters:");
			ModelicaFormatMessage("\t 1<=k<=5 : %s\n", ((1<=k) && (k<=5)) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t m>k : %s\n", (m>k) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t nest>2*k+2 : %s\n", (spl->nest>(2*k+2)) ? "OK" : "FAIL");
			for (s=i=0;i<m;i++) { if (w[i] <= 0)	s = 1; }
			ModelicaFormatMessage("\t w(i)>0 : %s\n", (s==0) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t x_lim[1]<=x[1] : %s\n", (x_lim[0]<=x[0]) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t x_lim[2]>=x[m] : %s\n", (x_lim[1]>=x[m-1]) ? "OK" : "FAIL");
			for (s=i=0;i<(m-1);i++) { if (x[i] >= x[i+1])	s = 1; }
			ModelicaFormatMessage("\t x[i]<x[i+1] : %s\n", (s==0) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t lwrk>=(k+1)*m+nest*(7+3*k) : %s\n",
				lwrk>=((k+1)*m+spl->nest*(7+3*k)) ? "OK" : "FAIL");
			if (iopt==-1) 
				ModelicaFormatMessage("\t 2*k+2<=n<=min(nest,m+k+1) : %s\n",
					((2*k+2) <= n) && (n <= min(spl->nest,m+k+1)) ?  "OK" : "FAIL");
		}
		fflush(NULL);
		ModelicaFormatError("\n core routine 'curfit' returns err: ier=%d \n", ierr);
	} 

	free(x);
	free(y);
	free(w);
	free(wrk);
	free(iwrk);

	return (void *) spl;

}

void *curve1dPeriodicNew(double *data, int m, int nc, int k, double s, double *t, int n) {

	Curve1d *spl;
	int i,iopt,ierr,lwrk,kwrk,*iwrk;
	double *x,*y,*w,*wrk;

	/* allocate data arrays */
	x = (double *) calloc(m,sizeof(double));
	y = (double *) calloc(m,sizeof(double));
	w = (double *) calloc(m,sizeof(double));
	if ((!x) || (!y) || (!w)) ModelicaError("Out of memory while allocating data arrays");

	/* as data is ordered row-major and read-only, we have to copy the whole data into arrays */ 
	switch (nc) {
		case 2:
			qsort(data, m, 2*sizeof(double), dblcmp);
			for(i=0;i<m;i++) {
				x[i] = data[2*i];
				y[i] = data[2*i + 1];
				w[i] = 1.0;
			}
			break;
		case 3:
			qsort(data, m, 3*sizeof(double), dblcmp);
			for(i=0;i<m;i++) {
				x[i] = data[3*i];
				y[i] = data[3*i + 1];
				w[i] = data[3*i + 2];
			}
			break;
		default:
			ModelicaFormatError("data array must have 2 or 3 columns but it has %d", nc);
	}

	/* scale s with standard deviation of x values */
	//s = s*standardDeviation(x,m);

	/* allocate curve1d structure */
	spl = (Curve1d *) calloc(1,sizeof(Curve1d));
	if (!spl) ModelicaError("Out of memory while allocating spline struct");

	/* set degree of spline */
	spl->k = k;

	/* mark spline as periodic */
	spl->isPeriodic = 1;

	/* calculate array sizes */
	spl->nest = m + 2*k ;
	lwrk = (m*(k+1)+spl->nest*(8+5*k));
	kwrk = spl->nest;

	/* allocate result arrays */
	spl->c = (double *) calloc(spl->nest,sizeof(double));
	spl->t = (double *) calloc(spl->nest,sizeof(double));
	if ((!spl->c) || (!spl->t)) ModelicaError("Out of memory while allocating result arrays");

	if(n>0){
		/* if knots are given, we generate a least square spline */
		iopt = -1;
		spl->n = n;
		for (i=0;i<n;i++) 
			spl->t[i] = t[i];
	} else {
		/* knots are not given, we generate a smoothing splines */
		iopt = 0;
		spl->n = 0;
	}

	/* copy interpolation interval boundaries */
	spl->x_lim[0] = x[0];
	spl->x_lim[1] = x[m-1];

	/* allocate working local arrays */
	wrk = (double *) calloc(lwrk,sizeof(double));
	iwrk = (int *)   calloc(kwrk,sizeof(int));
	if ((!wrk) || (!iwrk)) ModelicaError("Out of memory while allocating work arrays");

	/* allocate working array needed in splder */
	spl->wrk = (double *) calloc(spl->nest,sizeof(double));
	if (!spl->wrk) ModelicaError("Out of memory while allocating work array for splder");

	/* the real thing: call the FORTRAN subroutine */
core:	
	percur_(&iopt, &m, x, y, w, &k, &s, &(spl->nest), &(spl->n),
		spl->t, spl->c, &(spl->fp), wrk, &(lwrk), iwrk, &ierr); 

	if (ierr>0){
		int i;		

		if ((ierr<=3) && (s>0)) {
			s = 2*s;
			ModelicaFormatMessage("\n\npercur: s probably too small, trying s=%f\n", s);
			goto core;
		}

		ModelicaFormatMessage("\n\npercur: iopt=%d\n", iopt);
		if ((iopt==-1) && (n<2*k+2) || (n>min(spl->nest,m+k+1))) {
			ModelicaFormatMessage("percur: n=%d (range: %d ... %d)\n", n, 2*k+2, min(spl->nest,m+k+1));
		} else {
			ModelicaFormatMessage("percur: n=%d\n", spl->n);
		}
		ModelicaFormatMessage("percur: m=%d\n",  m);
		ModelicaFormatMessage("percur: k=%d\n", k);
		ModelicaFormatMessage("percur: nest=%d\n", spl->nest);
		ModelicaFormatMessage("percur: fp=%f\n", spl->fp);
		for(i=0;i<m;i++) {
			ModelicaFormatMessage("percur: x[%d]=%f  y[%d]=%f  w[%d]=%f\n", i, x[i], i, y[i], i, w[i]);
		}
		ModelicaFormatMessage("percur: n=%d\n",  spl->n);
		for(i=0;i<spl->n;i++) {
			ModelicaFormatMessage("percur: c[%d]=%f  t[%d]=%f\n", i, spl->c[i], i, spl->t[i]);
		}
		if (ierr==10) {
			ModelicaFormatMessage("%s\n", "check input paramters:");
			ModelicaFormatMessage("\t 1<=k<=5 : %s\n", ((1<=k) && (k<=5)) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t m>k : %s\n", (m>k) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t nest>2*k+2 : %s\n", (spl->nest>(2*k+2)) ? "OK" : "FAIL");
			for (s=i=0;i<m;i++) { if (w[i] <= 0)	s = 1; }
			ModelicaFormatMessage("\t w(i)>0 : %s\n", (s==0) ? "OK" : "FAIL");
			for (s=i=0;i<(m-1);i++) { if (x[i] >= x[i+1])	s = 1; }
			ModelicaFormatMessage("\t x[i]<x[i+1] : %s\n", (s==0) ? "OK" : "FAIL");
			ModelicaFormatMessage("\t lwrk>=(k+1)*m+nest*(7+3*k) : %s\n",
				lwrk>=((k+1)*m+spl->nest*(7+3*k)) ? "OK" : "FAIL");
			if (iopt==-1) 
				ModelicaFormatMessage("\t 2*k+2<=n<=min(nest,m+k+1) : %s\n",
					((2*k+2) <= n) && (n <= min(spl->nest,m+k+1)) ?  "OK" : "FAIL");
		}

		fflush(NULL);
		ModelicaFormatError("\n core routine 'percur' returns err: ier=%d \n", ierr);
	} 

	free(x);
	free(y);
	free(w);
	free(wrk);
	free(iwrk);

	return (void *) spl;

}


void *curve1dNew(int isPeriodic, double *data, int m, int nc, double *x_lim, int k, double s, double *t, int n) {
	if (isPeriodic) {
		return curve1dPeriodicNew(data, m, nc, k, s, t, n);
	} else {
		return curve1dNonPeriodicNew(data, m, nc, x_lim, k, s, t, n);
	}
	return NULL; /* never reached */
}

void curve1dDel(void *obj) {

	Curve1d *spl = (Curve1d *) obj;

	if (spl == NULL) return;

	if(spl->c) free(spl->c);
	if(spl->t) free(spl->t);
	if(spl->wrk) free(spl->wrk);
	free(spl);

	return;
}


int curve1dGetNumberOfKnots(void *obj) {
	Curve1d *spl = (Curve1d *) obj;
	return (int) spl->n;
}


double curve1dEval(void *obj, double x) {
	Curve1d *spl = (Curve1d *) obj;
	double xx,y;
	int ierr,m;
	m = 1;

	/* for period splines we map the actual argument into the base period to get a periodic answer */
	if (spl->isPeriodic) {
		xx = mapPeriodic(spl->x_lim, x);
	} else {
		xx = x;
	}

	splev_(spl->t, &(spl->n), spl->c, &(spl->k), &xx, &y, &m, &ierr);

	if (ierr!=0) 
		ModelicaFormatError("\n\nsplev returns err: ierr=%d \n", ierr);

	return y;
}

double curve1dDer(void *obj, double x, double ddx) {
	Curve1d *spl = (Curve1d *) obj;
	double xx,y;
	int ierr,m,nu;
	m = 1;
	nu = 1;

	/* for period splines we map the actual argument into the base period to get a periodic answer */
	if (spl->isPeriodic) {
		xx = mapPeriodic(spl->x_lim, x);
	} else {
		xx = x;
	}

	splder_(spl->t, &(spl->n), spl->c, &(spl->k), &nu, &xx, &y, &m, spl->wrk, &ierr);

	if (ierr!=0) 
		ModelicaFormatError("\n\nsplev returns err: ierr=%d \n", ierr);

	return ddx*y;
}



double curve1dPeriodicEval(void *obj, double x) {
	/**
	DO NOT USE THIS ANY MORE USE curve1dEval instead 
	**/

	Curve1d *spl = (Curve1d *) obj;
	double xx,y;
	int ierr,m;
	m = 1;

	/* we map the actual argument into the base period to get a periodic answer */
	xx = mapPeriodic(spl->x_lim, x);

	splev_(spl->t, &(spl->n), spl->c, &(spl->k), &xx, &y, &m, &ierr);

	if (ierr!=0) 
		ModelicaFormatError("\n\nsplev returns err: ierr=%d \n", ierr);

	return y;
}

double curve1dPeriodicDer(void *obj, double x, double ddx) {
	/**
	DO NOT USE THIS ANY MORE USE curve1dEval instead 
	**/

	Curve1d *spl = (Curve1d *) obj;
	double xx,y;
	int ierr,m,nu;
	m = 1;
	nu = 1;

	/* we map the actual argument into the base period to get a periodic answer */
	xx = mapPeriodic(spl->x_lim, x);

	splder_(spl->t, &(spl->n), spl->c, &(spl->k), &nu, &xx, &y, &m, spl->wrk, &ierr);

	if (ierr!=0) 
		ModelicaFormatError("\n\nsplev returns err: ierr=%d \n", ierr);

	return ddx*y;
}

double mapPeriodic(double *x_lim, double x) {
	if(x>x_lim[1]) {
		x = fmod(x - x_lim[0], x_lim[1] - x_lim[0]) + x_lim[0];
	} else if(x<x_lim[0]) {
		x = -fmod(x_lim[1] - x, x_lim[1] - x_lim[0]) + x_lim[1];
	} 
	return x;
}



void *surf2dScatNew(double *data, int m, int nc, double *x_lim, double *y_lim, 
					int kx, int ky, double s, double *tx, double *ty, int nx, int ny) {

	Surf2d *surf;
	int i,iopt,ierr,nmax,lwrk1,lwrk2,kwrk,*iwrk;
	double *x,*y,*z,*w,*wrk1,*wrk2;
	double eps = sqrt(DBL_EPSILON);

	/* allocate data arrays */
	x = (double *) calloc(m,sizeof(double));
	y = (double *) calloc(m,sizeof(double));
	z = (double *) calloc(m,sizeof(double));
	w = (double *) calloc(m,sizeof(double));
	if ((!x) || (!y) || (!z) || (!w)) ModelicaError("Out of memory while allocating data arrays");

	/* as data is ordered row-major and read-only, we have to copy the whole data into arrays */ 
	switch (nc) {
		case 3: 
			for(i=0;i<m;i++) {
				x[i] = data[3*i];
				y[i] = data[3*i + 1];
				z[i] = data[3*i + 2];
				w[i] = 1.0;
			}
			break;
		case 4:
			for(i=0;i<m;i++) {
				x[i] = data[4*i];
				y[i] = data[4*i + 1];
				z[i] = data[4*i + 2];
				w[i] = data[4*i + 3];
			}
			break;
		default:
			ModelicaFormatError("data array must have 2 or 3 columns but it has %d", nc);
	}

	/* allocate curve1d structure */
	surf = (Surf2d *) calloc(1,sizeof(Surf2d));
	if (!surf) ModelicaError("Out of memory while allocating spline struct");

	/* set degree of spline */
	surf->kx = kx;
	surf->ky = ky;

	/* calculate array sizes */
	surf->nxest = kx+1+(int)(ceil(sqrt(m/2)));
	surf->nyest = ky+1+(int)(ceil(sqrt(m/2))) ;
	nmax = max(surf->nxest, surf->nyest);
    
	/* calculate local working array sizes */
	lwrk1 = surf2dLengthOfWrk1(surf->nxest, surf->nyest, kx, ky, m) ;
	lwrk2 = surf2dLengthOfWrk2(surf->nxest, surf->nyest, kx, ky, m) ;
	kwrk  = m+(surf->nxest - 2*kx - 1)*(surf->nyest - 2*ky - 1);


	/* allocate result arrays */
	surf->c  = (double *) calloc((surf->nxest - kx - 1)*(surf->nyest - ky - 1),sizeof(double));
	surf->tx = (double *) calloc(nmax,sizeof(double));
	surf->ty = (double *) calloc(nmax,sizeof(double));
	if ((!surf->c) || (!surf->tx) || (!surf->ty)) 
		ModelicaError("Out of memory while allocating result arrays");

	if(nx*ny>0){
		/* if knots are given, we generate a least square spline */
		iopt = -1;
		surf->nx = nx;
		surf->ny = ny;
		for (i=0;i<nx;i++) surf->tx[i] = tx[i];
		for (i=0;i<ny;i++) surf->ty[i] = ty[i];
	} else {
		/* knots are not given, we generate a smoothing splines */
		iopt = 0;
		surf->nx = 0;
		surf->ny = 0;
	}

	/* copy interpolation interval boundaries */
	surf->x_lim[0] = x_lim[0];
	surf->x_lim[1] = x_lim[1];
	surf->y_lim[0] = y_lim[0];
	surf->y_lim[1] = y_lim[1];

	/* allocate local working arrays */
	wrk1 = (double *) calloc(lwrk1,sizeof(double));
	wrk2 = (double *) calloc(lwrk2,sizeof(double));
	iwrk = (int *)   calloc(kwrk,sizeof(int));
	if ((!wrk1) || (!wrk2) || (!iwrk)) 
		ModelicaFormatError("Out of memory while allocating work arrays: lwrk1=%d  lwrk2=%d  kwrk=%d", 
			lwrk1, lwrk2, kwrk);


	/* the real thing: call the FORTRAN subroutine */
	surfit_(&iopt, &m, x, y, z, w, x_lim, x_lim+1, y_lim, y_lim+1, &(surf->kx), &(surf->ky), 
		&s, &(surf->nxest), &(surf->nyest), &nmax, &eps, &(surf->nx), surf->tx, &(surf->ny), 
		surf->ty, surf->c, &(surf->fp), wrk1, &(lwrk1), wrk2, &(lwrk2), 
		iwrk, &(kwrk), &ierr);

	if (ierr<-2) {
		ModelicaFormatMessage("\n\nsurfit: solution may be inaccurate due to rank deficiency system. Rank deficiency is %d\n",
			(surf->nx-surf->kx-1)*(surf->ny-surf->ky-1)+ierr);
	}

	if (ierr==1) {
		ModelicaFormatMessage("\n\nsurfit: approximation returned is the weighted least-squares sline according to the current set of knots.\n"
			"However, smoothing parameter s=%f cannot be fulfilled. Weighted sum of squared residuals fp=%f\n",
			s, surf->fp);
	}

	if (ierr==2) {
		ModelicaFormatMessage("\n\nsurfit: a theoretically impossible result was found during the iteration.\n"
			"Weighted sum of squared residuals does not satisfy the condition abs(fp-s)/s < 0.001: s=%f, fp=%f\n",
			s, surf->fp);
	
	}

	if (ierr==3) {
		ModelicaFormatMessage("\n\nsurfit: maximal number of iterations maxit=20 been reached.\n"
			"Weighted sum of squared residuals does not satisfy the condition abs(fp-s)/s < 0.001: s=%f, fp=%f\n",
			s, surf->fp);
	}

	if (ierr==4) {
		ModelicaFormatMessage("\n\nsurfit: no more knots can be added because the additional knot would (quasi) coincide with an old one.\n"
			"smoothing parameter s=%f cannot be fulfilled. Weighted sum of squared residuals fp=%f\n",
			s, surf->fp);
	}

	if (ierr==5) {
		ModelicaFormatMessage("\n\nsurfit: no more knots can be added because the number of b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds the number of data points m.\n"
			"smoothing parameter s=%f cannot be fulfilled. Weighted sum of squared residuals fp=%f\n",
			s, surf->fp);
	}

	if (ierr==10) {
		int s1;
		ModelicaFormatMessage("\n\nsurfit: iopt=%d\n", iopt);
		ModelicaFormatMessage("surfit: nx=%d  ny=%d\n", surf->nx, surf->ny);
		ModelicaFormatMessage("surfit: m=%d\n",  m);
		ModelicaFormatMessage("surfit: kx=%d  ky=%d\n", kx, ky);
		ModelicaFormatMessage("surfit: nxest=%d  nyest=%d\n", surf->nxest, surf->nyest);
		ModelicaFormatMessage("surfit: fp=%f\n", surf->fp);
		for(i=0;i<m;i++) {
			ModelicaFormatMessage("surfit: x[%d]=%f  y[%d]=%f  z[%d]=%f  w[%d]=%f\n", i, x[i], i, y[i], i, z[i], i, w[i]);
		}
		ModelicaFormatMessage("surfit: x_b=%f  x_e=%f  y_b=%f  y_e=%f \n", x_lim[0], x_lim[1], y_lim[0], y_lim[1]);

		ModelicaFormatMessage("\n\nsurfit: the input data are not valid:\n");
		ModelicaFormatMessage("\t 1<=kx<=5 : %s\n", ((1<=kx) && (kx<=5)) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t 1<=ky<=5 : %s\n", ((1<=ky) && (ky<=5)) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t m>=(kx+1)*(ky+1) : %s\n", (m>=(kx+1)*(ky+1)) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t nxest>=2*kx+2 : %s\n", (surf->nxest>=2*kx+2) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t nyest>=2*ky+2 : %s\n", (surf->nyest>=2*ky+2) ? "OK" : "FAIL");
		for (s1=i=0;i<m;i++) { if (w[i] <= 0)	s1 = 1; }
		ModelicaFormatMessage("\t w(i)>0 : %s\n", (s1==0) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t x_lim[1]<=x[1] : %s\n", (x_lim[0]<=x[0]) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t x_lim[2]>=x[m] : %s\n", (x_lim[1]>=x[m-1]) ? "OK" : "FAIL");
		for (s1=i=0;i<(m-1);i++) { if (x[i] >= x[i+1])	s1 = 1; }
		ModelicaFormatMessage("\t x[i]<x[i+1] : %s\n", (s1==0) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t y_lim[1]<=y[1] : %s\n", (y_lim[0]<=y[0]) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t y_lim[2]>=y[m] : %s\n", (y_lim[1]>=y[m-1]) ? "OK" : "FAIL");
		for (s1=i=0;i<(m-1);i++) { if (y[i] >= y[i+1])	s1 = 1; }
		ModelicaFormatMessage("\t y[i]<y[i+1] : %s\n", (s1==0) ? "OK" : "FAIL");
		ModelicaFormatError("core routine 'surfit' returned err: ier=%d \n", ierr);
	} 

	if (ierr>10) {
		ModelicaFormatError("\n\nsurfit: internal error: lwrk2 too small: lwrk2%d \n", lwrk2);
	}

	free(x);
	free(y);
	free(z);
	free(w);
	free(wrk1);
	free(wrk2);
	free(iwrk);

	/* calculate working arrays for derivative calculation */
	surf->lwrk = MX_MAX*(kx)+MY_MAX*(ky)+(nx-kx-1)*(ny-ky-1);
	surf->kwrk = MX_MAX+MY_MAX;

	/* allocate persistence working arrays for derivative calculation*/
	surf->wrk = (double *) calloc(surf->lwrk,sizeof(double));
	surf->iwrk = (int *)   calloc(surf->kwrk,sizeof(int));
	if ((!surf->wrk) || (!surf->iwrk)) 
		ModelicaFormatError("Out of memory while allocating persistence work arrays: lwrk=%d  kwrk=%d", 
			surf->lwrk, surf->kwrk);

	return (void *) surf;

}
void *surf2dRectNew(double *data, int mx, int my, double *x_lim, double *y_lim, 
					int kx, int ky, double s, double *tx, double *ty, int nx, int ny) {

	Surf2d *surf;
	int i,iopt,ierr,nmax,lwrk,kwrk,*iwrk;
	double *x,*y,*z,*wrk;
	double eps = sqrt(DBL_EPSILON);

	/* mx, my are the sizes of the Modelica data array but we need the lenghts of vectors x and y */
	mx--;my--;

	/* allocate data arrays */
	x = (double *) calloc(mx,sizeof(double));
	y = (double *) calloc(my,sizeof(double));
	z = (double *) calloc(mx*my,sizeof(double));
	if ((!x) || (!y) || (!z)) ModelicaError("Out of memory while allocating data arrays");

	/* as data is ordered row-major and read-only, we have to copy the whole data into arrays */ 
	separateMdcArray(data, mx, my, x, y, z);

	/* allocate curve1d structure */
	surf = (Surf2d *) calloc(1,sizeof(Surf2d));
	if (!surf) ModelicaError("Out of memory while allocating spline struct");

	/* set degree of spline */
	surf->kx = kx;
	surf->ky = ky;

	/* calculate array sizes */
	surf->nxest = mx+kx+1;
	surf->nyest = my+ky+1;
	nmax = max(surf->nxest, surf->nyest);
    
    


	/* allocate result arrays */
	surf->c  = (double *) calloc((surf->nxest - kx - 1)*(surf->nyest - ky - 1),sizeof(double));
	surf->tx = (double *) calloc(nmax,sizeof(double));
	surf->ty = (double *) calloc(nmax,sizeof(double));
	if ((!surf->c) || (!surf->tx) || (!surf->ty)) 
		ModelicaError("Out of memory while allocating result arrays");

	if(nx*ny>0){
		/* if knots are given, we generate a least square spline */
		iopt = -1;
		surf->nx = nx;
		surf->ny = ny;
		for (i=0;i<nx;i++) surf->tx[i] = tx[i];
		for (i=0;i<ny;i++) surf->ty[i] = ty[i];
	} else {
		/* knots are not given, we generate a smoothing splines */
		iopt = 0;
		surf->nx = 0;
		surf->ny = 0;
	}

	/* copy interpolation interval boundaries */
	surf->x_lim[0] = x_lim[0];
	surf->x_lim[1] = x_lim[1];
	surf->y_lim[0] = y_lim[0];
	surf->y_lim[1] = y_lim[1];

	/* calculate temporary working array sizes */
	lwrk = 4 + surf->nxest*(my + 2*kx + 5) + surf->nyest*(2*ky + 5) + 
		mx*(kx + 1) + my*(ky + 1) + max(my, surf->nxest);
	kwrk  = 3 + mx + my + surf->nxest + surf->nyest;

	/* allocate working arrays */
	wrk = (double *) calloc(lwrk,sizeof(double));
	iwrk = (int *)   calloc(kwrk,sizeof(int));
	if ((!wrk) || (!iwrk)) 
		ModelicaFormatError("Out of memory while allocating work arrays: lwrk=%d  kwrk=%d", 
			lwrk, kwrk);

	/* the real thing: call the FORTRAN subroutine */
	regrid_(&iopt, &mx, x, &my, y, z, x_lim, x_lim+1, y_lim, y_lim+1, &(surf->kx), &(surf->ky), 
		&s, &(surf->nxest), &(surf->nyest), &(surf->nx), surf->tx, &(surf->ny), 
		surf->ty, surf->c, &(surf->fp), wrk, &(lwrk),  
		iwrk, &(kwrk), &ierr);


	if (ierr==1) {
		ModelicaFormatMessage("\n\nregrid: approximation returned is the weighted least-squares sline according to the current set of knots.\n"
			"However, smoothing parameter s=%f cannot be fulfilled. Weighted sum of squared residuals fp=%f\n",
			s, surf->fp);
	}

	if (ierr==2) {
		ModelicaFormatMessage("\n\nregrid: a theoretically impossible result was found during the iteration.\n"
			"Weighted sum of squared residuals does not satisfy the condition abs(fp-s)/s < 0.001: s=%f, fp=%f\n",
			s, surf->fp);
	
	}

	if (ierr==3) {
		ModelicaFormatMessage("\n\nregrid: maximal number of iterations maxit=20 been reached.\n"
			"Weighted sum of squared residuals does not satisfy the condition abs(fp-s)/s < 0.001: s=%f, fp=%f\n",
			s, surf->fp);
	}

	if (ierr==10) {
		int s1,i,j;
		ModelicaFormatMessage("\n\nregrid: iopt=%d\n", iopt);
		ModelicaFormatMessage("regrid: nx=%d  ny=%d\n", surf->nx, surf->ny);
		ModelicaFormatMessage("regrid: mx=%d  my=%d\n",  mx, my);
		ModelicaFormatMessage("regrid: kx=%d  ky=%d\n", kx, ky);
		ModelicaFormatMessage("regrid: nxest=%d  nyest=%d\n", surf->nxest, surf->nyest);
		ModelicaFormatMessage("regrid: fp=%f\n", surf->fp);
		ModelicaFormatMessage("%s","x={");
		for (i=0;i<mx;i++) 
			ModelicaFormatMessage(" %f", x[i]);
		ModelicaFormatMessage("%s"," }\n");
		ModelicaFormatMessage("%s","y={");
		for (i=0;i<my;i++) 
			ModelicaFormatMessage(" %f", y[i]);
		ModelicaFormatMessage("%s"," }\n");
		ModelicaFormatMessage("%s","z={\n");
		for (i=0;i<mx;i++) {
			ModelicaFormatMessage("%s","{");
			for (j=0;j<my;j++) { 
				ModelicaFormatMessage(" %f", z[IDXC(i,j,mx,my)]);
			}
			ModelicaFormatMessage("%s"," }\n");
		}
		ModelicaFormatMessage("%s"," }\n");
		ModelicaFormatMessage("regrid: x_b=%f  x_e=%f  y_b=%f  y_e=%f \n", x_lim[0], x_lim[1], y_lim[0], y_lim[1]);

		ModelicaFormatMessage("\n\nregrid: the input data are not valid:\n");
		ModelicaFormatMessage("\t 1<=kx<=5 : %s\n", ((1<=kx) && (kx<=5)) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t 1<=ky<=5 : %s\n", ((1<=ky) && (ky<=5)) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t mx>=kx : %s\n", (mx>kx) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t my>=ky : %s\n", (my>ky) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t nxest>=2*kx+2 : %s\n", (surf->nxest>=2*kx+2) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t nyest>=2*ky+2 : %s\n", (surf->nyest>=2*ky+2) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t kwrk>=3+mx+my+nxest+nyest : %s\n", 
			(kwrk>=3+mx+my+surf->nxest+surf->nyest) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t lwrk>=4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+max(my,nxest) : %s\n", 
			(lwrk>=4+surf->nxest*(my+2*kx+5)+surf->nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)+max(my,surf->nxest)) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t x_lim[1]<=x[1] : %s\n", (x_lim[0]<=x[0]) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t x_lim[2]>=x[mx] : %s\n", (x_lim[1]>=x[mx-1]) ? "OK" : "FAIL");
		for (s1=i=0;i<(mx-1);i++) { if (x[i] >= x[i+1])	s1 = 1; }
		ModelicaFormatMessage("\t x[i]<x[i+1] : %s\n", (s1==0) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t y_lim[1]<=y[1] : %s\n", (y_lim[0]<=y[0]) ? "OK" : "FAIL");
		ModelicaFormatMessage("\t y_lim[2]>=y[my] : %s\n", (y_lim[1]>=y[my-1]) ? "OK" : "FAIL");
		for (s1=i=0;i<(my-1);i++) { if (y[i] >= y[i+1])	s1 = 1; }
		ModelicaFormatMessage("\t y[i]<y[i+1] : %s\n", (s1==0) ? "OK" : "FAIL");
		ModelicaFormatError("core routine 'regrid' returned err: ier=%d \n", ierr);

		if (iopt==-1) {
			ModelicaFormatMessage("\t 2*kx+2<=nx : %s", 
				(2*kx+2)<=surf->nx ? "OK" : "FAIL");
			ModelicaFormatMessage("\t nx<=min(nxest,mx+kx+1) : %s", 
				surf->nx<=min(surf->nxest,mx+kx+1) ? "OK" : "FAIL");
		}
	} 

	free(x);
	free(y);
	free(z);
	free(wrk);
	free(iwrk);

	/* calculate working arrays for derivative calculation */
	surf->lwrk = MX_MAX*(kx)+MY_MAX*(ky)+(nx-kx-1)*(ny-ky-1);
	surf->kwrk = MX_MAX+MY_MAX;

	/* allocate persistence working arrays for derivative calculation*/
	surf->wrk = (double *) calloc(surf->lwrk,sizeof(double));
	surf->iwrk = (int *)   calloc(surf->kwrk,sizeof(int));
	if ((!surf->wrk) || (!surf->iwrk)) 
		ModelicaFormatError("Out of memory while allocating persistence work arrays: lwrk=%d  kwrk=%d", 
			surf->lwrk, surf->kwrk);



	return (void *) surf;

}


void surf2dDel(void *obj) {
	Surf2d *s = (Surf2d *) obj;

	if (s == NULL) return;

	if(s->c) free(s->c);
	if(s->tx) free(s->tx);
	if(s->ty) free(s->ty);
	if(s->wrk) free(s->wrk);
	if(s->iwrk) free(s->iwrk);

	free(s);

	return;
}


int surf2dLengthOfWrk1(int nxest, int nyest, int kx, int ky, int m) {
	int u,v,km,ne,bx,by,b1,b2;
	u = nxest-kx-1; 
	v = nyest-ky-1;
	km = max(kx,ky)+1;
    ne = max(nxest,nyest);
	bx = kx*v+ky+1;
	by = ky*u+kx+1;
	if(bx<=by) { 
		b1 = bx;
		b2 = b1+v-ky;
	}
	if(bx>by) {
		b1 = by;
		b2 = b1+u-kx;
	}
    return u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1;
}

int surf2dLengthOfWrk2(int nxest, int nyest, int kx, int ky, int m) {
	int u,v,km,ne,bx,by,b1,b2;
	u = nxest-kx-1; 
	v = nyest-ky-1;
	km = max(kx,ky)+1;
    ne = max(nxest,nyest);
	bx = kx*v+ky+1;
	by = ky*u+kx+1;
	if(bx<=by) { 
		b1 = bx;
		b2 = b1+v-ky;
	}
	if(bx>by) {
		b1 = by;
		b2 = b1+u-kx;
	}
    return u*v*(b2+1)+b2;
}

double surf2dEval(void *obj, double x, double y) {
	Surf2d *s = (Surf2d *) obj;
	double z;
	int ierr = 0;
	int mx = 1;
	int my = 1;

	bispev_(s->tx, &(s->nx), s->ty, &(s->ny), s->c, &(s->kx), &(s->ky),&x, &mx, &y, &my, &z, 
		s->wrk, &(s->lwrk), s->iwrk, &(s->kwrk), &ierr);

	if (ierr!=0) 
		ModelicaFormatError("\n\nsplev returns err: ierr=%d \n", ierr);

	return z;
}


double surf2dDer(void *obj, double x, double y, double dx, double dy) {
	Surf2d *s = (Surf2d *) obj;
	double dzdx,dzdy;
	int ierr;
	int nux, nuy, mx, my;
	mx = my = 1;

	nux = 1; nuy = 0;
	parder_(s->tx, &(s->nx), s->ty, &(s->ny), s->c, &(s->kx), &(s->ky), &nux, &nuy, &x, &mx,
		&y, &my, &dzdx, s->wrk, &(s->lwrk), s->iwrk, &(s->kwrk), &ierr);

	nux = 0; nuy = 1;
	parder_(s->tx, &(s->nx), s->ty, &(s->ny), s->c, &(s->kx), &(s->ky), &nux, &nuy, &x, &mx,
		&y, &my, &dzdy, s->wrk, &(s->lwrk), s->iwrk, &(s->kwrk), &ierr);

	if (ierr!=0) 
		ModelicaFormatError("\n\nparder returns err: ierr=%d \n", ierr);
	
	return dx*dzdx + dy*dzdy;
}

void surf2dGetNumberOfKnots(void *obj, int *nx, int *ny) {
	Surf2d *s = (Surf2d *) obj;
	*nx = s->nx;
	*ny = s->ny;
}
int surf2dGetNumberOfKnotsX(void *obj) {
	Surf2d *s = (Surf2d *) obj;
	return s->nx;
}

int surf2dGetNumberOfKnotsY(void *obj) {
	Surf2d *s = (Surf2d *) obj;
	return s->ny;
}

void *surf2dNew(int isRect, double *data, int mx, int my, double *x_lim, double *y_lim, 
				int kx, int ky, double s, double *tx, double *ty, int nx, int ny) {
					if (isRect) {
						return surf2dRectNew(data, mx, my, x_lim, y_lim, kx, ky, s, tx, ty, nx, ny);
					} else {
						return surf2dScatNew(data, mx, my, x_lim, y_lim, kx, ky, s, tx, ty, nx, ny);
					}
					return NULL; /* never reached */
}
