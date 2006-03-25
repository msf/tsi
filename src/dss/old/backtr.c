#include <math.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


/* -----------------------------------------------------------------------
 *
 * Evaluate the standard normal cdf given a normal deviate x.  gcum is 
 * the area under a unit normal curve to the left of x.  The results are 
 * accurate only to about 5 decimal places. 
 *
 * ----------------------------------------------------------------------- */
double gcum(float z)
{
	/* System generated locals */
	/* float ret_val; */
	double ret_val;



	double t;
    double zd;
        
        zd = z;
        

	if (zd < 0) {   /* unpredictable branch */
		if (zd >= -6) {   /* unpredictable branch */
			t = 1 / (zd * -.2316419 + 1);
			ret_val = t * (t * (t * (t * (t * 1.330274429 - 1.821255978) + 
							1.781477937) - .356563782) + .31938153);
			ret_val = exp(-zd * zd / 2) * .3989422803 * ret_val;
		} else {
			ret_val = 0;
		}
	} else {
		if (zd <= 6) {   /* unpredictable branch */
			t = 1 / (zd * .2316419 + 1);
			ret_val = t * (t * (t * (t * (t * 1.330274429 - 1.821255978) + 
							1.781477937) - .356563782) + .31938153);
			ret_val = 1 - exp(-zd * zd / 2) * .3989422803 * ret_val;
		} else {
			ret_val = 1;
		}
	}
	return ( ret_val);

} /* gcum_ */




/* ----------------------------------------------------------------------- 
 *
 * Given an array "xx" of length "n", and given a value "x", this routine 
 * returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx 
 * must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is 
 * returned to indicate that x is out of range. 
 *
 * Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90. 
 * ----------------------------------------------------------------------- */
int locate(float *xx, int *n, int *is, int *ie, float *x, int *j)
{
	int jl, jm, ju;


	/* Initialize lower and upper methods: */

	/* Parameter adjustments */
	--xx;

	/* Function Body */
	if (*is <= 0) {
		*is = 1;
	}
	jl = *is - 1;
	ju = *ie;
	if (xx[*n] <= *x) {
		*j = *ie;
		return 0;
	}

	/* If we are not done then compute a midpoint: */
	while (ju - jl > 1) {
		jm = (ju + jl) / 2;

		/* Replace the lower or upper limit with the midpoint: */
		if ((xx[*ie] > xx[*is]) == (*x > xx[jm])) {
			jl = jm;
		} else {
			ju = jm;
		}
	}

	/* Return with the array index: */
	*j = jl;
	return 0;
} /* locate_ */





/* ----------------------------------------------------------------------- 
 *
 * Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh) 
 *                 for a value of x and a power pow. 
 *
 * ----------------------------------------------------------------------- */
double powint(float *xlow, float *xhigh, float *ylow, float *yhigh,
		float *xval, float *power)
{
	/* System generated locals */
	float ret_val;

	/* Builtin functions */

	if (*xhigh - *xlow < 1e-20f) {
		ret_val = (double) ((*yhigh + *ylow) / 2.f);
	} else {
		ret_val = (double) (*ylow + (*yhigh - *ylow) *
				(float) pow((double) ((*xval - *xlow) / (*xhigh - *xlow)),
							(double) (*power)));
	}
	return ret_val;

} /* powint_ */



/* ----------------------------------------------------------------------- 
 *
 *           Back Transform Univariate Data from Normal Scores 
 *           ************************************************* 
 *
 * This subroutine backtransforms a standard normal deviate from a 
 * specified back transform table and option for the tails of the 
 * distribution.  Call once with "first" set to true then set to false 
 * unless one of the options for the tail changes. 
 *
 * INPUT VARIABLES: 
 *
 *   vrgs             normal score value to be back transformed 
 *   nt               number of values in the back transform tbale 
 *   vr(nt)           original data values that were transformed 
 *   vrg(nt)          the corresponding transformed values 
 *   zmin,zmax        limits possibly used for linear or power model 
 *   ltail            option to handle values less than vrg(1): 
 *   ltpar            parameter required for option ltail 
 *   utail            option to handle values greater than vrg(nt): 
 *   utpar            parameter required for option utail 
 *
 *
 * ----------------------------------------------------------------------- */
double backtr(float *vrgs, int *nt, float *vr, float *vrg, float *zmin, 
		float *zmax, int *ltail, float *ltpar, int *utail, float *utpar)
{
        /* Table of constant values */
        float c_b2 = 0;
        float c_b3 = 1;
        int one = 1;

	/* System generated locals */
	int i__1, i__2;
	float ret_val;
	double d__1, d__2;

	/* Local variables */
	int j;
	float cpow, cdfhi, cdfbt, cdflo;
	double lambda;

	float a, b;
	/* parameter(EPSLON=1.0e-20) */

	/* Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid): */


	/* Parameter adjustments */
	--vrg;
	--vr;

	b = *vrgs;
	/* Function Body */
	if (*vrgs <= vrg[1]) {
		ret_val = vr[1];
		a = vrg[1];
		cdflo = (float) gcum(a);
		cdfbt = (float) gcum(b);
		if (*ltail == 1) {
			ret_val = powint(&c_b2, &cdflo, zmin, &vr[1], &cdfbt, &c_b3);
		} else if (*ltail == 2) {
			cpow = 1 / *ltpar;
			ret_val = powint(&c_b2, &cdflo, zmin, &vr[1], &cdfbt, &cpow);
		}

		/* Value in the upper tail?     1=linear, 2=power, 4=hyperbolic: */

	} else if (*vrgs >= vrg[*nt]) {
		ret_val = vr[*nt];
		a = vrg[*nt];
		cdfhi = (float) gcum(a);
		cdfbt = (float) gcum(b);
		if (*utail == 1) {
			ret_val = powint(&cdfhi, &c_b3, &vr[*nt], zmax, &cdfbt, &c_b3);
		} else if (*utail == 2) {
			cpow = 1.f / *utpar;
			ret_val = powint(&cdfhi, &c_b3, &vr[*nt], zmax, &cdfbt, &cpow);
		} else if (*utail == 4) {
			d__1 = (double) vr[*nt];
			d__2 = (double) (*utpar);
			lambda =  pow(d__1, d__2) * (1 - gcum(a));
			d__1 = (double) (lambda / (1 - gcum(b)));
			d__2 = (double) (1 / *utpar);
			ret_val = pow(d__1, d__2);
		}
	} else {

		/* Value within the transformation table: */

		locate(&vrg[1], nt, &one, nt, vrgs, &j);
		/* Computing MAX */
		/* Computing MIN */
		i__2 = *nt - 1;
		i__1 = MIN(i__2,j); /* min(i__2,j) */
		j = MAX(i__1,1); /* max(i__1,1); */
		ret_val = powint(&vrg[j], &vrg[j + 1], &vr[j], &vr[j + 1], vrgs, & c_b3);
	}
	return ret_val;
} /* backtr_ */

