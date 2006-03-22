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
float gcum(float z)
{
	/* System generated locals */
	float ret_val;



	/* Local variables */
	/*    static float t, z, e2; */
	float t;


	/* LPL: old code */
	/*    z = *x;

		  if (z < 0.f) {    //ub
		  z = -z;
		  }

		  t = 1.f / (z * .2316419f + 1.f);
		  ret_val = t * (t * (t * (t * (t * 1.330274429f - 1.821255978f) + 
		  1.781477937f) - .356563782f) + .31938153f);
		  e2 = 0.f;
		  */
	/*  6 standard deviations out gets treated as infinity: */

	/*    if (z <= 6.f) {       //ub
		  e2 = exp(-z * z / 2.f) * .3989422803f;
		  }

		  ret_val = 1.f - e2 * ret_val;

		  if (*x >= 0.f) {          //ub
		  return ret_val;
		  }
		  ret_val = 1.f - ret_val;
		  return ret_val;
		  */

	/* LPL: new code */
	if (z < 0.f) {   /* unpredictable branch */
		if (z >= -6.f) {   /* unpredictable branch */
			t = 1.f / (z * -.2316419f + 1.f);
			ret_val = t * (t * (t * (t * (t * 1.330274429f - 1.821255978f) + 
							1.781477937f) - .356563782f) + .31938153f);
			ret_val = exp(-z * z / 2.f) * .3989422803f * ret_val;
		} else {
			ret_val = 0.f;
		}
	} else {
		if (z <= 6.f) {   /* unpredictable branch */
			t = 1.f / (z * .2316419f + 1.f);
			ret_val = t * (t * (t * (t * (t * 1.330274429f - 1.821255978f) + 
							1.781477937f) - .356563782f) + .31938153f);
			ret_val = 1.f - exp(-z * z / 2.f) * .3989422803f * ret_val;
		} else {
			ret_val = 1.f;
		}
	}
	return ret_val;

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
	/* double d__1, d__2;*/

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
        float c_b2 = 0.f;
        float c_b3 = 1.f;
        int one = 1;

	/* System generated locals */
	int i__1, i__2;
	float ret_val;
	double d__1, d__2;

	/* Local variables */
	int j;
	float cpow, cdfhi, cdfbt, cdflo, lambda;

	/* parameter(EPSLON=1.0e-20) */

	/* Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid): */


	/* Parameter adjustments */
	--vrg;
	--vr;

	/* Function Body */
	if (*vrgs <= vrg[1]) {
		ret_val = vr[1];
		cdflo = gcum(vrg[1]);
		cdfbt = gcum(*vrgs);
		if (*ltail == 1) {
			ret_val = powint(&c_b2, &cdflo, zmin, &vr[1], &cdfbt, &c_b3);
		} else if (*ltail == 2) {
			cpow = 1.f / *ltpar;
			ret_val = powint(&c_b2, &cdflo, zmin, &vr[1], &cdfbt, &cpow);
		}

		/* Value in the upper tail?     1=linear, 2=power, 4=hyperbolic: */

	} else if (*vrgs >= vrg[*nt]) {
		ret_val = vr[*nt];
		cdfhi = gcum(vrg[*nt]);
		cdfbt = gcum(*vrgs);
		if (*utail == 1) {
			ret_val = powint(&cdfhi, &c_b3, &vr[*nt], zmax, &cdfbt, &c_b3);
		} else if (*utail == 2) {
			cpow = 1.f / *utpar;
			ret_val = powint(&cdfhi, &c_b3, &vr[*nt], zmax, &cdfbt, &cpow);
		} else if (*utail == 4) {
			d__1 = (double) vr[*nt];
			d__2 = (double) (*utpar);
			lambda = pow(d__1, d__2) * (1.f - gcum(vrg[*nt]));
			d__1 = (double) (lambda / (1.f - gcum(*vrgs)));
			d__2 = (double) (1.f / *utpar);
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

