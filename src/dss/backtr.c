#include <math.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


float gcum(float );

int locate(float *xx, int n, int is, int ie, float x);

float powint(float xlow, float xhigh, float ylow, float yhigh, float xval, float power);

/* -----------------------------------------------------------------------
 *
 * Evaluate the standard normal cdf given a normal deviate x.  gcum is 
 * the area under a unit normal curve to the left of x.  The results are 
 * accurate only to about 5 decimal places. 
 *
 * ----------------------------------------------------------------------- */
float gcum(float z)
{
	double ret_val;

	double t;
    double zd;
        
    zd = (double) z;
        

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
	return (float)ret_val;

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
int locate(float *xx, int n, int is, int ie, float x)
{
	int jl, jm, ju;

	int i;

	/* Initialize lower and upper methods: */

	/* Parameter adjustments */
	--xx;
	i = is;
	/* Function Body */
	if (i <= 0) {
		i = 1;
	}
	jl = i - 1;
	ju = ie;
	if (xx[n] <= x) {
		/* out of range */
		return ie;
	}

	/* If we are not done then compute a midpoint: */
	while (ju - jl > 1) {
		jm = (ju + jl) / 2;

		/* Replace the lower or upper limit with the midpoint: */
		if ((xx[ie] > xx[i]) == (x > xx[jm])) {
			jl = jm;
		} else {
			ju = jm;
		}
	}

	/* Return the array index: */
	return jl;
} /* locate_ */





/* ----------------------------------------------------------------------- 
 *
 * Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh) 
 *                 for a value of x and a power pow. 
 *
 * ----------------------------------------------------------------------- */
float powint(float xlow, float xhigh, float ylow, float yhigh, float xval, float power)
{
	/* System generated locals */
	double ret_val;

	/* Builtin functions */

	if ((xhigh - xlow) < 1e-20f) {
		ret_val = (double) ((yhigh + ylow) / 2.f);
	} else {
		ret_val = (double) (ylow + (yhigh - ylow) * pow((double) ((xval - xlow) / (xhigh - xlow)), (double) (power)));
	}
	return (float)ret_val;

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
 *   min_value,max_value        limits possibly used for linear or power model 
 *
 *
 * ----------------------------------------------------------------------- */
float backtr(float vrgs, int nt, float *vr, float *vrg, float min_value, float max_value)
{
	int j;
	float cdfhi, cdfbt, cdflo;
	float ret_val;

	float a, b;
	/* parameter(EPSLON=1.0e-20) */

	/* Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid): */


	/* Parameter adjustments */
	--vrg;
	--vr;

	b = vrgs;
	/* Function Body */
	if (vrgs <= vrg[1]) {
		ret_val = vr[1];
		a = vrg[1];
		cdflo = (float) gcum(a);
		cdfbt = (float) gcum(b);
		/* 
		if (LTAIL == 1) {
			ret_val = powint(0, cdflo, min_value, vr[1], cdfbt, 1);
		} else if (LTAIL == 2) {
			cpow = 1 / min_value;
			ret_val = powint(0, cdflo, min_value, vr[1], cdfbt, cpow);
		}
		*/
		ret_val = powint(0, cdflo, min_value, vr[1], cdfbt, 1);

		/* Value in the upper tail?     1=linear, 2=power, 4=hyperbolic: */

	} else if (vrgs >= vrg[nt]) {
		ret_val = vr[nt];
		a = vrg[nt];
		cdfhi = (float) gcum(a);
		cdfbt = (float) gcum(b);
		/*
		if (UTAIL == 1) {
			ret_val = powint(cdfhi, 1, vr[nt], max_value, cdfbt, 1);
		} else if (UTAIL == 2) {
			float cpow;
			cpow = 1.f / max_value;
			ret_val = powint(cdfhi, 1, vr[nt], max_value, cdfbt, cpow);
		} else if (UTAIL == 4) {
			double lambda;
			d__1 = (double) vr[nt];
			d__2 = (double) (max_value);
			lambda =  pow(d__1, d__2) * (1 - gcum(a));
			d__1 = (double) (lambda / (1 - gcum(b)));
			d__2 = (double) (1 / max_value);
			ret_val = (float) pow(d__1, d__2);
		} 
		*/
		ret_val = powint(cdfhi, 1, vr[nt], max_value, cdfbt, 1);
	} else {

		/* Value within the transformation table: */
		int i,t;

		j = locate(&vrg[1], nt, 1, nt, vrgs);
		i = nt - 1;
		t = MIN(i,j); /* min(i__2,j) */
		j = MAX( t, 1); /* max(i__1,1); */
		ret_val = powint( vrg[j], vrg[j + 1], vr[j], vr[j + 1], vrgs, 1);
	}
	return ret_val;
} /* backtr_ */

