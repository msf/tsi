#include <math.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* ----------------------------------------------------------------------- */
/* Computes the inverse of the standard normal cumulative distribution */
/* function with a numerical approximation from : Statistical Computing, */
/* by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95. */
/* INPUT/OUTPUT: */
/*   p    = double precision cumulative probability value: ] 0, 1 [  */
/*   xp   = G^-1 (p) in single precision */
/*   ierr = 1 - then error situation (p out of range), 0 - OK */
/* ----------------------------------------------------------------------- */


int gauinv(double p, float *xp)
{
	int rc;
	/* Initialized data */

	static double lim = 1e-10;
	static double q3 = .10353775285;
	static double q4 = .0038560700634;
	static double p0 = -.322232431088;
	static double p1 = -1.;
	static double p2 = -.342242088547;
	static double p3 = -.0204231210245;
	static double p4 = -4.53642210148e-5;
	static double q0 = .099348462606;
	static double q1 = .588581570495;
	static double q2 = .531103462366;


	/* Local variables */
	double y, pp;


	/* Coefficients of approximation: */


	/* Check for an error situation: */

	rc = 1;
	if (p < lim) {
		*xp = -1e10f;
		return rc;
	}
	if (p > 1 - lim) {
		*xp = 1e10f;
		return rc;
	}
	rc = 0;

	/* Get k for an error situation: */

	pp = p;
	if (p > .5f) {
		pp = 1 - pp;
	}
	*xp = 0.f;
	if (p == .5f) {
		return rc;
	}

	/* Approximate the function: */

	y = sqrt(log(1.f / (pp * pp)));
	*xp = (float) (y + ((((y * p4 + p3) * y + p2) * y + p1) * y + p0) / ((((y *
							q4 + q3) * y + q2) * y + q1) * y + q0));
	if ((p) ==  pp) {
		*xp = -(*xp);
	}

	/* Return with G^-1(p): */

	return rc;
} /* gauinv_ */


