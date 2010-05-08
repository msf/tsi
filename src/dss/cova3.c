#include <math.h>
#include "dss.h"
#include "dss_legacy.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)
#define PI	3.14159265


/* ----------------------------------------------------------------------- */
/*                    Covariance Between Two Points */
/*                    ***************************** */
/* model specified by a nugget effect and nested varigoram structures. */
/* The anisotropy definition can be different for each nested structure. */

/* INPUT VARIABLES: */
/*   x1,y1,z1         coordinates of first point */
/*   x2,y2,z2         coordinates of second point */
/*   nst(ivarg)       number of nested structures (maximum of 4) */
/*   ivarg            variogram number (set to 1 unless doing cokriging */
/*                       or indicator kriging) */
/*   MAXNST           size of variogram parameter arrays */
/*   c0(ivarg)        isotropic nugget constant */
/*   it(i)            type of each nested structure: */
/*                      1. spherical model of range a; */
/*                      2. exponential model of parameter a; */
/*                           i.e. practical range is 3a */
/*                      3. gaussian model of parameter a; */
/*                           i.e. practical range is 3a */
/*                      4. power model of power a (a must be gt. 0  and */
/*                           lt. 2).  if linear model, a=1,c=slope. */
/*                      5. hole effect model */
/*   cc(i)            multiplicative factor of each nested structure. */
/*                      (sill-c0) for spherical, exponential,and gaussian */
/*                      slope for linear model. */
/*   aa(i)            parameter "a" of each nested structure. */
/*   irot             index of the rotation matrix for the first nested */
/*                    structure (the second nested structure will use */
/*                    irot+1, the third irot+2, and so on) */
/*   MAXROT           size of rotation matrix arrays */
/*   rotmat           rotation matrices */

/* OUTPUT VARIABLES: */
/*   cmax             maximum covariance */
/*   cova             covariance between (x1,y1,z1) and (x2,y2,z2) */

/* EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance */
/*                      rotmat    computes rotation matrix for distance */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 * sqdist
 */

/** CUBOS utilizados
 * -
 */ 

/** CUBOS _nao_ utilizados
 * sim
 * tmp
 * order
 * clc
 * lvm
 * e os cov*
 */

/** structs globais utilizadas:
 * krigev_1.rotmat
 * cova3d_1
 */

/* these variables used to be passed has arguments, but have _allways_ these values */
#define MAXROT	5

/* Calculate the maximum covariance value (used for zero distances and */
/* for power model covariance): */
double cova3(float x1, float y1, float z1, float x2, float y2, float z2,
		int varnum, float nugget, variogram_t *variogram,
		double rotmat[5][3][3], double *cmax)
{
	/* Local variables */
	double cova;
	double h, hr;
	int ir, is;
	double hsqd;
	int istart;


	/* Parameter adjustments */
	--variogram;


	/* Function Body */
	istart = 1;
	*cmax = nugget;
	for (is = 1; is <= varnum; ++is) {
		if (variogram[is].type == 4) {
			*cmax += 999.f;
		} else {
			*cmax += variogram[is].cov;
		}
	}

	/* Check for "zero" distance, return with cmax if so: */

	hsqd = sqdist(x1, y1, z1, x2, y2, z2, 1, rotmat);
	if ((float) hsqd < 1e-10f) {
		cova = *cmax;
		return cova;
	}

	/* Loop over all the structures: */

	cova = 0;
	for (is = 1; is <= varnum; ++is) {

		/* Compute the appropriate distance: */

		if (is != 1) {
			/* Computing MIN */
			ir = MIN(is,MAXROT); /* min(temp2,maxrot); */
			hsqd = sqdist(x1, y1, z1, x2, y2, z2, ir, rotmat);
		}
		h = sqrt(hsqd);

		/* Spherical Variogram Model? */

		// FIXME: XXX BUG ON SPHERICAL TYPE
		if (variogram[is].type == VARIOGRAM_TYPE_SPHERICAL) {
			hr = h / variogram[is].aa;
			if (hr < 1) {
				cova += variogram[is].cov * (1 - hr * (1.5 - hr * 0.5 * hr));
			}

		/* Exponential Variogram Model? */
		} else if (variogram[is].type == VARIOGRAM_TYPE_EXPONENCIAL) {
			cova +=  variogram[is].cov * exp(h * -3 /  variogram[is].aa);

		/* Gaussian Variogram Model? */
		} else if ( variogram[is].type == VARIOGRAM_TYPE_GAUSSIAN) {
			/* Computing 2nd power */
			double tmp = h /  variogram[is].aa;
			cova +=  variogram[is].cov * exp(tmp * tmp * -3);

		/* Power Variogram Model? */
		} else if (variogram[is].type == VARIOGRAM_TYPE_POWER) {
			double t = variogram[is].aa;
			cova += *cmax - variogram[is].cov * pow(h, t);

		/* Hole Effect Model? */
		} else if (variogram[is].type == VARIOGRAM_TYPE_HOLE) {
			/*                 d = 10.0 * aa(ist) */
			/*                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI) */
			cova += variogram[is].cov * cos(h / variogram[is].aa * PI );
		}
	}

	/* Finished: */

	return cova;
} /* cova3 */

