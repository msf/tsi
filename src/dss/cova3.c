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

double cova3(float x1, float y1, float z1, float x2, float y2, float z2,
		int nst, float c0, int *it, float *cc,
		float *aa, double rotmat[5][3][3], double *cmax)
/* Calculate the maximum covariance value (used for zero distances and */
/* for power model covariance): */
{
	/* System generated locals */
	int temp;
	
	/* Local variables */
	double cova;
	double h, hr;
	int ir, is;
	double hsqd;
	int istart;


	/* Parameter adjustments */
	--it;
	--cc;
	--aa;


	/* Function Body */
	istart = 1;
	*cmax = c0;
	temp = nst;
	for (is = 1; is <= nst; ++is) {
		if (it[is] == 4) {
			*cmax += 999.f;
		} else {
			*cmax += cc[is];
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
	for (is = 1; is <= nst; ++is) {

		/* Compute the appropriate distance: */

		if (is != 1) {
			/* Computing MIN */
			ir = MIN(is,MAXROT); /* min(temp2,maxrot); */
			hsqd = sqdist(x1, y1, z1, x2, y2, z2, ir, rotmat);
		}
		h = sqrt(hsqd);

		/* Spherical Variogram Model? */

		if (it[is] == 1) {
			hr = h / aa[is];
			if (hr < 1) {
				cova += cc[is] * (1 - hr * (1.5 - hr * 0.5 * hr));
			}

			/* Exponential Variogram Model? */

		} else if (it[is] == 2) {
			cova += cc[is] * exp(h * -3 / aa[is]);

			/* Gaussian Variogram Model? */

		} else if (it[is] == 3) {
			/* Computing 2nd power */
			double tmp = h / aa[is];
			cova += cc[is] * exp(tmp * tmp * -3);

			/* Power Variogram Model? */

		} else if (it[is] == 4) {
			double t = aa[is];
			cova += *cmax - cc[is] * pow(h, t);

			/* Hole Effect Model? */

		} else if (it[is] == 5) {
			/*                 d = 10.0 * aa(ist) */
			/*                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI) */
			cova += cc[is] * cos(h / aa[is] * PI );
		}
	}

	/* Finished: */

	return cova;
} /* cova3 */

