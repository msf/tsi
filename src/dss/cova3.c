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
int ivarg = 1;
int maxnst = 4;
int irot = 1;
int maxrot = 5;

float cova3(float x1, float y1, float z1, float x2, float y2, float z2,
		int *nst, float *c0, int *it, float *cc,
		float *aa, double rotmat[][][], float *cmax)
/* Calculate the maximum covariance value (used for zero distances and */
/* for power model covariance): */
{
	/* System generated locals */
	int rotmat_dim1, rotmat_offset, temp, temp2;
	float r__1;
	double d__1, d__2;


	/* Local variables */
	float cova;
	float h, hr;
	int ir, is, ist;
	double hsqd;
	int istart;


	/* Parameter adjustments */
	--nst;
	--c0;
	--it;
	--cc;
	--aa;


	/* Function Body */
	istart = (ivarg - 1) * maxnst + 1;
	*cmax = c0[ivarg];
	temp = nst[ivarg];
	for (is = 1; is <= temp; ++is) {
		ist = istart + is - 1;
		if (it[ist] == 4) {
			*cmax += 999.f;
		} else {
			*cmax += cc[ist];
		}
	}

	/* Check for "zero" distance, return with cmax if so: */

	hsqd = sqdist(x1, y1, z1, x2, y2, z2, 1, rotmat);
	if ((float) hsqd < 1e-10f) {
		cova = *cmax;
		return cova;
	}

	/* Loop over all the structures: */

	cova = 0.f;
	temp = nst[ivarg];
	for (is = 1; is <= temp; ++is) {
		ist = istart + is - 1;

		/* Compute the appropriate distance: */

		if (ist != 1) {
			/* Computing MIN */
			temp2 = irot + is - 1;
			ir = MIN(temp2, maxrot); /* min(temp2,maxrot); */
			hsqd = sqdist(x1, y1, z1, x2, y2, z2, ir, rotmat);
		}
		h = (float) sqrt(hsqd);

		/* Spherical Variogram Model? */

		if (it[ist] == 1) {
			hr = h / aa[ist];
			if (hr < 1.f) {
				cova += cc[ist] * (1.f - hr * (1.5f - hr * .5f * hr));
			}

			/* Exponential Variogram Model? */

		} else if (it[ist] == 2) {
			cova += cc[ist] * exp(h * -3.f / aa[ist]);

			/* Gaussian Variogram Model? */

		} else if (it[ist] == 3) {
			/* Computing 2nd power */
			r__1 = h / aa[ist];
			cova += cc[ist] * exp(r__1 * r__1 * -3.f);

			/* Power Variogram Model? */

		} else if (it[ist] == 4) {
			d__1 = (double) h;
			d__2 = (double) aa[ist];
			cova += *cmax - cc[ist] * pow(d__1, d__2);

			/* Hole Effect Model? */

		} else if (it[ist] == 5) {
			/*                 d = 10.0 * aa(ist) */
			/*                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI) */
			cova += cc[ist] * cos(h / aa[ist] * PI );
		}
	}

	/* Finished: */

	return cova;
} /* cova3_ */

