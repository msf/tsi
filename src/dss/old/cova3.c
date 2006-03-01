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

#include <math.h>

#include "profile.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))

extern double sqdist(float *, float *, float *, float *, float *, float *, int *, int *, double *);

int cova3(float *x1, float *y1, float *z1, float *x2, float *y2, float *z2,
		int *ivarg, int *nst, int *maxnst, float *c0, int *it, float *cc,
		float *aa, int *irot, int *maxrot, double *rotmat, float *cmax,
		float *cova)
/* Calculate the maximum covariance value (used for zero distances and */
/* for power model covariance): */
{
	/* System generated locals */
	int rotmat_dim1, rotmat_offset, i__1, i__2;
	float r__1;
	double d__1, d__2;


	/* Local variables */
	float h__, hr;
	int ir, is, ist;
	double hsqd;
	int istart;

#ifdef PROFILE
	profile.cova3++;
#endif

	/* Parameter adjustments */
	--nst;
	--c0;
	--it;
	--cc;
	--aa;
	rotmat_dim1 = *maxrot;
	rotmat_offset = 1 + (rotmat_dim1 << 2);
	rotmat -= rotmat_offset;

	/* Function Body */
	istart = (*ivarg - 1) * *maxnst + 1;
	*cmax = c0[*ivarg];
	i__1 = nst[*ivarg];
	for (is = 1; is <= i__1; ++is) {
		ist = istart + is - 1;
		if (it[ist] == 4) {
			*cmax += 999.f;
		} else {
			*cmax += cc[ist];
		}
	}

	/* Check for "zero" distance, return with cmax if so: */

	hsqd = sqdist(x1, y1, z1, x2, y2, z2, irot, maxrot, &rotmat[
			rotmat_offset]);
	if ((float) hsqd < 1e-10f) {
		*cova = *cmax;
		return 0;
	}

	/* Loop over all the structures: */

	*cova = 0.f;
	i__1 = nst[*ivarg];
	for (is = 1; is <= i__1; ++is) {
		ist = istart + is - 1;

		/* Compute the appropriate distance: */

		if (ist != 1) {
			/* Computing MIN */
			i__2 = *irot + is - 1;
			ir = MIN(i__2, *maxrot); //min(i__2,*maxrot);
			hsqd = sqdist(x1, y1, z1, x2, y2, z2, &ir, maxrot, &rotmat[
					rotmat_offset]);
		}
		h__ = (float) sqrt(hsqd);

		/* Spherical Variogram Model? */

		if (it[ist] == 1) {
			hr = h__ / aa[ist];
			if (hr < 1.f) {
				*cova += cc[ist] * (1.f - hr * (1.5f - hr * .5f * hr));
			}

			/* Exponential Variogram Model? */

		} else if (it[ist] == 2) {
			*cova += cc[ist] * exp(h__ * -3.f / aa[ist]);

			/* Gaussian Variogram Model? */

		} else if (it[ist] == 3) {
			/* Computing 2nd power */
			r__1 = h__ / aa[ist];
			*cova += cc[ist] * exp(r__1 * r__1 * -3.f);

			/* Power Variogram Model? */

		} else if (it[ist] == 4) {
			d__1 = (double) h__;
			d__2 = (double) aa[ist];
			*cova = *cova + *cmax - cc[ist] * pow(d__1, d__2);

			/* Hole Effect Model? */

		} else if (it[ist] == 5) {
			/*                 d = 10.0 * aa(ist) */
			/*                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI) */
			*cova += cc[ist] * cos(h__ / aa[ist] * 3.14159265f);
		}
	}

	/* Finished: */

	return 0;
} /* cova3_ */

