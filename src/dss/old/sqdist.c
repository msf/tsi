#include "dss.h"

#undef PROFILE

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* ----------------------------------------------------------------------- */

/*    Squared Anisotropic Distance Calculation Given Matrix Indicator */
/*    *************************************************************** */

/* This routine calculates the anisotropic distance between two points */
/*  given the coordinates of each point and a definition of the */
/*  anisotropy. */


/* INPUT VARIABLES: */

/*   x1,y1,z1         Coordinates of first point */
/*   x2,y2,z2         Coordinates of second point */
/*   ind              The rotation matrix to use */
/*   MAXROT           The maximum number of rotation matrices dimensioned */
/*   rotmat           The rotation matrices */



/* OUTPUT VARIABLES: */

/*   sqdist           The squared distance accounting for the anisotropy */
/*                      and the rotation of coordinates (if any). */


/* NO EXTERNAL REFERENCES */


/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 * -
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
 */


double sqdist(float *x1, float *y1, float *z1, float *x2, float *y2,
		float *z2, int *ind, int *maxrot, double *rotmat)
{
	/* System generated locals */
	int rotmat_dim1, rotmat_offset;
	double ret_val;

	/* Local variables */
	int i__;
	double dx, dy, dz, cont;

#ifdef PROFILE
	profile.sqdist++;
#endif

	/* Compute component distance vectors and the squared distance: */

	/* Parameter adjustments */
	rotmat_dim1 = *maxrot;
	rotmat_offset = 1 + (rotmat_dim1 << 2);
	rotmat -= rotmat_offset;

	/* Function Body */
	dx = (double) (*x1 - *x2);
	dy = (double) (*y1 - *y2);
	dz = (double) (*z1 - *z2);
	ret_val = 0.f;
	for (i__ = 1; i__ <= 3; ++i__) {
		cont = rotmat[*ind + (i__ + 3) * rotmat_dim1] * dx + rotmat[*ind + (
				i__ + 6) * rotmat_dim1] * dy + rotmat[*ind + (i__ + 9) * 
			rotmat_dim1] * dz;
		ret_val += cont * cont;
	}
	return ret_val;
} /* sqdist_ */

