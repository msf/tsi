#include <math.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* ----------------------------------------------------------------------- */
/*              Sets up an Anisotropic Rotation Matrix */
/*              ************************************** */
/* Sets up the matrix to transform cartesian coordinates to coordinates */
/* accounting for angles and anisotropy (see manual for a detailed */
/* definition): */

/* INPUT PARAMETERS: */

/*   ang1             Azimuth angle for principal direction */
/*   ang2             Dip angle for principal direction */
/*   ang3             Third rotation angle */
/*   anis1            First anisotropy ratio */
/*   anis2            Second anisotropy ratio */
/*   ind              matrix indicator to initialize */
/*   MAXROT           maximum number of rotation matrices dimensioned */
/*   rotmat           rotation matrices */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 */

/** CUBOS utilizados
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 */


int setrot(float *ang1, float *ang2, float *ang3, float *anis1, float *anis2,
		int *ind, int *maxrot, double *rotmat)
{
	/* System generated locals */
	int rotmat_dim1, rotmat_offset;

	/* Local variables */
	float beta;
	double cosa, cosb, sina, sinb, cost, afac1, afac2, sint;
	float alpha, theta;

#ifdef PROFILE
	profile.setrot++;
	profBegin("setrot");
#endif

	/* Converts the input angles to three angles which make more */
	/*  mathematical sense: */

	/*         alpha   angle between the major axis of anisotropy and the */
	/*                 E-W axis. Note: Counter clockwise is positive. */
	/*         beta    angle between major axis and the horizontal plane. */
	/*                 (The dip of the ellipsoid measured positive down) */
	/*         theta   Angle of rotation of minor axis about the major axis */
	/*                 of the ellipsoid. */

	/*      if(ang1.ge.0.0.and.ang1.lt.270.0) then */
	/*            alpha = (90.0   - ang1) * DEG2RAD */
	/*      else */
	/*            alpha = (450.0  - ang1) * DEG2RAD */
	/*      endif */
	/*      beta  = -1.0 * ang2 * DEG2RAD */

	/* Parameter adjustments */
	rotmat_dim1 = *maxrot;
	rotmat_offset = 1 + (rotmat_dim1 << 2);
	rotmat -= rotmat_offset;

	/* Function Body */
	alpha = *ang1 * .017453292522222223f;
	beta = *ang2 * .017453292522222223f;
	theta = *ang3 * .017453292522222223f;

	/* Get the required sines and cosines: */

	sina = (double) sin(alpha);
	sinb = (double) sin(beta);
	sint = (double) sin(theta);
	cosa = (double) cos(alpha);
	cosb = (double) cos(beta);
	cost = (double) cos(theta);

	/* Construct the rotation matrix in the required memory: */

	afac1 = 1.f / (double) MAX(*anis1, 1e-20f); /* dmax(*anis1,1e-20f); */
	afac2 = 1.f / (double) MAX(*anis2, 1e-20f); /* dmax(*anis2,1e-20f); */
	rotmat[*ind + (rotmat_dim1 << 2)] = cosb * cosa;
	rotmat[*ind + rotmat_dim1 * 7] = cosb * sina;
	rotmat[*ind + rotmat_dim1 * 10] = -sinb;
	rotmat[*ind + rotmat_dim1 * 5] = afac1 * (-cost * sina + sint * sinb * 
			cosa);
	rotmat[*ind + (rotmat_dim1 << 3)] = afac1 * (cost * cosa + sint * sinb * 
			sina);
	rotmat[*ind + rotmat_dim1 * 11] = afac1 * (sint * cosb);
	rotmat[*ind + rotmat_dim1 * 6] = afac2 * (sint * sina + cost * sinb * 
			cosa);
	rotmat[*ind + rotmat_dim1 * 9] = afac2 * (-sint * cosa + cost * sinb * 
			sina);
	rotmat[*ind + rotmat_dim1 * 12] = afac2 * (cost * cosb);

#ifdef PROFILE
	profEnd("setrot");
#endif

	/* Return to calling program: */
	return 0;
} /* setrot_ */


