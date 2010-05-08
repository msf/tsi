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
int setrot(float ang1, float ang2, float ang3, float anis1, float anis2,
		int ind, double rotmat[5][3][3])
{
	/* System generated locals */

	/* Local variables */
	float beta;
	double cosa, cosb, sina, sinb, cost, afac1, afac2, sint;
	float alpha, theta;


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


	/* Function Body */
	alpha = ang1 * .017453292522222223f;
	beta = ang2 * .017453292522222223f;
	theta = ang3 * .017453292522222223f;

	/* Get the required sines and cosines: */

	sina = (double) sin(alpha);
	sinb = (double) sin(beta);
	sint = (double) sin(theta);
	cosa = (double) cos(alpha);
	cosb = (double) cos(beta);
	cost = (double) cos(theta);

	/* Construct the rotation matrix in the required memory: */

	afac1 = 1.f / (double) MAX(anis1, 1e-20f);
	afac2 = 1.f / (double) MAX(anis2, 1e-20f);


	if(ind < 0 || ind > 4) {
		fprintf(stderr,"setrot() - ERROR, invalid rotation matrix: %d, only 0-4 are valid\n",ind);
		return 1;
	}

	rotmat[ind][0][0] = cosb * cosa;
	rotmat[ind][0][1] = cosb * sina;
	rotmat[ind][0][2] = -sinb;
	rotmat[ind][1][0] = afac1 * (-cost * sina + sint * sinb * cosa);
	rotmat[ind][1][1] = afac1 * (cost * cosa + sint * sinb * sina);
	rotmat[ind][1][2] = afac1 * (sint * cosb);
	rotmat[ind][2][0] = afac2 * (sint * sina + cost * sinb * cosa);
	rotmat[ind][2][1] = afac2 * (-sint * cosa + cost * sinb * sina);
	rotmat[ind][2][2] = afac2 * (cost * cosb);


	/* Return to calling program: */
	return 0;
} /* setrot_ */


