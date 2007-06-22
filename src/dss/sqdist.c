#include "dss.h"


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


double sqdist(float x1, float y1, float z1, 
		float x2, float y2, float z2, 
		int ind, double rotmat[5][3][3])
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	int i;
	double dx, dy, dz, cont;


	/* Compute component distance vectors and the squared distance: */

	/* Parameter adjustments */
	ind--;

	/* Function Body */
	dx = (double) (x1 - x2);
	dy = (double) (y1 - y2);
	dz = (double) (z1 - z2);

	ret_val = 0.f;

	if(ind < 0 || ind > 4) {
		fprintf(stderr,"sqdist() - ERROR, invalid rotation matrix: %d, only 0-4 (1-5) are valid\n",ind);
		return ret_val;
	}

	for (i = 0; i < 3; ++i) {
		cont =  rotmat[ind][i][0] * dx;
		cont += rotmat[ind][i][1] * dy;
		cont += rotmat[ind][i][2] * dz;
		
		ret_val += cont * cont;
	}
	return ret_val;
} /* sqdist_ */

