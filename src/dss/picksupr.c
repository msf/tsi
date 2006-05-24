#include "dss.h"
#include "dss_legacy.h"

#undef PROFILE

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


/* ----------------------------------------------------------------------- */
/*             Establish Which Super Blocks to Search */
/*             ************************************** */
/* This subroutine establishes which super blocks must be searched given */
/* that a point being estimated/simulated falls within a super block */
/* centered at 0,0,0. */

/* INPUT VARIABLES: */
/*   nxsup,xsizsup    Definition of the X super block grid */
/*   nysup,ysizsup    Definition of the Y super block grid */
/*   nzsup,zsizsup    Definition of the Z super block grid */
/*   irot             index of the rotation matrix for searching */
/*   rotmat           rotation matrices */
/*   radsqd           squared search radius */

/* OUTPUT VARIABLES: */
/*   nsbtosr          Number of super blocks to search */
/*   ixsbtosr         X offsets for super blocks to search */
/*   iysbtosr         Y offsets for super blocks to search */
/*   izsbtosr         Z offsets for super blocks to search */

/* EXTERNAL REFERENCES: */
/*   sqdist           Computes anisotropic squared distance */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 * sqdist
 */

/** CUBOS utilizados
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 */


int picksup(int nxsup, float xsizsup, int nysup, float ysizsup,
		int nzsup, float zsizsup, int irot,
		double rotmat[5][3][3], float radsqd, int *nsbtosr, int * ixsbtosr,
		int *iysbtosr, int *izsbtosr)
{
        /* Table of constant values */


	/* System generated locals */
	int rotmat_dim1, rotmat_offset, t1, t2, t3;

	/* Local variables */
	double shortest;
	int t, j, k, i1, j1, k1, i2, j2, k2;
	float xo, yo, zo;
	double hsqd;
	float xdis, ydis, zdis;


	/* MAIN Loop over all possible super blocks: */

	/* Parameter adjustments */
	--ixsbtosr;
	--iysbtosr;
	--izsbtosr;

	/* Function Body */
	*nsbtosr = 0;
	t1 = nxsup - 1;
	for (t = -(nxsup - 1); t <= t1; ++t) {
		t2 = nysup - 1;
		for (j = -(nysup - 1); j <= t2; ++j) {
			t3 = nzsup - 1;
			for (k = -(nzsup - 1); k <= t3; ++k) {
				xo = (float) t * xsizsup;
				yo = (float) j * ysizsup;
				zo = (float) k * zsizsup;

				/* Find the closest distance between the corners of the super blocks: */

				shortest = 1e21f;
				for (i1 = -1; i1 <= 1; ++i1) {
					for (j1 = -1; j1 <= 1; ++j1) {
						for (k1 = -1; k1 <= 1; ++k1) {
							for (i2 = -1; i2 <= 1; ++i2) {
								for (j2 = -1; j2 <= 1; ++j2) {
									for (k2 = -1; k2 <= 1; ++k2) {
										if (i1 != 0 && j1 != 0 && k1 != 0 && 
												i2 != 0 && j2 != 0 && k2 != 0)
										{
											xdis = (float) (i1 - i2) * .5f * xsizsup + xo;
											ydis = (float) (j1 - j2) * .5f * ysizsup + yo;
											zdis = (float) (k1 - k2) * .5f * zsizsup + zo;
											hsqd = sqdist(0, 0, 0, xdis, ydis, zdis, irot, rotmat);
											if (hsqd < shortest) {
												shortest = hsqd;
											}
										}
									}
								}
							}
						}
					}
				}

				/* Keep this super block if it is close enoutgh: */

				if ((float) shortest <= radsqd) {
					++(*nsbtosr);
					ixsbtosr[*nsbtosr] = t;
					iysbtosr[*nsbtosr] = j;
					izsbtosr[*nsbtosr] = k;
				}
			}
		}
	}

	/* Finished: */

	return 0;
} /* picksup_ */

