/* ----------------------------------------------------------------------- */
/*              Search Within Super Block Search Limits */
/*              *************************************** */

/* This subroutine searches through all the data that have been tagged in */
/* the super block subroutine.  The close data are passed back in the */
/* index array "close".  An octant search is allowed. */

/* INPUT VARIABLES: */
/*   xloc,yloc,zloc   location of point being estimated/simulated */
/*   radsqd           squared search radius */
/*   irot             index of the rotation matrix for searching */
/*   MAXROT           size of rotation matrix arrays */
/*   rotmat           rotation matrices */
/*   nsbtosr          Number of super blocks to search */
/*   ixsbtosr         X offsets for super blocks to search */
/*   iysbtosr         Y offsets for super blocks to search */
/*   izsbtosr         Z offsets for super blocks to search */
/*   noct             If >0 then data will be partitioned into octants */
/*   nd               Number of data */
/*   x(nd)            X coordinates of the data */
/*   y(nd)            Y coordinates of the data */
/*   z(nd)            Z coordinates of the data */
/*   tmp(nd)          Temporary storage to keep track of the squared */
/*                      distance associated with each data */
/*   nisb()           Array with cumulative number of data in each */
/*                      super block. */
/*   nxsup,xmnsup,xsizsup  Definition of the X super block grid */
/*   nysup,ymnsup,ysizsup  Definition of the X super block grid */
/*   nzsup,zmnsup,zsizsup  Definition of the X super block grid */

/* OUTPUT VARIABLES: */
/*   nclose           Number of close data */
/*   close()          Index of close data */
/*   infoct           Number of informed octants (only computes if */
/*                      performing an octant search) */

/* EXTERNAL REFERENCES: */
/*   sqdist           Computes anisotropic squared distance */
/*   sortem           Sorts multiple arrays in ascending order */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 * sortem
 * sqdist
 */

/** CUBOS utilizados
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 */
#include "profile.h"

extern double sqdist(float *, float *, float *, float *, float *, float *, int *, int *, double *);
extern int sortem(int *, int *, float *, int *, float *, float *, float *, float *, float *, float *, float *); 
extern int getindx(int *, float *, float *, float *, int *, int *);
/* Table of constant values */

static int one = 1;

int srchsupr(float *xloc, float *yloc, float *zloc, float * radsqd,
		int *irot, int *maxrot, double *rotmat, int * nsbtosr,
		int *ixsbtosr,  int *iysbtosr, int *izsbtosr, int *noct,
		int *nd, float *x, float *y, float *z__, float *tmp, int *nisb,
		int *nxsup, float *xmnsup, float *xsizsup, int * nysup,
		float *ymnsup, float *ysizsup,  int *nzsup, float *zmnsup,
		float *zsizsup, int *nclose, float *close, int *infoct)
{
	/* System generated locals */
	int rotmat_dim1, rotmat_offset, i__1, i__2;

	/* Local variables */
	float c__, d__, e, f, g, h__;
	int i__, j, na, ii, iq;
	float dx, dy, dz;
	int ix, iy, iz, nt;
	double hsqd;
	int isup, nums, inoct[8], ixsup, iysup, izsup;
	int inflag;

#ifdef PROFILE
	profile.srchsupr++;
#endif

	/* Determine the super block location of point being estimated: */

	/* Parameter adjustments */
	rotmat_dim1 = *maxrot;
	rotmat_offset = 1 + (rotmat_dim1 << 2);
	rotmat -= rotmat_offset;
	--ixsbtosr;
	--iysbtosr;
	--izsbtosr;
	--x;
	--y;
	--z__;
	--tmp;
	--nisb;
	--close;

	/* Function Body */
	getindx(nxsup, xmnsup, xsizsup, xloc, &ix, &inflag);
	getindx(nysup, ymnsup, ysizsup, yloc, &iy, &inflag);
	getindx(nzsup, zmnsup, zsizsup, zloc, &iz, &inflag);

	/* Loop over all the possible Super Blocks: */

	*nclose = 0;
	i__1 = *nsbtosr;
	for (isup = 1; isup <= i__1; ++isup) {

		/* Is this super block within the grid system: */

		ixsup = ix + ixsbtosr[isup];
		iysup = iy + iysbtosr[isup];
		izsup = iz + izsbtosr[isup];
		if (ixsup <= 0 || ixsup > *nxsup || iysup <= 0 || iysup > *nysup || 
				izsup <= 0 || izsup > *nzsup) {
			continue;
		}

		/* Figure out how many samples in this super block: */

		ii = ixsup + (iysup - 1) * *nxsup + (izsup - 1) * *nxsup * *nysup;
		if (ii == 1) {
			nums = nisb[ii];
			i__ = 0;
		} else {
			nums = nisb[ii] - nisb[ii - 1];
			i__ = nisb[ii - 1];
		}

		/* Loop over all the data in this super block: */

		i__2 = nums;
		for (ii = 1; ii <= i__2; ++ii) {
			++i__;

			/* Check squared distance: */

			hsqd = sqdist(xloc, yloc, zloc, &x[i__], &y[i__], &z__[i__], 
					irot, maxrot, &rotmat[rotmat_offset]);
			if ((float) hsqd > *radsqd) {
				continue;
			}

			/* Accept this sample: */

			++(*nclose);
			close[*nclose] = (float) i__;
			tmp[*nclose] = (float) hsqd;
		}
	}

	/* Sort the nearby samples by distance to point being estimated: */

	sortem(&one, nclose, &tmp[1], &one, &close[1], &c__, &d__, &e, &f, &g, 
			&h__);

	/* If we aren't doing an octant search then just return: */

	if (*noct <= 0) {
		return 0;
	}

	/* PARTITION THE DATA INTO OCTANTS: */

	for (i__ = 1; i__ <= 8; ++i__) {
		inoct[i__ - 1] = 0;
	}

	/* Now pick up the closest samples in each octant: */

	nt = *noct << 3;
	na = 0;
	i__1 = *nclose;
	for (j = 1; j <= i__1; ++j) {
		i__ = (int) close[j];
		h__ = tmp[j];
		dx = x[i__] - *xloc;
		dy = y[i__] - *yloc;
		dz = z__[i__] - *zloc;
		if (dz < 0.f) {
			iq = 8;
			if (dx <= 0.f && dy > 0.f) {
				iq = 5;
			}
			if (dx > 0.f && dy >= 0.f) {
				iq = 6;
			}
			if (dx < 0.f && dy <= 0.f) {
				iq = 7;
			}
		}
		iq = 4;
		if (dx <= 0.f && dy > 0.f) {
			iq = 1;
		}
		if (dx > 0.f && dy >= 0.f) {
			iq = 2;
		}
		if (dx < 0.f && dy <= 0.f) {
			iq = 3;
		}
		++inoct[iq - 1];

		/* Keep this sample if the maximum has not been exceeded: */

		if (inoct[iq - 1] <= *noct) {
			++na;
			close[na] = (float) i__;
			tmp[na] = h__;
			if (na == nt) {
				break;
			}
		}
	}

	/* End of data selection. Compute number of informed octants and return: */

	*nclose = na;
	*infoct = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
		if (inoct[i__ - 1] > 0) {
			++(*infoct);
		}
	}

	/* Finished: */

	return 0;
} /* srchsupr_ */

