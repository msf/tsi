#include "dss.h"
#include "dss_legacy.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


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
/*   rotmat           rotation matrices */
/*   nsbtosr          Number of super blocks to search */
/*   ixsbtosr         X offsets for super blocks to search */
/*   iysbtosr         Y offsets for super blocks to search */
/*   izsbtosr         Z offsets for super blocks to search */
/*   noct             If >0 then data will be partitioned into octants */
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
 * getIndex
 */

/** CUBOS utilizados
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 */

int srchsupr(float *xloc, float *yloc, float *zloc,
             float *radsqd,
		     int *irot, double rotmat[][][],
             int *nsbtosr, int *ixsbtosr, int *iysbtosr, int *izsbtosr,
             int *noct,
             float *x, float *y, float *z,
             float *tmp,
             int *nisb,
             int *nxsup, float *xmnsup, float *xsizsup,
             int *nysup, float *ymnsup, float *ysizsup,
             int *nzsup, float *zmnsup, float *zsizsup,
             int *nclose, float *close,
             int *infoct)
{
	/* Table of constant values */

	/* System generated locals */
	int rotmat_dim1, rotmat_offset, i1, i2;

	/* Local variables */
	int i, j, na, ii, iq, h;
	float dx, dy, dz;
	int ix, iy, iz, nt;
	double hsqd;
	int isup, nums, inoct[8], ixsup, iysup, izsup;


	/* Determine the super block location of point being estimated: */

	/* Parameter adjustments */
	--ixsbtosr;
	--iysbtosr;
	--izsbtosr;
	--x;
	--y;
	--z;
	--tmp;
	--nisb;
	--close;

	/* Function Body */
	ix = getIndex(*xmnsup, *xsizsup, *xloc);
	iy = getIndex(*ymnsup, *ysizsup, *yloc);
	iz = getIndex(*zmnsup, *zsizsup, *zloc);

	/* Loop over all the possible Super Blocks: */

	*nclose = 0;
	i1 = *nsbtosr;
	for (isup = 1; isup <= i1; ++isup) {

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
			i = 0;
		} else {
			nums = nisb[ii] - nisb[ii - 1];
			i = nisb[ii - 1];
		}

		/* Loop over all the data in this super block: */

		i2 = nums;
		for (ii = 1; ii <= i2; ++ii) {
			++i;

			/* Check squared distance: */

			/* FIXME: ?loc nao devem ser pointrs */
			hsqd = sqdist(*xloc, *yloc, *zloc, x[i], y[i], z[i], *irot, rotmat);
			if ((float) hsqd > *radsqd) {
				continue;
			}

			/* Accept this sample: */

			++(*nclose);
			close[*nclose] = (float) i;
			tmp[*nclose] = (float) hsqd;
		}
	}

	/* Sort the nearby samples by distance to point being estimated: */
	sort_permute_float(1, *nclose, &tmp[1], &close[1]);
	/* sortem(&one, nclose, &tmp[1], &one, &close[1], &c__, &d__, &e, &f, &g, &h);
	 */

	/* If we aren't doing an octant search then just return: */

	if (*noct <= 0) {
		return 0;
	}

	/* PARTITION THE DATA INTO OCTANTS: */

	for (i = 1; i <= 8; ++i) {
		inoct[i - 1] = 0;
	}

	/* Now pick up the closest samples in each octant: */

	/* nt = *noct << 3;  //LPL */
	nt = *noct * 8;  /* LPL */
	na = 0;
	i1 = *nclose;
	for (j = 1; j <= i1; ++j) {
		i = (int) close[j];
		h = tmp[j];
		dx = x[i] - *xloc;
		dy = y[i] - *yloc;
		dz = z[i] - *zloc;
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
			close[na] = (float) i;
			tmp[na] = h;
			if (na == nt) {
				break;
			}
		}
	}

	/* End of data selection. Compute number of informed octants and return: */

	*nclose = na;
	*infoct = 0;
	for (i = 1; i <= 8; ++i) {
		if (inoct[i - 1] > 0) {
			++(*infoct);
		}
	}

	/* Finished: */

	return 0;
} /* srchsupr_ */

