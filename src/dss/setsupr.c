#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)






/* ----------------------------------------------------------------------- */
/*           Establish Super Block Search Limits and Sort Data */
/*           ************************************************* */

/* This subroutine sets up a 3-D "super block" model and orders the data */
/* by super block number.  The limits of the super block is set to the */
/* minimum and maximum limits of the grid; data outside are assigned to */
/* the nearest edge block. */

/* The idea is to establish a 3-D block network that contains all the */
/* relevant data.  The data are then sorted by their index location in */
/* the search network, i.e., the index location is given after knowing */
/* the block index in each coordinate direction (ix,iy,iz): */
/*          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix */
/* An array, the same size as the number of super blocks, is constructed */
/* that contains the cumulative number of data in the model.  With this */
/* array it is easy to quickly check what data are located near any given */
/* location. */

/* INPUT VARIABLES: */
/*   nx,xmn,xsiz      Definition of the X grid being considered */
/*   ny,ymn,ysiz      Definition of the Y grid being considered */
/*   nz,zmn,zsiz      Definition of the Z grid being considered */
/*   nd               Number of data */
/*   x(nd)            X coordinates of the data */
/*   y(nd)            Y coordinates of the data */
/*   z(nd)            Z coordinates of the data */
/*   vr(nd)           Variable at each location. */
/*   tmp(nd)          Temporary storage to keep track of the super block */
/*                      index associated to each data (uses the same */
/*                      storage already allocated for the simulation) */
/*   nsec             Number of secondary variables to carry with vr */
/*   sec1(nd)         First secondary variable (if nsec >= 1) */
/*   sec2(nd)         Second secondary variable (if nsec >= 2) */
/*   sec3(nd)         Third secondary variable (if nsec = 3) */
/*   MAXSB[X,Y,Z]     Maximum size of super block network */

/* OUTPUT VARIABLES: */
/*   nisb()                Array with cumulative number of data in each */
/*                           super block. */
/*   nxsup,xmnsup,xsizsup  Definition of the X super block grid */
/*   nysup,ymnsup,ysizsup  Definition of the Y super block grid */
/*   nzsup,zmnsup,zsizsup  Definition of the Z super block grid */

/* EXTERNAL REFERENCES: */
/*   sortem           Sorting routine to sort the data */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 * sortem
 * getindx
 */

/** CUBOS utilizados
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 */

int setsupr(int *nx, float *xmn, float *xsiz, int *ny, float *ymn,
		float *ysiz, int *nz, float *zmn, float *zsiz, int * nd,
		float *x, float *y, float *z__, float *vr, float *tmp, int *nsec,
		float *sec1, float *sec2, float *sec3, int *maxsbx, int *maxsby,
		int *maxsbz, int *nisb, int *nxsup, float *xmnsup,
		float * xsizsup, int *nysup, float *ymnsup, float *ysizsup,
		int *nzsup, float *zmnsup, float *zsizsup)
{
        /* Table of constant values */
        int one = 1;

	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, ii, ix, iy, iz, nsort;
	int inflag;

#ifdef PROFILE
	profile.setsupr++;
#endif

	/* Establish the number and size of the super blocks: */

	/* Parameter adjustments */
	--nisb;
	--sec3;
	--sec2;
	--sec1;
	--tmp;
	--vr;
	--z__;
	--y;
	--x;

	/* Function Body */
	*nxsup = MIN(*nx, *maxsbx); /* min(*nx,*maxsbx); */
	*nysup = MIN(*ny, *maxsby); /* min(*ny,*maxsby); */
	*nzsup = MIN(*nz, *maxsbz); /* min(*nz,*maxsbz); */
	*xsizsup = (float) (*nx) * *xsiz / (float) (*nxsup);
	*ysizsup = (float) (*ny) * *ysiz / (float) (*nysup);
	*zsizsup = (float) (*nz) * *zsiz / (float) (*nzsup);
	*xmnsup = *xmn - *xsiz * .5f + *xsizsup * .5f;
	*ymnsup = *ymn - *ysiz * .5f + *ysizsup * .5f;
	*zmnsup = *zmn - *zsiz * .5f + *zsizsup * .5f;

	/* Initialize the extra super block array to zeros: */

	i__1 = *nxsup * *nysup * *nzsup;
	for (i__ = 1; i__ <= i__1; ++i__) {
		nisb[i__] = 0;
	}

	/* Loop over all the data assigning the data to a super block and */
	/* accumulating how many data are in each super block: */

	i__1 = *nd;
	for (i__ = 1; i__ <= i__1; ++i__) {
	/*
		getindx(nxsup, xmnsup, xsizsup, &x[i__], &ix, &inflag);
		getindx(nysup, ymnsup, ysizsup, &y[i__], &iy, &inflag);
		getindx(nzsup, zmnsup, zsizsup, &z__[i__], &iz, &inflag);
	*/
		GETINDX(nxsup, xmnsup, xsizsup, &x[i__], &ix);
		GETINDX(nysup, ymnsup, ysizsup, &y[i__], &iy);
		GETINDX(nzsup, zmnsup, zsizsup, &z__[i__], &iz);
		ii = ix + (iy - 1) * *nxsup + (iz - 1) * *nxsup * *nysup;
		tmp[i__] = (float) ii;
		++nisb[ii];
	}

	/* Sort the data by ascending super block number: */

	nsort = *nsec + 4;
	sortem(&one, nd,
               &tmp[1],
               &nsort, &x[1],
               &y[1], &z__[1], &vr[1], &sec1[1], &sec2[1], &sec3[1]);

	/* Set up array nisb with the starting address of the block data: */

	i__1 = *nxsup * *nysup * *nzsup - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		nisb[i__ + 1] = nisb[i__] + nisb[i__ + 1];
	}

	/* Finished: */

	return 0;
} /* setsupr_ */

