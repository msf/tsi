#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* ----------------------------------------------------------------------- */
/*               Search for nearby Simulated Grid nodes */

/* The idea is to spiral away from the node being simulated and note all */
/* the nearby nodes that have been simulated. */

/* INPUT VARIABLES: */
/*   ix,iy,iz        index of the point currently being simulated */
/*   sim             the realization so far */
/*   nodmax          the maximum number of nodes that we want */
/*   nlooku          the number of nodes in the look up table */
/*   i[x,y,z]node    the relative indices of those nodes. */
/*   [x,y,z]mn       the origin of the global grid netwrok */
/*   [x,y,z]siz      the spacing of the grid nodes. */

/* OUTPUT VARIABLES: */
/*   ncnode          the number of close nodes */
/*   icnode()        the number in the look up table */
/*   cnode[x,y,z]()  the location of the nodes */
/*   cnodev()        the values at the nodes */
/* ----------------------------------------------------------------------- */


/**
 * funcoes utilizadas
 * -
 */

/** CUBOS utilizados
 * sim
 * -
 */ 

/** CUBOS _nao_ utilizados
 * tmp
 * order
 * clc
 * lvm
 * e os cov*
 */

/** structs globais utilizadas:
 * general
 * search
 * covtable_lookup
 */



int srchnod(int *ix, int *iy, int *iz, float *sim,
		general_vars_t * general,
		search_vars_t * search,
		covtable_lookup_vars_t * covtable_lookup)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, j, k, il, iq, ind, idx, idy, idz, ninoct[8];

#ifdef PROFILE
	profile.srchnd++;
#endif

	/* Consider all the nearby nodes until enough have been found: */
	/* Parameter adjustments */
	--sim;

	/* Function Body */
	covtable_lookup->ncnode = 0;
	if (search->noct > 0) {
		for (i__ = 1; i__ <= 8; ++i__) {
			ninoct[i__ - 1] = 0;
		}
	}
	i__1 = covtable_lookup->nlooku;
	for (il = 2; il <= i__1; ++il) {
		if (covtable_lookup->ncnode == covtable_lookup->nodmax) {
			return 0;
		}
		i__ = *ix + (covtable_lookup->ixnode[il - 1] - covtable_lookup->nctx - 1);
		j = *iy + (covtable_lookup->iynode[il - 1] - covtable_lookup->ncty - 1);
		k = *iz + (covtable_lookup->iznode[il - 1] - covtable_lookup->nctz - 1);
		if (i__ < 1 || j < 1 || k < 1) {
			continue;
		}
		if (i__ > general->nx || j > general->ny || k > general->nz) {
			continue;
		}
		ind = i__ + (j - 1) * general->nx + (k - 1) * general->nxy;
		if (sim[ind] > general->nosvalue) {
			/* Check the number of data already taken from this octant: */
			if (search->noct > 0) {
				idx = *ix - i__;
				idy = *iy - j;
				idz = *iz - k;
				if (idz > 0) {
					iq = 4;
					if (idx <= 0 && idy > 0) {
						iq = 1;
					}
					if (idx > 0 && idy >= 0) {
						iq = 2;
					}
					if (idx < 0 && idy <= 0) {
						iq = 3;
					}
				} else {
					iq = 8;
					if (idx <= 0 && idy > 0) {
						iq = 5;
					}
					if (idx > 0 && idy >= 0) {
						iq = 6;
					}
					if (idx < 0 && idy <= 0) {
						iq = 7;
					}
				}
				++ninoct[iq - 1];
				if (ninoct[iq - 1] > search->noct) {
					continue;
				}
			}
			++covtable_lookup->ncnode;
			covtable_lookup->icnode[covtable_lookup->ncnode - 1] = il;
			covtable_lookup->cnodex[covtable_lookup->ncnode - 1] = general->xmn + (float) (i__ 
					- 1) * general->xsiz;
			covtable_lookup->cnodey[covtable_lookup->ncnode - 1] = general->ymn + (float) (j - 
					1) * general->ysiz;
			covtable_lookup->cnodez[covtable_lookup->ncnode - 1] = general->zmn + (float) (k - 
					1) * general->zsiz;
			covtable_lookup->cnodev[covtable_lookup->ncnode - 1] = sim[ind];
		}
	}
	/* 	Return to calling program: */
	return 0;
} /* srchnd_ */

