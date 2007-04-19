#include "dss.h"
#include "dss_legacy.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* ----------------------------------------------------------------------- */
/*               Search for nearby Simulated Grid nodes */

/* The idea is to spiral away from the node being simulated and note all */
/* the nearby nodes that have been simulated. */

/* INPUT VARIABLES: */
/*   *ix,*iy,iz        index of the point currently being simulated */
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



int srchnod(int ix, int iy, int iz, float *sim,
		general_vars_t * general,
		search_vars_t * search,
		covtable_lookup_vars_t * covtable_lookup)
{

	/* Local variables */
	int i, j, k, il, iq, ind, idx, idy, idz, ninoct[8];
	int nctx, ncty, nctz;


	/* Consider all the nearby nodes until enough have been found: */
	/* Parameter adjustments */
	--sim;

	/* Function Body */
	covtable_lookup->ncnode = 0;
	if (search->noct > 0) {
		for(i= 0; i < 8; i++)
			ninoct[i] = 0;
	}
	
	nctx = ix - covtable_lookup->nctx - 1;
	ncty = iy - covtable_lookup->ncty - 1;
	nctz = iz - covtable_lookup->nctz - 1;
	for (il = 1; il < covtable_lookup->nlooku; ++il) {
		if (covtable_lookup->ncnode == covtable_lookup->nodmax) {
			return 0;
		}
		i = nctx + covtable_lookup->ixnode[il];
		if( i < 1 || i > general->nx)
			continue;
		
		j = ncty + covtable_lookup->iynode[il];
		if( j < 1 || j > general->ny )
			continue;

		k = nctz + covtable_lookup->iznode[il];
		if( k < 1 || k > general->nz )
			continue;

		ind = getPos(i, j, k, general->nx, general->nxy);
		if (sim[ind] > general->nosvalue) {
			/* Check the number of data already taken from this octant: */
			if (search->noct > 0) {
				idx = ix - i;
				idy = iy - j;
				idz = iz - k;
				if (idz > 0) {
					iq = 3;
					if (idx <= 0 && idy > 0) {
						iq = 0;
					}
					else if (idx > 0 && idy >= 0) {
						iq = 1;
					}
					else if (idx < 0 && idy <= 0) {
						iq = 2;
					}
				} else {
					iq = 7;
					if (idx <= 0 && idy > 0) {
						iq = 4;
					}
					else if (idx > 0 && idy >= 0) {
						iq = 5;
					}
					else if (idx < 0 && idy <= 0) {
						iq = 6;
					}
				}
				++ninoct[iq];
				if (ninoct[iq] > search->noct) {
					continue;
				}
			}
			covtable_lookup->icnode[covtable_lookup->ncnode] = il +1;
			covtable_lookup->cnodex[covtable_lookup->ncnode] = general->xmn + (float) (i - 1) * general->xsiz;
			covtable_lookup->cnodey[covtable_lookup->ncnode] = general->ymn + (float) (j - 1) * general->ysiz;
			covtable_lookup->cnodez[covtable_lookup->ncnode] = general->zmn + (float) (k - 1) * general->zsiz;
			covtable_lookup->cnodev[covtable_lookup->ncnode] = sim[ind];
			++covtable_lookup->ncnode;
		}
	}
	/* 	Return to calling program: */
	return 0;
} /* srchnd_ */

