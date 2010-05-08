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
		covtable_lookup_vars_t * covtable_lookup,
		search_node_t	*nodes)
{

	/* Local variables */
	unsigned int i, j, k;
	int il, ind;
	int nctx, ncty, nctz;

	unsigned int count = 0;


	/* Consider all the nearby nodes until enough have been found: */
	/* Parameter adjustments */
	--sim;

	/* Function Body */
	nctx = ix - covtable_lookup->nctx - 1;
	ncty = iy - covtable_lookup->ncty - 1;
	nctz = iz - covtable_lookup->nctz - 1;
	for (il = 1; il < covtable_lookup->nlooku; ++il) {
		if (count == covtable_lookup->nodmax) {
			break;
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
		if (sim[ind] > general->nosim_value) {

			nodes[count].index = il;
			nodes[count].x = getAbsolutePos(general->xmn, general->xsiz, i);
			nodes[count].y = getAbsolutePos(general->ymn, general->ysiz, j);
			nodes[count].z = getAbsolutePos(general->zmn, general->zsiz, k);
			nodes[count].value = sim[ind];
			count++;
		}
	}
	/* 	Return to calling program: */
	return count;
} /* srchnd_ */

