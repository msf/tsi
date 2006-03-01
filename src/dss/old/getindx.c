/* ----------------------------------------------------------------------- */
/*     Gets the coordinate index location of a point within a grid */
/*     *********************************************************** */

/* n       number of "nodes" or "cells" in this coordinate direction */
/* min     origin at the center of the first cell */
/* siz     size of the cells */
/* loc     location of the point being considered */
/* index   output index within [1,n] */
/* inflag  true if the location is actually in the grid (false otherwise */
/*         e.g., if the location is outside then index will be set to */
/*         nearest boundary */
/* ----------------------------------------------------------------------- */

/** Funcoes utilizadas
 */

/** CUBOS utilizados
 * 
 */ 

/** CUBOS _nao_ utilizados
 * -
 */

/** structs globais utilizadas:
 */
#include "profile.h"

#define TRUE_ (1)
#define FALSE_ (0)

int getindx(int *n, float *min__, float *siz, float *loc, int *index, int *inflag)
{
#ifdef PROFILE
	profile.getindx++;
#endif

	/* Compute the index of "loc": */

	*index = (int) ((*loc - *min__) / *siz + 1.5f);

	/* Check to see if in or out: */

	if (*index < 1) {
		*index = 1;
		*inflag = FALSE_;
	} else if (*index > *n) {
		*index = *n;
		*inflag = FALSE_;
	} else {
		*inflag = TRUE_;
	}

	/* Return to calling program: */

	return 0;
} /* getindx_ */

