#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* ----------------------------------------------------------------------- */
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

int getindx(int *n, float *min__, float *siz, float *loc, int *index, int *inflag)
{

	/* Compute the index of "loc": */

	*index = (int) ((*loc - *min__) / *siz + 1.5f);

	/* Check to see if in or out: */

	if (*index < 1) {
		*index = 1;
		*inflag = FALSE;
	} else if (*index > *n) {
		*index = *n;
		*inflag = FALSE;
	} else {
		*inflag = TRUE;
	}

	/* Return to calling program: */

	return 0;
} /* getindx_ */



