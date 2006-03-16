#include <stdlib.h>
#include "dss.h"
#include "dss_legacy.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)

#define ACORNI_SIZE 17

/* ----------------------------------------------------------------------- 
 * Fortran implementation of ACORN random number generator of order less 
 * than or equal to 12 (higher orders can be obtained by increasing the 
 * parameter value MAXORD). 
 *
 * NOTES: 1. The variable idum is a dummy variable. The common block 
 *           IACO is used to transfer data into the function. 
 *        2. Before the first call to ACORN the common block IACO must 
 *           be initialised by the user, as follows. The values of 
 *           variables in the common block must not subsequently be 
 *           changed by the user. 
 *             KORDEI - order of generator required ( must be =< MAXORD) 
 *             MAXINT - modulus for generator, must be chosen small 
 *                      enough that 2*MAXINT does not overflow 
 *             ixv(1) - seed for random number generator 
 *                      require 0 < ixv(1) < MAXINT 
 *             (ixv(I+1),I=1,KORDEI) 
 *                    - KORDEI initial values for generator 
 *                      require 0 =< ixv(I+1) < MAXINT 
 *        3. After initialisation, each call to ACORN generates a single 
 *           random number between 0 and 1.
 *        4. An example of suitable values for parameters is
 *             KORDEI   = 10 
 *             MAXINT   = 2**30 
 *             ixv(1)   = an odd integer in the (approximate) range 
 *                        (0.001 * MAXINT) to (0.999 * MAXINT) 
 *             ixv(I+1) = 0, I=1,KORDEI 
 *
 * Author: R.S.Wikramaratna,                           Date: October 1990 
 * ----------------------------------------------------------------------- 
 */
int * new_acorni(int seed)
{
	int i, *ixv;

	/* we must guarantee that seed*2 < MAX_INT */
	--seed;
	seed /= 2;
        ixv = (int *) malloc(ACORNI_SIZE*sizeof(int));
	for(i = 1; i < ACORNI_SIZE; i++) ixv[i] = 0;
	ixv[0] = seed;

        for (i = 0; i < 1000; i++) acorni(ixv);
	return ixv;
} /* new_acorni */



double acorni(int *ixv)
{
	int i;

	for (i = 1; i < ACORNI_SIZE; ++i) {
		ixv[i] += ixv[i - 1];
		/* LPL: old code */
			if (ixv[i] >= 1073741824) {
			ixv[i] += -1073741824;
			}
		/* LPL: hack only for 32 bit ints */
		/* ixv[i] &= 1073741823; */
	}

	return (double) ((double)ixv[ACORNI_SIZE-1] / 1073741824.f);
} /* acorni */

