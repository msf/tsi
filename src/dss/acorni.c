#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


int ixv[17];    /* acorni */

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
void newAcorni(int seed)
{
	int i;

	/* we must guarantee that seed*2 < MAX_INT */
	--seed;
	seed /= 2; 


	for(i = 1; i <= 16; ++i) ixv[i] = 0;
	ixv[0] = seed;
} /* newAcorni_ */



double acorni()
{
	int i;


	for (i = 1; i <= 16; ++i) {
		ixv[i] += ixv[i - 1];
		/* LPL: old code */
			if (ixv[i] >= 1073741824) {
			ixv[i] += -1073741824;
			}
		/* LPL: hack only for 32 bit ints */
		/* ixv[i] &= 1073741823; */
	}

	return (double) ((double)ixv[16] / 1073741824.f);
} /* acorni_ */

