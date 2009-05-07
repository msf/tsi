#include <math.h>
#include "dss.h"

static 
float my_roundf(const float val)
{
    return ((val - floorf(val)) < 0.5) ? floorf(val) : ceilf(val);
}


int getPos(int x, int y, int z, int xlen, int xylen)
{	
	return (x + (xlen * (y -1)) + (xylen * (z-1)));
}


void get3Dcoords(int ind, int xlen, int xylen, int *x, int *y, int *z)
{
	int iy, iz;
	int t;
	
	ind--;

	iz = ind / xylen;
	*z = iz +1;

	t = (ind % xylen);
	iy = t / xlen;
	*y = iy +1;

	*x = (ind % xlen) + 1;
}


/*     Gets the coordinate index location of a point within a grid 
 min     origin at the center of the first cell 
 siz     size of the cells 
 loc     location of the point being considered 

 coordinate index should be checked to see if its inside given grid by the caller.
 */

int getIndex(float min, float siz, float loc)
{
	return (int) my_roundf( 1 + (loc - min) / siz );
}


/* !Calculo do equivalente valor gaussiano        (SDSIM) */
float compute_gaussian_equiv(float cmean, unsigned size, float *vrtr, float *vrgtr)
{
    unsigned low, i, j;

    float vmy;
    if (cmean <= vrtr[0]) {
        return	vrgtr[0];
    }
    if (cmean >= vrtr[size - 1]) {
        return vrgtr[size - 1];
    }

    low = 0;
    i = size/2;
    // binary search for value closer to cmean in global histogram 
    do {
        if(vrtr[i] < cmean) {
            j = (i-low)/2;
            low = i;
            i += j;
        } else if(vrtr[i] > cmean) {
            i -= (i-low)/2;

        } else
            break;
    } while( low + 2 < i);

    vmy = vrgtr[i - 1] + (cmean - vrtr[i - 1]) * 
        (vrgtr[i] - vrgtr[i - 1]) / (vrtr[i] - vrtr[i - 1]);

    return vmy;
}

