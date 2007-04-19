#include <math.h>
#include "dss.h"

float my_roundf(const float val)
{
    return ((val - floorf(val)) < 0.5) ? floorf(val) : ceilf(val);
}

int getPos(int x, int y, int z, int xlen, int xylen)
{	
	return (x + (xlen * (y -1)) + (xylen * (z-1)));
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
