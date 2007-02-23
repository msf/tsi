#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"



inline int getPos(int x, int y, int z, int xlen, int xylen)
{	
	return ((x-1) + (xlen * (y -1)) + (xylen * (z-1)));
}


/*     Gets the coordinate index location of a point within a grid 
 min     origin at the center of the first cell 
 siz     size of the cells 
 loc     location of the point being considered 

 coordinate index should be checked to see if its inside given grid by the caller.
 */
int getIndex(float min, float siz, float loc)
{
	return (int) ( (loc - min) / siz );
}
