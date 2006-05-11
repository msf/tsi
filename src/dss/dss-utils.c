#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"



inline int getPos(int x, int y, int z, int xlen, int xylen)
{	
	return ((x-1) + (xlen * (y -1)) + (xylen * (z-1)));
}


int getIndex(float min, float siz, float loc)
{
	return (int) ( ((loc - min) / siz) + .5);
}
