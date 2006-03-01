#include <stdio.h>
#include "dss.h"

int getPos(int x, int y, int z, int xlen, int xylen)
{	
	return ((x-1) + (xlen * (y -1)) + (xylen * (z-1)));
}
