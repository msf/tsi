#include <stdio.h>
#include <string.h>
#include "dss.h"
#include "memdebug.h"

double getPosd(double x, double y, double z, int xlen, int xylen)
{
	        return ((x-1) + (xlen * (y -1)) + (xylen * (z-1)));
}

/* wellsData may have points of outsize the current grid size
 * so, we must know how many are inside the grid,
 * alocate arrays for those, and copy them!
 *
 * all this should be done outside the dss, since its constant throughout the whole execution
 * 
 */
int readWellsData(general_vars_t * general, double * wellsData, unsigned int wellsDataSize) 
{

	double x, y, z;
	double origVal;
	unsigned int i,j,t;
	unsigned int *temp;
	double *vals;

	t = wellsDataSize / general->nvari;
	
	temp = (unsigned int *) my_malloc(sizeof(unsigned int) * t);
	vals = (double *) my_malloc(sizeof(double) * t);		

	if( wellsDataSize % general->nvari ) {
		printf("readWellsData: ERROR: incorrect number of values\n");
		return -1;
	}

	i = 0;
	j = 0;
	while( i < wellsDataSize) {
		x = wellsData[i++];
		y = wellsData[i++];
		z = wellsData[i++];
		origVal = wellsData[i++];
		if( origVal < general->tmin ||
			origVal > general->tmax) {
			continue;
		}
		
		t = (unsigned int) getPosd(x,y,y,general->nx,general->nxy);
		if(t >= (unsigned int) general->nxyz) {
			continue;
		}
		temp[j] = t;
		vals[j] = origVal;
		j++;
	}

	/* j is the size of the WellsData */
	general->wellsNPoints = j;
	general->wellsDataVal = (double *) my_malloc(sizeof(double) * j);
	general->wellsDataPos = (unsigned int *) my_malloc(sizeof(unsigned int) * j); 
	memcpy(general->wellsDataVal, vals, sizeof(double) * j);
	memcpy(general->wellsDataPos, temp, sizeof(unsigned int) * j);
	my_free(vals);
	my_free(temp);
			
	return 0;
}
