#include <stdio.h>
#include <string.h>
#include "dss.h"
#include "memdebug.h"

int getIndex(float min, float siz, float loc)
{
	return (int) ( ((loc - min) / siz) + .5);
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

	int x, y, z;
	float f;
	double origVal;
	unsigned int i,j,t;
	unsigned int *temp;
	double *vals;

	t = wellsDataSize / general->nvari;
	
	temp = (unsigned int *) tsi_malloc(sizeof(unsigned int) * t);
	vals = (double *) tsi_malloc(sizeof(double) * t);		

	if( wellsDataSize % general->nvari ) {
		printf("readWellsData: ERROR: incorrect number of values\n");
		return -1;
	}

	i = 0;
	j = 0;
	while( i < wellsDataSize) {
		f = (float) wellsData[i++];
		x = getIndex(general->xmn, general->xsiz, f); 
		if( x < 0 || x >= general->nx) {
			i += 3;
			continue;
		}
		f = (float) wellsData[i++];
		y = getIndex(general->ymn, general->ysiz, f);
		if( y < 0 || y >= general->ny) {
			i += 2;
			continue;
		}
		f = (float) wellsData[i++];
		z = getIndex(general->zmn, general->zsiz, f);
		if( z < 0 || z >= general->nz) {
			i++;
			continue;
		}
		origVal = wellsData[i++];
		
		t = (unsigned int) getPos(x,y,z,general->nx,general->nxy);
		if(t >= (unsigned int) general->nxyz) {
			continue;
		}
		temp[j] = t;
		vals[j] = origVal;
		j++;
	}

	/* j is the size of the WellsData */
	printf("readWellsData: %d points\n",j);
	general->wellsNPoints = j;
	general->wellsDataVal = (double *) tsi_malloc(sizeof(double) * j);
	general->wellsDataPos = (unsigned int *) tsi_malloc(sizeof(unsigned int) * j); 
	memcpy(general->wellsDataVal, vals, sizeof(double) * j);
	memcpy(general->wellsDataPos, temp, sizeof(unsigned int) * j);
	tsi_free(vals);
	tsi_free(temp);
			
	return 0;
}
