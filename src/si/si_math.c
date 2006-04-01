/* si_math.c */

#include <math.h>
#include <stdlib.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "si.h"


/* wavelet functions */
float point_value(si *s, int point);
float index_value(si *s, int index);


unsigned int getPoint(si *s, int x, int y, int z)
{

//	return (x-1) + ( (y-1) * s->xsize) + ( (z-1) * s->xsize * s->ysize);
	return ((z)* s->xsize * s->ysize + (y) * s->xsize + (x+1)) - 1;
}

int make_reflections_grid(si *s, float *AI, float *RG) 
{

	unsigned int x, y, z;
	double value;

	for(x = 0; x < s->xsize; x++)
		for(y = 0; y < s->ysize; y++)
			for(z = 0; z < s->zsize -1; z++) {
				 value = (AI[getPoint(s,x,y,z+1)]-AI[getPoint(s,x,y,z)]) /
					     (AI[getPoint(s,x,y,z+1)]+AI[getPoint(s,x,y,z)]);
	             RG[getPoint(s,x,y,z)] = value; 
			}
	for(x = 0; x < s->xsize; x++)
		for(y = 0; y < s->ysize; y++)
			RG[getPoint(s,x,y,s->zsize -1)] = 0;
	
    return 1;
} /* make_reflections_grid */


/**
 * creates a syntetic seismic grid.
 */
int make_synthetic_grid(si *s, float *RG, float *SY) {
	printf_dbg2("make_syntethic_grid(): called\n");
	int x, y, z;
	int wavelet_spots;

	wavelet_spots = s->wavelet_used_values / 2;

	memset(SY,0,s->grid_size);

	int j;
	long it = 0;
	for(x = 0; x < s->xsize; x++) 
		for(y = 0; y < s->ysize; y++) 
			for(z = 0; z < s->zsize; z++) {
				it = z+1;
				for(j=-wavelet_spots+1; j<=wavelet_spots; j++) {
					if ((it+j >= 0) && (it+j < s->zsize-1)) {
						SY[getPoint(s,x,y,it+j)] += RG[getPoint(s,x,y,z)] * point_value(s, j);
					}
				}
			}

	return 1;
} /* make_synthetic_grid */



/**
 * 	
 * 	calculates correlations by layers between two cubes.
 */
int make_correlations_grid(si *s, float *seismic, float *syntetic, float *corr_grid) 
{
	printf_dbg2("make_correlations_grid(): called\n");

	int x, y, z;
	int i, z1, z2;
	int temp;
	
	struct layers_t *layers = s->layers;

	double sum_xy;
	double sum_x,sum_x_2,sum_x2;
	double sum_y,sum_y2,sum_y_2;
	double denom;

	float r;

	int n;

	for(x = 0; x < s->xsize; x++)
		for(y = 0; y < s->ysize; y++) {
			z1 = 0;
			for(i = 0; i < layers->number; i++) {
				z2 = z1+layers->layer[i];

				sum_x = 0;
				sum_x_2 = 0;
				sum_xy = 0;
				sum_y = 0;
				sum_y2 = 0;
				n = layers->layer[i];

				for(z = z1; z < z2; z++) {
					temp = getPoint(s,x,y,z);
					
					sum_x += seismic[temp];
					sum_x2 += seismic[temp] * seismic[temp];
					sum_xy += seismic[temp] * syntetic[temp];
					sum_y += syntetic[temp];
					sum_y2 += syntetic[temp] * syntetic[temp];

				}

				sum_x_2 = sum_x * sum_x;
				sum_y_2 = sum_y * sum_y;

				denom = sqrt(
						((n * sum_x2 ) - sum_x_2 ) * 
						((n * sum_y2 ) - sum_y_2)
						);

				if(denom == 0)
					r = 0; 
				else 
					r = ((n * sum_xy) - ( sum_x *sum_y ))/denom;

				for(z = z1; z < z2; z++) {
					corr_grid[getPoint(s,x,y,z)] = r; 
				}
				
				// reposition Z pointer
				z1=z2;
			}
		}
	/* divide corrAvg for number of points to give the average. */

	
    return 1;
} /* make_correlations_grid */



int expand_correlations_grid(cm_grid *cmg, float *CM) {
	printf_dbg("expand_correlations_grid(): called but not implemented\n");
    return 1;
} /* expand_correlations_grid */



float point_value(si* s, int point) 
{
    printf_dbg2("Point %d is at array position %d.\n", point, (s->wavelet_used_values / 2) + point);
    return s->values[(s->wavelet_used_values / 2) + point];
} /* point_value */



float index_value(si *s, int index) 
{
    return s->values[index];
} /* index_value */


/**
 * Generator of randomLayers, randomizes:
 * - number of layers
 * - size of each layer
 * also assures minimum_size, minimum_number, and total size
 */
int generateRandomLayers(struct layers_t * l)
{

	unsigned short int min_size = l->minimum_size;
	unsigned short int number = l->minimum_number;
	unsigned short int max_size = l->total_size;

	/*
	printf("generateRandomLayers: min_size: %d, min_number: %d, total_size: %d\n",min_size, number, max_size);
	*/
	int delta = max_size - (min_size * number);
	unsigned short int *lays, *temp;
	unsigned int i, x, sum;
	
	if(delta < 0) { 
		/* too many layers, or layers too big */
		fprintf(stderr,"randomLayers: ERROR, too many layers, or layers too big\n");
		return -1;
	} else if(delta == 0) {
		/* all layers must be of min_size */
		lays = (unsigned short int *) tsi_malloc(sizeof(unsigned short int) * number);
		for(i = 0; i < number; i++)
			lays[i] = min_size;
		l->number = number;
		l->layer = lays;
		return 0;
	}
	/* proper random Layers above...
	 * number of layers: number < N < max_size / min_size
	 */
	
	unsigned int n = max_size / min_size;
	temp = (unsigned short int *) tsi_malloc(n*sizeof(unsigned short int));
	for(i = 0; i < n; i++)
		temp[i] = 0;

	/* first, we warrantee the minimum number of layers */
	do {
		x = random() % (min_size * 2);
		x += min_size;
		i = (unsigned int) random() % n;

		temp[i] = x;

		x = 0;
		for(i = 0; i < n; i++) {
			if(temp[i] != 0) {
				x++;
			}
		}

	} while(  x < number); 

	/* now we assure the max_size */
	do {
		x = random() % (min_size);
		i = random() % n;
		if(temp[i] == 0)
			temp[i] += min_size;
		temp[i] += x;

		sum = 0;
		for(i = 0; i < n; i++)
			sum += temp[i];
	} while( sum < max_size );

	while(sum > max_size) {
		/* lets find the biggest layer */
		sum = 0;
		for(i = 0; i < n; i++) {
			if(temp[i] > sum) {
				sum = temp[i];
				x = i;
			}
		}
		/* now we reduce its size */
		temp[x] -= min_size;
		temp[x] /= 2;
		temp[x] += min_size;
		/* lets see if its good now */
		sum = 0;
		for(i = 0; i < n; i++)
			sum += temp[i];
	}
	if( sum != max_size) {
		/* add the diference to max_size */
		i = max_size - sum;
		temp[x] += i;
	}
	/* now lets count the layers != 0 */
	sum = 0;
	for(i = 0; i < n; i++) {
		if(temp[i] != 0)
			sum++;
	}

	/* this is the number of layers we are going to return */
	lays = (unsigned short int *) tsi_malloc(sizeof(unsigned short int) * sum);

	x = 0;
	for(i = 0; i < n; i++) {
		if(temp[i] != 0)
			lays[x++] = temp[i];
	}
	tsi_free(temp);

	l->number = sum;
	l->layer = lays;

	return 0;
}


void printLayers(struct layers_t *lay)
{
	int i;
	for(i = 0; i < lay->number; i++)
		printf("layer[%d] = %d\n",i,lay->layer[i]);
}

/* end of file si_math.c */
