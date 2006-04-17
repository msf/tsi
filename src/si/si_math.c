/* si_math.c */

#include <math.h>
#include <stdlib.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "timer.h"
#include "si.h"


/* wavelet functions */
float point_value(si *s, int point);
float index_value(si *s, int index);


unsigned int getPoint(si *s, int x, int y, int z)
{
//	return (x-1) + ( (y-1) * s->xsize) + ( (z-1) * s->xsize * s->ysize);
//	return ((z)* s->xsize * s->ysize + (y) * s->xsize + (x+1)) - 1;
	return (z * s->xsize * s->ysize + y * s->xsize + x);
}

int make_reflections_grid(si *s, float *AI, float *RG) 
{
    unsigned int x, y, z;
    double value;
    struct timeval t1, t2;
        
    printf_dbg("make_reflections_grid(): called");
    getCurrTime(&t1);
    for (x = 0; x < s->xsize; x++)
        for (y = 0; y < s->ysize; y++)
            for (z = 0; z < s->zsize -1; z++) {
                 value =  AI[getPoint(s,x,y,z+1)] - AI[getPoint(s,x,y,z)];
                 value /= AI[getPoint(s,x,y,z+1)] + AI[getPoint(s,x,y,z)];
                 RG[getPoint(s,x,y,z)] = value;
            }
    for (x = 0; x < s->xsize; x++)
        for (y = 0; y < s->ysize; y++)
             RG[getPoint(s,x,y,s->zsize -1)] = 0;
    getCurrTime(&t2);
    printf_dbg(" (%f secs)\n", getElapsedTime(&t1, &t2));
	
    return 1;
} /* make_reflections_grid */


/**
 * creates a syntetic seismic grid.
 */
int make_synthetic_grid(si *s, float *RG, float *SY) {
    int x, y, z, nxy, nxyz, px;
    int wavelet_spots;
    float *wvals, wval;
    int j, it;
    struct timeval t1, t2;
    unsigned int aux;

    aux = 0;
    printf_dbg("make_synthetic_grid(): called");
    getCurrTime(&t1);
    wavelet_spots = s->wavelet_used_values / 2;
    nxy = s->xsize * s->ysize;
    nxyz = s->zsize * nxy;
    memset(SY,0,s->grid_size);


    it = 0;
    for (x = 0; x < s->xsize; x++) 
        for (y = 0; y < s->ysize; y++) 
            for (z = 0; z < s->zsize; z++) {
                it = z + 1;
                for (j = -wavelet_spots+1; j <= wavelet_spots; j++) {
                    if ((it+j >= 0) && (it+j < s->zsize-1)) {
                        SY[getPoint(s,x,y,it+j)] += RG[getPoint(s,x,y,z)] * point_value(s, j);
                        aux++;
                    }
                }
            }


/*
    it = 0;
    for (x = 0; x < s->xsize; x++) 
        for (y = 0; y < s->ysize; y++) 
            for (z = 0; z < s->zsize; z++) {
                it = z + 1;
                for (j = -wavelet_spots+1; j <= wavelet_spots; j++) {
                    if ((it+j >= 0) && (it+j < s->zsize-1)) {
                        SY[(it+j) * nxy + y * s->xsize + x] += RG[z * nxy + y * s->xsize + x] * s->values[wavelet_spots+j];
                        aux++;
                    }
                }
            }
*/

/*
    it = 0;
    wvals = s->values + wavelet_spots;
    for (j = -wavelet_spots + 1; j <= wavelet_spots; j++) {
        for (z = 0; z < s->zsize; z++) {
            it = z + 1 + j;
            if ((it > 0) && (it < s->zsize)) {
                wval = *(wvals + j);
                for (x = 0; x < s->xsize; x++) 
                    for (y = 0; y < s->ysize; y++) {
                        SY[getPoint(s,x,y,it)] += RG[getPoint(s,x,y,z)] * wval;
                        aux++;
                    }
            }
        }
    }
*/
    getCurrTime(&t2);
    printf_dbg(" (%f secs, %u calcs)\n", getElapsedTime(&t1, &t2), aux);

    return 1;
} /* make_synthetic_grid */



/**
 * 	
 * 	calculates correlations by layers between two cubes.
 */
int make_correlations_grid(si *s, float *seismic, float *synthetic, float *CM) 
{
    int x, y, z;
    int i, z1, z2;
    int temp;

    double sum_xy;
    double sum_x,sum_x_2,sum_x2;
    double sum_y,sum_y2,sum_y_2;
    double denom;

    float r, *corr_grid;

    int n, nlayers, *layer_size;
    struct timeval t1, t2;

    printf_dbg("make_correlations_grid(): called\n");
    getCurrTime(&t1);
    corr_grid = s->cmg->cg;   /* compressed grid */
    nlayers = s->cmg->nlayers;
    layer_size = s->cmg->layer_size;
        
    for (x = 0; x < s->xsize; x++)
        for (y = 0; y < s->ysize; y++) {
            z1 = 0;
            for (i = 0; i < nlayers; i++) {
                z2 = z1+layer_size[i];

                sum_x = 0;
                sum_x_2 = 0;
                sum_xy = 0;
                sum_y = 0;
                sum_y2 = 0;
                n = layer_size[i];

                for (z = z1; z < z2; z++) {
                    temp = getPoint(s,x,y,z);

                    sum_x += seismic[temp];
                    sum_x2 += seismic[temp] * seismic[temp];
                    sum_xy += seismic[temp] * synthetic[temp];
                    sum_y += synthetic[temp];
                    sum_y2 += synthetic[temp] * synthetic[temp];

                }

                sum_x_2 = sum_x * sum_x;
                sum_y_2 = sum_y * sum_y;

                denom = sqrt(((n * sum_x2 ) - sum_x_2 ) * ((n * sum_y2 ) - sum_y_2));
                if (denom > 0)
                    r = ((n * sum_xy) - ( sum_x *sum_y ))/denom;
                else 
                    r = 0; 

                corr_grid[getPoint(s,x,y,i)] = r; 
				
                // reposition Z pointer
                z1 = z2;
            }
        }
       
    expand_correlations_grid(s->cmg, CM);
    getCurrTime(&t2);
    printf_dbg("make_synthetic_grid(): finished (%f secs)\n", getElapsedTime(&t1, &t2));
	
    return 1;
} /* make_correlations_grid */



int expand_correlations_grid(cm_grid *cmg, float *CM) {
    int z, i, j, n, nxy;
    float *cg;
    float *t;

    printf_dbg("expand_correlations_grid(): called\n");
    nxy = cmg->nxy;
    z = 0;
    for (i = 0; i < cmg->nlayers; i++) {
        cg =  &cmg->cg[nxy*i]; // next layer
        n = cmg->layer_size[i];
        for (j = 0; j < n; j++) {
	    t = &CM[(nxy*z) + (nxy*j)];
	    printf_dbg2("expand...: z=%d, i=%d\n", z+j, i);
            memcpy(t, cg, (nxy * sizeof(float)));
        }
        z += n;
    }
    printf_dbg("expand_correlations_grid(): finished\n");
    return 0;
} /* expand_correlations_grid */



float point_value(si* s, int point) 
{
    printf_dbg2("Point %d is at array position %d.\n", point, (s->wavelet_used_values / 2) + point);
    return (s->values[(s->wavelet_used_values / 2) + point]);
} /* point_value */



float index_value(si *s, int index) 
{
    return s->values[index];
} /* index_value */



/* end of file si_math.c */
