/* si_math.c */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "log.h"
#include "registry.h"
#include "grid_heap.h"
#include "timer.h"
#include "si.h"
#include "si_math.h"

/* prototypes */
inline float point_value(si *s, int point);
inline float index_value(si *s, int idx);


inline unsigned int getPoint(si *s, int x, int y, int z)
{
//	return (x-1) + ( (y-1) * s->xsize) + ( (z-1) * s->xsize * s->ysize);
//	return ((z)* s->xsize * s->ysize + (y) * s->xsize + (x+1)) - 1;
	return (z * s->xsize * s->ysize + y * s->xsize + x);
}


/**
 * make Reflection Coeficients grid
 */
int make_reflections_grid(si *s, float *AI, float *RG) 
{
    int x, y, z, z_;
    double value;
    struct my_time t1, t2;
	unsigned int p1, p2;
        
    getCurrTime(&t1);

    for (x = 0; x < s->xsize; x++)
        for (y = 0; y < s->ysize; y++)
            for (z = 0; z < s->zsize -1; z++) {
				p1 = getPoint(s,x,y,z+1);
				p2 = getPoint(s,x,y,z);
                 value =  (AI[p1] - AI[p2]) / (AI[p1] + AI[p2]);
                 RG[p2] = value;
            }

    z_ = s->zsize-1;
    for (x = 0; x < s->xsize; x++)
        for (y = 0; y < s->ysize; y++) {
             RG[getPoint(s,x,y,z_)] = 0;
        }

    getCurrTime(&t2);

	log_action_time(s->l, 2, "make_reflections_grid() cpuTime",getElapsedTime(&t1,&t2));
	
    return 1;
} /* make_reflections_grid */



/**
 * creates a syntetic seismic grid.
 */
int make_synthetic_grid(si *s, float *RG, float *SY) {
    int x, y, z, nxy, nxyz;
    int wavelet_spots;
    int j, it, tmp;
    struct my_time t1, t2;
    unsigned int aux;

    aux = 0;
    printf_dbg2("make_synthetic_grid(): called");
    getCurrTime(&t1);
    wavelet_spots = s->wavelet_used_values / 2;
    nxy = s->xsize * s->ysize;
    nxyz = s->zsize * nxy;
    memset(SY, 0, s->grid_size * sizeof(float));

    /* original code */
    it = 0;
    for (x = 0; x < s->xsize; x++) 
        for (y = 0; y < s->ysize; y++) 
            for (z = 0; z < s->zsize; z++) {
                it = z + 1;
                for (j = -wavelet_spots+1; j <= wavelet_spots; j++) {
					tmp = it+j;
                    if ((tmp >= 0) && (tmp < s->zsize)) {
                        SY[getPoint(s, x, y, tmp)] += RG[getPoint( s, x, y, z)] * point_value(s, j);
                        aux++;
                    }
                }
            }

    getCurrTime(&t2);
	log_action_time(s->l, 2, "make_syntetic_grid() cpuTime",getElapsedTime(&t1,&t2));
    printf_dbg2("make_syntetic_grid (%f secs, %u calcs)\n", getElapsedTime(&t1, &t2), aux);

    return 1;
} /* make_synthetic_grid */



/**
 * 	
 * 	calculates correlations by layers between two cubes.
 */
int make_correlations_grid(si *s, float *seismic, float *synthetic, float *UNUSED)
{
    int x, y, z;
    int z1, z2;
    int temp;
    unsigned int i;

    double sum_xy;
    double sum_x,sum_x_2,sum_x2;
    double sum_y,sum_y2,sum_y_2;
    double denom;

    float r, *corr_grid;

	unsigned int nlayers, *layer_size;
    int n;
    struct my_time t1, t2;

    printf_dbg2("make_correlations_grid(): called\n");
    load_cmgrid(s->cmg);
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
                sum_y = 0;
                sum_x2 = 0;
                sum_y2 = 0;
                sum_xy = 0;
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

                denom = sqrt(
						((n * sum_x2 ) - sum_x_2 ) * 
						((n * sum_y2 ) - sum_y_2)
						);
                if (denom > 0)
                    r = ((n * sum_xy) - ( sum_x *sum_y ))/denom;
                else 
                    r = 0; 

				temp = getPoint(s,x,y,i);
				if (r < 0)
					corr_grid[temp] = 0; 
				else
					corr_grid[temp] = r; 
				
                // reposition Z pointer
                z1 = z2;
            }
        }
	
    //expand_correlations_grid(s->cmg, CM);
    getCurrTime(&t2);
	log_action_time(s->l, 2, "make_correlations_grid() cpuTime",getElapsedTime(&t1,&t2));
    dirty_cmgrid(s->cmg);
    return 1;
} /* make_correlations_grid */



int expand_correlations_grid(cm_grid *cmg, float *CM) {
    int z, j, n, nxy;
    float *cg;
    float *t;
    unsigned int i;

    printf_dbg2("expand_correlations_grid(): called\n");
    load_cmgrid(cmg);
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
    clear_cmgrid(cmg);
    printf_dbg2("expand_correlations_grid(): finished\n");
    return 0;
} /* expand_correlations_grid */


inline float point_value(si* s, int point) 
{
   // printf_dbg2("Point %d is at array position %d.\n", point, (s->wavelet_used_values / 2) + point);
    return (s->values[(s->wavelet_used_values / 2) + point]);
} /* point_value */


inline float index_value(si *s, int idx) 
{
    return s->values[idx];
} /* index_value */



/* end of file si_math.c */
