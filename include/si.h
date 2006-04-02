#ifndef _SI_H
#define _SI_H

#include "debug.h"
#include "registry.h"
#include "grid_heap.h"

typedef struct cm_grid_type {
    int nlayers;      /* number of layers */
    int *layer_size;  /* size of each layer */
    float *cg;        /* compressed correlations grid */
    int nxy;          /* size of each xy grid */
} cm_grid;


typedef struct si_type {
    registry  *reg;
    grid_heap *heap;

	/* grid size */
    int xsize,
        ysize,
        zsize;
    unsigned int grid_size;   /* = xsize * ysize * zsize */

    /* wavelet data */
    int   *points;
    float *values;
    int   wavelet_used_values,
          max_values;  /* max number of wavelet values (for point > 0) */
          
	/* layers data */
	int min_size;
	int min_number;

    /* correlations grid */
    cm_grid *cmg;
} si;


si *new_si(registry *r, grid_heap *h);

int setup_si(si *s);

int run_si(si *s, float *AI, float *seismic, float *CM, float *SY);

void delete_si(si *s);


cm_grid *new_cmgrid(si *s);

cm_grid *load_cmgrid(si *s);

void store_cmgrid(si *s, cm_grid *g);

void delete_cmgrid(cm_grid *g);

#endif /* _SI_H */

