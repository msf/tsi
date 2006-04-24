#ifndef _SI_H
#define _SI_H

#include "debug.h"
#include "registry.h"
#include "grid_heap.h"

typedef struct cm_grid_type {
    int nlayers;      /* number of layers */
    int nxy;          /* size of each xy grid */
    int *layer_size;  /* size of each layer */
    float *cg;        /* compressed correlations grid */
} cm_grid;


typedef struct si_type {
    registry  *reg;
    grid_heap *heap;
    log_t *l;                   /* reference to the log */

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
    int random;
    int min_size;
    int min_number;

    /* correlations grid */
    cm_grid *cmg;
    
    /* execution log */
    int dump_sy,
        dump_rg;
    char *dump_path;
} si;


si *new_si(registry *r, grid_heap *h, log_t *l);

int setup_si(si *s);

int run_si(si *s, float *AI, float *seismic, float *CM, float *SY);

void delete_si(si *s);


cm_grid *new_cmgrid(si *s, int empty);

int build_cmgrid(cm_grid *g, int nlayers, int *layer_size);

cm_grid *load_cmgrid(si *s);

int get_nlayers(cm_grid *g);

int *get_layers(cm_grid *g);

int store_cmgrid(si *s, cm_grid *g);

void delete_cmgrid(cm_grid *g);

void print_layers(cm_grid *g);

#endif /* _SI_H */

