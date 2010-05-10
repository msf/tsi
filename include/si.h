#ifndef _SI_H
#define _SI_H

#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "log.h"

typedef struct cm_grid_type {
    unsigned int nxy;          /* size of each xy grid */
    unsigned int nlayers;      /* number of layers */
    unsigned int *layer_size;  /* size of each layer */
    float *cg;        /* compressed correlations grid */
    int cg_idx;
    grid_heap *heap;
} cm_grid;


typedef struct si_type {
    registry  *reg;
    grid_heap *heap;
    log_t *l;                   /* reference to the log */

    int n_procs,
        proc_id;

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
    char *dump_path,
         empty_path;
    int dump_file;        /* dump file type */
} si;


/* SI methods */
si *new_si(registry *r, grid_heap *h, log_t *l, int n_procs, int proc_id);

int run_si(si *s, float *AI, float *seismic, float *CM, float *SY, int it, int sim);

void delete_si(si *s);

/* compressed correlations grid methods */
cm_grid *new_cmgrid(si *s, int empty_flag);

void delete_cmgrid(cm_grid *g);

int load_cmgrid(cm_grid *g);

int dirty_cmgrid(cm_grid *g);

int clear_cmgrid(cm_grid *g);

int store_cmgrid(si *s, cm_grid *g);

int build_cmgrid(cm_grid *g, unsigned int nlayers, unsigned int *layer_size);

cm_grid *clone_cmgrid(cm_grid *g);

cm_grid *get_cmgrid(si *s);

unsigned int get_nlayers(cm_grid *g);

unsigned int *get_layers(cm_grid *g);

void print_layers(log_t *l, cm_grid *g);



#endif /* _SI_H */
