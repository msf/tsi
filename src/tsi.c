#include <stdlib.h>
#include <stdio.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
//#include "dss_wrapper.h"
//#include "si.h"
#include "tsi.h"
#ifdef TSI_MPI
#include "tsi_parallel.h"
#endif /* TSI_MPI */


int load_ascii_grid_from_file(float *grid, char *filename);

int load_binary_grid_from_file(float *grid, char *filename);

int load_segy_grid_from_file(float *grid, char *filename);


tsi *new_tsi(registry *reg) {
    tsi *t;
    reg_key *k;
    int n_grids;
    char *filename;

    t = (tsi *) malloc(sizeof(tsi));
    t->reg = reg;
    t->heap = NULL;

    /* get initial data from registry */
    k = get_key(reg, "GLOBAL", "ITERATIONS");
    if (k)
       t->iterations = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get number of iterations from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GLOBAL", "SIMULATIONS");
    if (k)
       t->simulations = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get number of simulations from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "XNUMBER");
    if (k)
       t->xsize = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get size of X from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "YNUMBER");
    if (k)
       t->ysize = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get size of Y from the registry!\n");
       delete_tsi(t);
       return NULL;
    }
    k = get_key(reg, "GRID", "ZNUMBER");
    if (k)
       t->zsize = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get size of Z from the registry!\n");
       delete_tsi(t);
       return NULL;
    }

    /* get heap data */
    k = get_key(reg, "GLOBAL", "GRIDS");
    if (k)
       t->maxgrids = get_int(k);
    else {
       printf_dbg("new_tsi(): failed to get max number of grids from the registry!\n");
       delete_tsi(t);
       return NULL;
    }

    t->heap = new_heap(11, t->maxgrids, t->xsize, t->ysize, t->zsize);
    if (!t->heap) {
       printf_dbg("new_tsi(): failed to start heap\n");
       delete_tsi(t);
       return NULL;
    }

/*
    k = get_key(reg, "SEISMIC", "FILENAME");
    if (k)
       filename = get_string(k);
    else {
    }
    t->seismic_idx = new_grid(t->heap);
    if (t->seismic_idx) {
        t->seismic = load_grid(t->heap, t->seismic_idx);
        if (load_ascii_grid_from_file(t->seismic, filename)) {
            dirty_grid(t->heap, t->seismic_idx);
        } else {
            printf_dbg("Failed to load seismic file!\n");
            delete_tsi(t);
            return NULL;
        }
    } else {
       printf_dbg("Failed to start the seismic grid!\n");
       delete_tsi(t);
       return NULL;
    }
*/

    return t;
} /* new_tsi */



void delete_tsi(tsi *t) {
} /* delete_tsi */



int start_tsi(tsi *t) {
    return 0;
} /* start_tsi */



int load_ascii_grid_from_file(float *grid, char *filename) {
    return 1;
} /* load_ascii_grid_from_file */



int load_binary_grid_from_file(float *grid, char *filename) {
    return 1;
} /* load_binary_grid_from_file */

/* end of file */
