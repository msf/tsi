#ifndef TSI_H
#define TSI_H

#include "registry.h"
#include "grid_heap.h"
#include "dss_wrapper.h"
#include "si_wrapper.h"

typedef struct tsi_type {
    /* auxiliar objects */
    registry *reg;            /* reference to the registry */
    grid_heap *heap;          /* reference to the grid heap */

    /* geo-stats engines */
    dss *dss_eng;             /* reference to the DSS and CoDSS engines */
    si  *si_eng;              /* reference to the SI engine */

    /* execution parameters*/
    int iterations,           /* number of iterations */
        simulations;          /* number of simulations for each iteration */

    /* grid size */
    int xsize,
        ysize,
        zsize;
    unsigned int grid_size;


    /* random layers */
    int layers_min,
        size_layers_min;

    /* grid heap references */
    int seismic_idx,
        bestAI_idx,
        currBAI_idx,
        currBCM_idx,
        nextBAI_idx,
        nextBCM_idx,
        ai_idx,
        cm_idx,
        sy_idx;
        
    /* data pointers */
    float *seismic,
          *bestAI,
          *currBAI,
          *currBCM,
          *nextBAI,
          *nextBCM,
          *ai,
          *cm,
          *sy;

    /* parallel execution related data */
    unsigned int grid_segsize;   /* grid segment size = (int)grid_size/n_procs */
    int haveBestAI,
        n_procs,              /* number of processes running */
        proc_id;              /* process ID */
} tsi;

tsi *new_tsi(registry *reg);

void delete_tsi(tsi *t);

int run_tsi(tsi *t);

#endif /* TSI_H */
