#ifndef TSI_H
#define TSI_H

#include "registry.h"
#include "grid_heap.h"
#include "dss_wrapper.h"
#include "si_wrapper.h"

#ifdef TSI_MPI
#include "tsi_parallel.h"
#else

typedef struct tsi_type {
    /* auxiliar objects */
    registry *reg;
    grid_heap *heap;

    /* geo-stats engines */
    dss *dss_eng;
    si  *si_eng;

    /* execution parameters*/
    int iterations,
        simulations,
        usefs;

    /* grid size */
    int xsize,
        ysize,
        zsize;
    unsigned int grid_size;

    /* ramdom layers */
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
} tsi;

tsi *new_tsi(registry *reg);

void delete_tsi(tsi *t);

int run_tsi(tsi *t);

#endif /* TSI_MPI */
#endif /* TSI_H */
