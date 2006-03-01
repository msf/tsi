#ifndef _TSI_H
#define _TSI_H

#include "registry.h"
#include "grid_heap.h"

typedef struct tsi_type {
    /* auxiliar objects */
    registry *reg;
    grid_heap *heap;

    /* execution parameters*/
    int iterations,
        simulations,
        maxgrids;

    /* grid size */
    int xsize,
        ysize,
        zsize;

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
        cm_idx;
        
    /* data pointers */
    float *seismic,
          *bestAI,
          *currBAI,
          *currBCM,
          *nextBAI,
          *nextBCM,
          *cm;
} tsi;

tsi *new_tsi(registry *reg);

void delete_tsi(tsi *t);

int start_tsi(tsi *t);

#endif /* _TSI_H */
