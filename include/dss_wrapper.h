#include "debug.h"
#include "grid_heap.h"
#include "registry.h"

typedef struct dss_type {
    registry  *reg;
    grid_heap *heap;

    int g1_idx,
        g2_idx,
        g3_idx,
        g4_idx,
        g5_idx;

    float *g1,
          *g2,
          *g3,
          *g4,
          *g5;
} dss;


dss *new_dss(registry *r, grid_heap *h);

int dss_simulation(dss *d, float *AI);

int codss_simulation(dss *d, float *currBAI, float *currBCM, float *AI);

void delete_dss(dss *d);

