#include "debug.h"
#include "registry.h"
#include "grid_heap.h"

typedef struct si_type {
    registry  *reg;
    grid_heap *heap;

    int g;
} si;


si *new_si(registry *r, grid_heap *h);

int si_simulation(si *s, float *AI, float *seismic, float *CM);

void delete_si(si *s);


