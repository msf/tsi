#include <stdlib.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "si_wrapper.h"


si *new_si(registry *r, grid_heap *h) {
    si *s;
    
    s = (si *) malloc(sizeof(si));
    if (s) {
        s->reg = r;
        s->heap = h;
        return s;
    }

    return NULL;
} /* new_si */



int si_simulation(si *s, float *AI, float *seismic, float *CM, float *SY) {
    /* do calcs */
    return 1;
} /* si */



void delete_si(si *s) {
    /* TODO */
} /* delete_si */


