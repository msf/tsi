#include <stdlib.h>
#include "debug.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss_wrapper.h"


dss *new_dss(registry *r, grid_heap *h) {
    dss *d;
    
    d = (dss *) malloc(sizeof(dss));
    if (d) {
        d->reg = r;
        d->heap = h;
        return d;
    }
    return NULL;
} /* new_dss */



int dss_simulation(dss *d, float *AI) {
/*    d->g1_idx = new_grid(d->heap);
    d->g2_idx = new_grid(d->heap);
    d->g3_idx = new_grid(d->heap);
    d->g4_idx = new_grid(d->heap);
    d->g5_idx = new_grid(d->heap);
    d->g1 = load_grid(d->heap, d->g1_idx);
    d->g2 = load_grid(d->heap, d->g1_idx);
    d->g3 = load_grid(d->heap, d->g1_idx);
    d->g4 = load_grid(d->heap, d->g1_idx);
    d->g5 = load_grid(d->heap, d->g1_idx);
*/   
    /* SIMULATION */
    
/*    delete_grid(d->heap, d->g1_idx);
    delete_grid(d->heap, d->g2_idx);
    delete_grid(d->heap, d->g3_idx);
    delete_grid(d->heap, d->g4_idx);
    delete_grid(d->heap, d->g5_idx);
*/
} /* dss_simulation */



int codss_simulation(dss *d, float *currBAI, float *currBCM, float *AI) {
/*    d->g1_idx = new_grid(d->heap);
    d->g2_idx = new_grid(d->heap);
    d->g3_idx = new_grid(d->heap);
    d->g4_idx = new_grid(d->heap);
    d->g5_idx = new_grid(d->heap);
    d->g1 = load_grid(d->heap, d->g1_idx);
    d->g2 = load_grid(d->heap, d->g1_idx);
    d->g3 = load_grid(d->heap, d->g1_idx);
    d->g4 = load_grid(d->heap, d->g1_idx);
    d->g5 = load_grid(d->heap, d->g1_idx);
*/    
    /* SIMULATION */
    
/*    delete_grid(d->heap, d->g1_idx);
    delete_grid(d->heap, d->g2_idx);
    delete_grid(d->heap, d->g3_idx);
    delete_grid(d->heap, d->g4_idx);
    delete_grid(d->heap, d->g5_idx);
*/
} /* codss_simulation */



void delete_dss(dss *d) {
    /* TODO */
} /* delete_dss */

