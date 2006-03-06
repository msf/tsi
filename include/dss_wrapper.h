#ifndef DSS_WRAPPER_H
#define DSS_WRAPPER_H

#include "debug.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss.h"

typedef struct dss_type {
    registry  *reg;
    grid_heap *heap;
	general_vars_t *general;
	search_vars_t *search;
	simulation_vars_t * simulation;
	covariance_vars_t *covariance;
	covtable_lookup_vars_t *covtable_lookup;
	krige_vars_t *krige;


    int g1_idx,
        g2_idx,
        g3_idx,
        g4_idx,
        g5_idx,
		g6_idx;

    float *g1,
          *g2,
          *g3,
          *g4,
          *g5,
		  *g6;
} dss;


dss *new_dss(registry *r, grid_heap *h);

int dss_simulation(dss *d, float *AI);

int codss_simulation(dss *d, float *currBAI, float *currBCM, float *AI);

void delete_dss(dss *d);

#endif /* DSS_WRAPPER_H */
