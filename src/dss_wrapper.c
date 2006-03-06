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
		d->general = (general_vars_t *) malloc(sizeof(general_vars_t));
		d->search = (search_vars_t *) malloc(sizeof(search_vars_t));
		d->simulation = (simulation_vars_t *) malloc(sizeof(simulation_vars_t));
		d->covariance = (covariance_vars_t *) malloc(sizeof(covariance_vars_t));
		d->covtable_lookup = (covtable_lookup_vars_t *) malloc(sizeof(covtable_lookup_vars_t));
		d->krige = (krige_vars_t *) malloc(sizeof(krige_vars_t));

        return d;
    }
    return NULL;
} /* new_dss */



int run_dss(dss *d, float *AI) {
	/* Cubos:
	 * sim - aiCube - copia do ponteiro de output_data
	 * order - ordem da simulacao

	covtable_lookup.covtab = (float *) malloc(general.nxyz * sizeof(float));
	covtable_lookup.ixnode = (int *) malloc(general.nxyz * sizeof(int));
	covtable_lookup.iynode = (int *) malloc(general.nxyz * sizeof(int));
	covtable_lookup.iznode = (int *) malloc(general.nxyz * sizeof(int));
	 * 
	 */

	float *params,*variograms;
	double *wellsData;
	int wellsDataSize;

	params = d->r->getParams();
	variograms = d->r->getVariograms();
	wellsData = d->r->getWellsData();
	wellsDataSize = d->r->getWellsDataSize();
	
    d->g1_idx = new_grid(d->heap);
    d->g2_idx = new_grid(d->heap);
    d->g3_idx = new_grid(d->heap);
    d->g4_idx = new_grid(d->heap);
    d->g5_idx = new_grid(d->heap);
    d->g6_idx = new_grid(d->heap);


    d->g1 = load_grid(d->heap, d->g1_idx);
    d->g2 = load_grid(d->heap, d->g2_idx);
    d->g3 = load_grid(d->heap, d->g3_idx);
    d->g4 = load_grid(d->heap, d->g4_idx);
    d->g5 = load_grid(d->heap, d->g5_idx);
    d->g6 = load_grid(d->heap, d->g6_idx);
   
    /* SIMULATION */
    
    delete_grid(d->heap, d->g1_idx);
    delete_grid(d->heap, d->g2_idx);
    delete_grid(d->heap, d->g3_idx);
    delete_grid(d->heap, d->g4_idx);
    delete_grid(d->heap, d->g5_idx);

} /* dss_simulation */



int codss_simulation(dss *d, float *currBAI, float *currBCM, float *AI) {
    d->g1_idx = new_grid(d->heap);
    d->g2_idx = new_grid(d->heap);
    d->g3_idx = new_grid(d->heap);
    d->g4_idx = new_grid(d->heap);
    d->g5_idx = new_grid(d->heap);
    d->g1 = load_grid(d->heap, d->g1_idx);
    d->g2 = load_grid(d->heap, d->g2_idx);
    d->g3 = load_grid(d->heap, d->g3_idx);
    d->g4 = load_grid(d->heap, d->g4_idx);
    d->g5 = load_grid(d->heap, d->g5_idx);
    
    /* SIMULATION */
    
    delete_grid(d->heap, d->g1_idx);
    delete_grid(d->heap, d->g2_idx);
    delete_grid(d->heap, d->g3_idx);
    delete_grid(d->heap, d->g4_idx);
    delete_grid(d->heap, d->g5_idx);

} /* codss_simulation */



void delete_dss(dss *d) {
	
    /* TODO */
	free(d->general);
	free(d->search);
	free(d->simulation);
	free(d->covariance);
	free(d->covtable_lookup);
	free(d->krige);
} /* delete_dss */

