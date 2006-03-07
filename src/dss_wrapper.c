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
	int *order;
	
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

	/* 4 of the 6 grids are int * */
	order = (int *) d->g3;	
	d->covtable_lookup.covtab = d->g3;
	d->covtable_lookup.ixnode = (int *) d->g4;
	d->covtable_lookup.iynode = (int *) d->g5;
	d->covtable_lookup.iznode = (int *) d->g6;
   
    /* SIMULATION */
	readparam(params, variograms, 
			d->general, d->search, d->simulation, 
			d->covariance, d->covtable_lookup);
	readdata( NULL, wellsData, wellsDataSize, 
			d->general, d->search, d->simulation);
    
	/* read wells data */
	d->general.wellsNPoints = hard_data_size / general.nvari;
	d->general.wellsDataPos = (int *) malloc(general.wellsNPoints * sizeof(int));
	d->general.wellsDataVal = (double *) malloc(general.wellsNPoints * sizeof(double));
	readWellsData(d->general, wellsData, wellsDataSize);

	/* Call sdsim for the simulation */
	sdsim(sim, bestAICube, bestCorrCube, order, NULL,
			&general, &search, &simulation, &covariance, 
			&covtable_lookup, &krige_vars);

	free(general.wellsDataPos);
	free(general.wellsDataVal);
	
    delete_grid(d->heap, d->g1_idx);
    delete_grid(d->heap, d->g2_idx);
    delete_grid(d->heap, d->g3_idx);
    delete_grid(d->heap, d->g4_idx);
    delete_grid(d->heap, d->g5_idx);

} /* dss_simulation */



int codss_simulation(dss *d, float *currBAI, float *currBCM, float *AI) {
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
	int *order;
	float *sim;
	
	params = d->r->getParams();
	variograms = d->r->getVariograms();
	wellsData = d->r->getWellsData();
	wellsDataSize = d->r->getWellsDataSize();
	
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
    
	/* 4 of the 6 grids are int * */
	sim = d->g1;
	order = (int *) d->g3;	
	d->covtable_lookup.covtab = d->g3;
	d->covtable_lookup.ixnode = (int *) d->g4;
	d->covtable_lookup.iynode = (int *) d->g5;
	d->covtable_lookup.iznode = (int *) d->g6;
   
    /* SIMULATION */
	readparam(params, variograms, 
			d->general, d->search, d->simulation, 
			d->covariance, d->covtable_lookup);
	readdata( currBAI, wellsData, wellsDataSize, 
			d->general, d->search, d->simulation);
    
	/* read wells data */
	d->general.wellsNPoints = WellsDataSize / d->general->nvari;
	d->general.wellsDataPos = (int *) malloc(general.wellsNPoints * sizeof(int));
	d->general.wellsDataVal = (double *) malloc(general.wellsNPoints * sizeof(double));
	d->readWellsData(&general, wellsData, wellsDataSize);

	/* Call sdsim for the simulation */
	sdsim(sim, currBAI, currBCM, order, NULL,
			d->general, d->search, d->simulation, d->covariance, 
			d->covtable_lookup, d->krige_vars);

	free(d->general.wellsDataPos);
	free(d->general.wellsDataVal);
    
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

