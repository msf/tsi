#include <stdio.h>
#include "dss.h"
#include "memdebug.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)
	
/*
int dss(float *params, float *models, double *hardData, int *hDataSize, float *output_data)
{
	return coDss(params, models, hardData, hDataSize, NULL, NULL, output_data);
}
*/

int dsslib(float *params,
           float *models,
           double * hard_data, unsigned int hard_data_size,
		   float *bcm_data, float *bai_data,
           float *output_data)
{
	/* Cubos:
	 * bestCorrCube - BCM - copia do ponteiro de: bcm_data
	 * sim - aiCube - copia do ponteiro de output_data
	 * bestAICube - BAI
 	 * tmp - temporario, allocado e deallocado no covtable.
	 * order - ordem da simulacao
	 */
	float * bestCorrCube, * sim, * bestAICube;
	int * order;
	
	general_vars_t general;
	search_vars_t search;
	simulation_vars_t simulation;
	covariance_vars_t covariance;
	covtable_lookup_vars_t covtable_lookup;
	krige_vars_t	krige_vars;

	/* Function Body */
	general.lout = 2;


	/* copy pointers of aiCube, BCM and BAI */
	sim = output_data;
	bestCorrCube = bcm_data;
	bestAICube = bai_data;

	
	
	/* Read the parameters */
	readparam(params, models, &general, &search, &simulation, &covariance, &covtable_lookup);

	printf("dsslib(); nxyz: %d\n",general.nxyz); 

	/* read wells data */
	readdata(bestAICube, hard_data, hard_data_size, &general, &search, &simulation);
	readWellsData(&general, hard_data, hard_data_size);

	covtable_lookup.covtab = NULL;
	covtable_lookup.ixnode = covtable_lookup.iynode = covtable_lookup.iznode = order = NULL;
	
	covtable_lookup.covtab = (float *) my_malloc(general.nxyz * sizeof(float));
	covtable_lookup.ixnode = (int *) my_malloc(general.nxyz * sizeof(int));
	covtable_lookup.iynode = (int *) my_malloc(general.nxyz * sizeof(int));
	covtable_lookup.iznode = (int *) my_malloc(general.nxyz * sizeof(int));
	order = (int *) my_malloc(general.nxyz * sizeof(int));
	if(!covtable_lookup.covtab ||
			!covtable_lookup.ixnode ||
			!covtable_lookup.iynode ||
			!covtable_lookup.iznode ||
			!order) {
		printf("dsslib(): ERROR, unable to alocate grids\n");
		return 1;
	}

	/* Call sdsim for the simulation */
	sdsim(sim, bestAICube, bestCorrCube, order, NULL,
			&general, &search, &simulation, &covariance, 
			&covtable_lookup, &krige_vars);



	my_free(order);
	my_free(covtable_lookup.covtab);
	my_free(covtable_lookup.ixnode);
	my_free(covtable_lookup.iynode);
	my_free(covtable_lookup.iznode);

	/* these are alocated on readWellsData */
	my_free(general.wellsDataPos);
	my_free(general.wellsDataVal);

	/* !Finished: */
	return 0;
} /* dssdll_ */



