#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)

int dsslib(float *params,
           float *models,
           double * hard_data, int *hard_data_size,
	   float *bcm_data, float *bai_data,
           float *output_data)
{
   return dssdll(params, models, hard_data, hard_data_size, bcm_data, bai_data, NULL, output_data);
}

int dssdll(float *params,
           float *models,
           double * hard_data, int *hard_data_size,
	   float *bcm_data, float *bai_data,
           int * mask_data,
           float *output_data)
{
	/* Cubos:
	 * bestCorrCube - BCM - copia do ponteiro de: bcm_data
	 * sim - aiCube - copia do ponteiro de output_data
	 * bestAICube - BAI
 	 * tmp - temporario, allocado e deallocado no covtable.
	 * mask - mascara, nao esta a ser usado
	 * order - ordem da simulacao
	 */
	float * bestCorrCube, * sim, * bestAICube;
	/* static float * mask; */
	int * order;
	
	general_vars_t general;
	search_vars_t search;
	simulation_vars_t simulation;
	covariance_vars_t covariance;
	covtable_lookup_vars_t covtable_lookup;
	krige_vars_t	krige_vars;

	/* Parameter adjustments */
	--output_data;
	if(mask_data) --mask_data;
	--bai_data;
	--bcm_data;
	--hard_data;
	--models;
	--params;


	/* Function Body */
	general.lout = 2;



#ifdef PROFILE
	newProfile();
	profile.dssdll++;
	profBegin("dssdll");
#endif

	/* copy pointers of aiCube, BCM and BAI */
	sim = &output_data[1];
	bestCorrCube = &bcm_data[1];
	bestAICube = &bai_data[1];

	/*      !ldbg = 3 */
	/*      !lbestAICube = 4 */

	/* Read the parameters */
	readparam(&params[1], &models[1], &general, &search, &simulation, &covariance, &covtable_lookup);

	readdata(bestAICube, &hard_data[1], hard_data_size, &general, &search, &simulation);

	covtable_lookup.covtab = (float *) malloc(general.nxyz * sizeof(float));
	covtable_lookup.ixnode = (int *) malloc(general.nxyz * sizeof(int));
	covtable_lookup.iynode = (int *) malloc(general.nxyz * sizeof(int));
	covtable_lookup.iznode = (int *) malloc(general.nxyz * sizeof(int));
	/* mask = (float *) malloc(general.nxyz * sizeof(float)); */
	order = (int *) malloc(general.nxyz * sizeof(int));

	/* Call sdsim for the simulation */
	/* sdsim(sim, bestAICube, bestCorrCube, order, &mask_data[1], */
	sdsim(sim, bestAICube, bestCorrCube, order, NULL,
			&general, &search, &simulation, &covariance, 
			&covtable_lookup, &krige_vars);

	free(order);
	/* free(mask); */
	free(covtable_lookup.covtab);
	free(covtable_lookup.ixnode);
	free(covtable_lookup.iynode);
	free(covtable_lookup.iznode);

#ifdef PROFILE
	profEnd("dssdll");
	showResults();
#endif

	/* !Finished: */
	return 0;
} /* dssdll_ */


