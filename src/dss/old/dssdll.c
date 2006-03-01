/*                Sequential Gaussian Simulation */

/* The program is executed with no command line arguments.  The user */
/* will be prompted for the name of a parameter file.  The parameter */
/* file is described in the documentation (see the example sgsim.par) */

/* The output file will be a GEOEAS file containing the simulated values */
/* The file is ordered by x,y,z, and then simulation (i.e., x cycles */
/* fastest, then y, then z, then simulation number).  The values will be */
/* backtransformed to the original data values if a normal scores */
/* transform was performed. */
/* ----------------------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>

#include "dss.h"
#include "profile.h"


extern int readdata(float *, double *, int *,
		general_vars_t *, search_vars_t *, simulation_vars_t *);	
extern int readparam(float *, float *, 
		general_vars_t *, search_vars_t *, simulation_vars_t *,
		covariance_vars_t *, covtable_lookup_vars_t *);
extern int sdsim(float *, float *, float *, int *, int *,
		general_vars_t *,
		search_vars_t *,
		simulation_vars_t *,
		covariance_vars_t *,
		covtable_lookup_vars_t *,
		krige_vars_t *);

int dssdll(float *params, float *models, double * hard_data, int *hard_data_size,
		float *bcm_data, float *bai_data, int * mask_data, float *output_data)
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
//	static float * mask;
	int * order;
	
	general_vars_t general;
	search_vars_t search;
	simulation_vars_t simulation;
	covariance_vars_t covariance;
	covtable_lookup_vars_t covtable_lookup;
	krige_vars_t	krige_vars;

	/* Parameter adjustments */
	--output_data;
	--mask_data;
	--bai_data;
	--bcm_data;
	--hard_data;
	--models;
	--params;


	/* Function Body */
	general.lout = 2;


	newProfile();

#ifdef PROFILE
	profile.dssdll++;
	profBegin("dssdll");
#endif

	//copy pointers of aiCube, BCM and BAI
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
	//mask = (float *) malloc(general.nxyz * sizeof(float));
	order = (int *) malloc(general.nxyz * sizeof(int));

	/* Call sdsim for the simulation */
	sdsim(sim, bestAICube, bestCorrCube, order, &mask_data[1],
			&general, &search, &simulation, &covariance, 
			&covtable_lookup, &krige_vars);

	free(order);
	//free(mask);
	free(covtable_lookup.covtab);
	free(covtable_lookup.ixnode);
	free(covtable_lookup.iynode);
	free(covtable_lookup.iznode);

#ifdef PROFILE
	showResults();
	profEnd("dssdll");
#endif

	/* !Finished: */
	return 0;
} /* dssdll_ */

