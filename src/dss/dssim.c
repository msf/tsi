/* -----------------------------------------------------------------------
   Direct Sequential Simulation and co-Simulation of a 3-D Rectangular Grid
 * ----------------------------------------------------------------------- */


/* This subroutine generates 3-D realizations of a Gaussian process with */
/* a given autocovariance model, and conditional to input Gaussian data. */
/* The conditional simulation is achieved by sequential simulation of all */
/* the nodes visited by a random path. */

/* PROGRAM NOTES: */

/*  1. The three dimensional anisotropy parameters, i.e., of the search */
/*     ellipse and variogram ranges are described in section 2.3 of the */
/*     manual.   The variogram parameters are described in the same place */

/*  2. The original data and previously simulated grid nodes can be */
/*     searched separately.  There can be a different maximum number of */
/*     each and a minimum number of original data can be specified */
/*     to restrict simulation beyond the limits of the data.  The */
/*     closeness of previously simulated grid nodes is measured according */
/*     to the variogram structural distance. */

/* INPUT VARIABLES: */
/*   nd               Number of data (no missing values) */
/*   x,y,z(nd)        coordinates of the data */
/*   vr(nd)           gaussian data (normal scores) */
/*   nx,ny,nz         Number of blocks in X,Y, and Z */
/*   xmn,ymn,zmn      Coordinate at the center of the first Block */
/*   xsiz,ysiz,zsiz   spacing of the grid nodes (block size) */
/*   nsim             number of simulations */
/*   ktype            =1, ordinary kriging; =0, simple kriging */
/*   sim              the current realization */
/*   idbg             integer debugging level (0=none,2=normal,4=serious) */
/*   ldbg             unit number for the debugging output */
/*   lout             unit number for the output */
/*   radius           Maximum search radius */
/*   sang1            Azimuth angle of the principal search direction */
/*   sang2            Dip angle of the principal search direction */
/*   sang3            Third rotation angle of the search ellipse */
/*   sanis1           Anisotropy for the dip angle */
/*   sanis2           Anisotropy for the plunge angle */
/*   ndmin            Minimum number of data required before sim */
/*   ndmax            Maximum number of samples for simulation */
/*   noct             Maximum number per octant if an octant search is */
/*                      desired (if <= 0, then no octant search) */
/*   nodmax           Maximum number of previously simulated grid nodes */
/*                      to consider in the simulation.  The structural */
/*                      variogram distance is used to identify close ones */
/*   c0               Nugget constant (isotropic). */
/*   cc(nst)          Multiplicative factor of each nested structure. */
/*   aa(nst)          Parameter "a" of each nested structure. */
/*   it(nst)          Type of nested structures (1=sph,2=exp,3=gau,4=pow) */
/*   ang1(nst)        Azimuth angle for the principal direction */
/*   ang2(nst)        Dip angle for the principal direction */
/*   ang3(nst)        Third rotation angle to rotate the two minor */
/*                      directions around the principal direction */
/*   anis1(nst)       Anisotropy (radius in minor direction at 90 */
/*                      degrees from "ang1" divided by the principal */
/*                      radius in direction "ang1") */
/*   anis2(nst)       Anisotropy (radius in minor direction at 90 degrees */
/*                      vertical from "ang1" divided by the principal */
/*                      radius in direction "ang1") */

/* OUTPUT VARIABLES:  Simulated Values are written to "lout"
 *   sim 			The current simulation is the output of sdsim */

/* EXTERNAL REFERENCES: */
/*   super            Sets up the super block search of original data */
/*   search           Search for nearby data values */
/*   covtable           Builds a covariance table and "spiral" search */
/*   srchnd           Search for nearby simulated grid nodes */
/*   sqdist           computes anisotropic squared distance */
/*   sortem           sorts multiple arrays in ascending order (separate) */
/*   krige            Sets up and solves either the SK or OK system */
/*   ksol             Linear system solver using Gaussian elimination */
/* Concepts taken from F. Alabert and E. Isaaks */
/* ----------------------------------------------------------------------- */

/** Funcoes utilizadas
 * covtable (covtable)
 * gauinv() /
 * getindx
 * krige
 * setrot
 * sortem
 * srchsupr
 * srchnod
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "math_random.h"
#include "dss.h"
#include "dss_legacy.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


/* Table of constant values */


int dssim(float *sim, float *bestAICube, float *bestCorrCube, int *order, int ktype,
        mask_t                  *   mask,
		general_vars_t			*	general,
		harddata_t				*	harddata,
		search_vars_t			*	search,
		covariance_vars_t		*	covariance,
		covtable_lookup_vars_t	*	covtable_lookup,
		krige_vars_t			*	krige_vars)
{

	unsigned int toSim;
	unsigned int i;


	/* Local variables */
	int kinicial;
	double p, zvariance;
	int in;
	int ix, iz, iy;
	float xp, xx, yy, zz;
	float vmy = 0, vms;
	int ierr;
	float cmean, global_mean;

	double zmean;
	int index;
	float std_deviation;
	int lktype;
	float clcorr;
	float simval = 0;


	/* Parameter adjustments */
	--order;
	--bestCorrCube;
	--bestAICube;
	--sim;

	printf_dbg2("\tdssim() called\n");

	/* Function Body */
	for (i = 0; i < covariance->varnum; ++i) {
		setrot( covariance->variogram[i].ang1,
				covariance->variogram[i].ang2,
				covariance->variogram[i].ang3,
				covariance->variogram[i].anis1,
				covariance->variogram[i].anis2,
				i, krige_vars->rotmat);
	}
	covariance->isrot = 4;
	setrot(search->sang1, search->sang2, search->sang3,
			search->sanis1, search->sanis2, covariance->isrot,
			krige_vars->rotmat);

	printf_dbg2("dssim() calling covtable\n");

	/* !Set up the covariance table and the spiral search: */
    covtable_lookup->covtab = (float *) tsi_malloc( sizeof(float) * general->nxyz);
    covtable_lookup->ixnode = (short *) tsi_malloc( sizeof(short) * general->nxyz);
    covtable_lookup->iynode = (short *) tsi_malloc( sizeof(short) * general->nxyz);
    covtable_lookup->iznode = (short *) tsi_malloc( sizeof(short) * general->nxyz);
	covtable(general, search, covariance, covtable_lookup, krige_vars);

    order = (int *) tsi_malloc( sizeof(int) * general->nxyz);

	/* prepare the random path */
	/* init the sim table with NOSVALUE (no simulated value) */
	for(i = 1; i <= general->nxyz; ++i) {
		order[i] = i;
		sim[i] = general->nosim_value;
	}

	toSim = general->nxyz;
	index = in = 0;

	/* copy hard data to simulation grid */
	for(i = 0; i < harddata->in_grid_count; i++) {
		in = harddata->in_grid[i].index;
		sim[in] = harddata->in_grid[i].value;

		// mark point has simulated
		order[in] = order[toSim];
		order[toSim] = in;
		toSim--;
	}


	printf_dbg2("\tdssim(): Starting simulation now\n");
	/* !MAIN LOOP OVER ALL THE NODES: */
	ierr = 0;
	in = 0;
	zmean = 0.f;
	zvariance = 0.f;
	while( toSim > 0 ) {

		/* generate point to simulate */
		in = ((int) tsi_random() % toSim) +1; /* +1 because of fortran offsets */
		index = order[in];  /* point to be simulated */

		/* mark has simulated */
		order[in] = order[toSim];
		order[toSim] = index;
		toSim--;


#ifdef TSI_DEBUG
		if(toSim == (general->nxyz * 0.75)) {
			printf_dbg("\tdsssim(): 1/4 completed.\n");
		} else if(toSim == (general->nxyz / 2)) {
			printf_dbg("\tdsssim(): 1/2 completed.\n");
		} else if(toSim == (general->nxyz / 4)) {
			printf_dbg("\tdsssim(): 3/4 completed.\n");
 		}

		if(index < 1)
			printf("dssim(): ERROR, INDEX(%d) < 1\n", index);
		/* if value has a value allready (like in the case of a hard data guiven point), skip simulation */
		if( sim[index] >= harddata->min_value &&
		    sim[index] <= harddata->max_value) {
			++ierr;
			printf("dssim(): ERROR: index(%d) with valid data: %f\n", index, sim[index]);
		}

		if (sim[index] > general->nosim_value + 1e-20f ||
				 sim[index] < general->nosim_value * 2.f) {
			++ierr;
			printf("dssim(): ERROR: index(%d) allready has data: %f\n",index,sim[index]);
			continue;
		}
#endif


		/* skip points that belong to mask */
		if( mask && mask_isset(mask, index) )
			continue;


		/* get relative x,y,z from index */
		get3Dcoords(index, general->nx, general->nxy, &ix, &iy, &iz);

		/* getting absolute coords from relative coords */
		xx = getAbsolutePos(general->xmn, general->xsiz, ix);
		yy = getAbsolutePos(general->ymn, general->ysiz, iy);
		zz = getAbsolutePos(general->zmn, general->zsiz, iz);

		/* !Now, we'll simulate the point ix,iy,iz.  First, get the close data */
		/* !and make sure that there are enough to actually simulate a value, */
		/* !we'll only keep the closest "ndmax" data, and look for previously */
		/* !simulated grid nodes: */

		covtable_lookup->ncnode = srchnod(ix, iy, iz, &sim[1],
				general,
				covtable_lookup,
				covtable_lookup->node);
		/* !WARNING:Para ter em atencao; bai c/ NOSIMVALUE, do simple kriging*/
		kinicial = ktype;
		if (ktype == CO_KRIG && bestAICube[index] == general->nosim_value) {
			ktype = SIMPLE_KRIG;
		}
		/* !Calculate the conditional mean and standard deviation.  This will be */
		/* !done with kriging if there are data, otherwise, the global mean and */
		/* !standard deviation will be used: */
		if ( ktype == CO_KRIG) {
			global_mean = bestAICube[index];
		} else {
			global_mean = harddata->average;
		}
		if (covtable_lookup->ncnode < 1) {
			cmean = harddata->average;
			std_deviation = sqrt(harddata->variance);
		} else {
			/* !Perform the kriging.  Note that if there are fewer than four data */
			/* !then simple kriging is prefered so that the variance of the */
			/* !realization does not become artificially inflated: */
			lktype = ktype;
			if ( ktype == ORDINARY_KRIG && covtable_lookup->ncnode < 4) {
				lktype = SIMPLE_KRIG;
			}
			/* !Estimacao em xo (SDSIM) */
			/* aceder ao bestCorrCube, apenas se ktype >= 4 */
			if(ktype == CO_KRIG) {
				clcorr = bestCorrCube[index];
			}

			krige(ix, iy, iz, xx, yy, zz, lktype, global_mean,
					&cmean, &std_deviation, // these are the output vars of krige
					&bestAICube[1], clcorr,
					general, harddata,
					covariance, covtable_lookup, krige_vars, covtable_lookup->node);

		}
		general->ktype = kinicial;
		vmy = compute_gaussian_equiv( cmean, harddata->point_count, harddata->point);

		/* !Gera um valor aleatorio com distribuicao Gaussiana */
		vms = 0;
		for (i = 0; i < covtable_lookup->ntry; ++i) {
			p = tsi_random_real();
			gauinv(p, &xp);
			xp = xp * std_deviation + vmy;
			/* !Transformada inversa final (iv)               (SDSIM) */
			simval = backtr(xp, harddata->point_count, harddata->point,
					harddata->min_value, harddata->max_value);
			vms += simval;
		}
		vms /= covtable_lookup->ntry;
		/* !Reter o valor simulado  (SDSIM) */
		if (covtable_lookup->ntry > 1) {
			sim[index] = simval + (cmean - vms);
		} else {
			sim[index] = simval;
		}
		if (sim[index] < harddata->min_value) {
			sim[index] = harddata->min_value;
		} else if (sim[index] > harddata->max_value) {
			sim[index] = harddata->max_value;
		}
	}	/* !END MAIN LOOP OVER NODES: */

	/* free aux grids */
	tsi_free(order);
	tsi_free(covtable_lookup->ixnode);
	tsi_free(covtable_lookup->iynode);
	tsi_free(covtable_lookup->iznode);
	tsi_free(covtable_lookup->covtab);

	printf_dbg("dssim(): DEBUG: SKIP points: %d\n",ierr);

#ifdef TSI_DEBUG
	ierr = 0;
	for(i = 0; i < harddata->in_grid_count; i++) {
		in = harddata->in_grid[i].index;
		simval = sim[in] - harddata->in_grid[i].value;
		if(simval != 0 ) {
			printf_dbg("sim[%d] - harddata = %f\n",in,simval);
			ierr++;
		}
	}
	printf_dbg("dssim(): sim grid disrespects %d wells data points\n",ierr);
#endif



	/* !Return to the main program: */
	return 0;
} /* sdsim_ */

