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
 * picksup (apenas se search->sstrat == 0)
 * setsupr (apenas se search->sstrat == 0)
 * setrot
 * sortem
 * srchsupr
 * srchnod
 */

/** CUBOS utilizados
 * sim
 * bestCorrCube
 * bestAICube
 * order
 * ??
 * 
 */ 

/** CUBOS _nao_ utilizados
 * tmp - e' apenas passado para covtable(). 
 * -
 */

/** structs globais utilizadas:
 * generl_1
 * search_1
 * simula_1
 * cova3d_1
 * clooku_1
 * krige_varsv_1
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
		general_vars_t			*	general,
		search_vars_t			*	search,
		simulation_vars_t		*	simulation,
		covariance_vars_t		*	covariance,
		covtable_lookup_vars_t	*	covtable_lookup,
		krige_vars_t			*	krige_vars)
{

	unsigned int toSim;
	int i;


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
	for (i = 1; i <= covariance->varnum; ++i) {
		setrot(	covariance->variogram[i-1].ang1,
				covariance->variogram[i-1].ang2,
				covariance->variogram[i-1].ang3,
				covariance->variogram[i-1].anis1,
				covariance->variogram[i-1].anis2,
				i, krige_vars->rotmat);
	}
	covariance->isrot = 5;
	setrot(search->sang1, search->sang2, search->sang3,
			search->sanis1, search->sanis2, covariance->isrot, 
			krige_vars->rotmat);
	/* !Set up the super block search: */
    /* TSI NOTE: SSTRAT IS ALLWAYS 1 */
    /*
	if (search->sstrat == 0) {
		nsec = 1;
		setsupr(&general->nx, &general->xmn, &general->xsiz, &general->ny,
				&general->ymn, &general->ysiz, &general->nz, &general->zmn,
				&general->zsiz, &general->nd, general->x, general->y, 
				general->z, general->vr, general->wt, &nsec, general->sec,
				&sec2, &sec3, &five, &five, &five, nisb, &nxsup, &xmnsup,
				&xsizsup, &nysup, &ymnsup, &ysizsup, &nzsup, &zmnsup, &zsizsup)
			;
		picksup(nxsup, xsizsup, nysup, ysizsup, nzsup, zsizsup,
				covariance->isrot, krige_vars->rotmat, search->radsqd,
				&nsbtosr, ixsbtosr, iysbtosr, izsbtosr);
	}
    */
	printf_dbg2("dssim() calling covtable\n");
	/* !Set up the covariance table and the spiral search: */
	covtable(&order[1], &sim[1], general, search, covariance, covtable_lookup, krige_vars);

	/* codigo reescrito ..................................... */

	/* prepare the random path */
	/* init the sim table with NOSVALUE (no simulated value) */
	for(i = 1; i <= general->nxyz; ++i) {
		order[i] = i;
		sim[i] = general->nosim_value;
	}

	toSim = general->nxyz;
	index = in = 0;

	/* copy wells data to simulation grid */
	for(i = 0; i < (int) general->wellsNPoints; i++) {
		in = general->wellsDataPos[i];
		/* there can be several wells datapoints for the same grid point */
		if( sim[in] != general->nosim_value )
			continue;
		sim[in] = (float) general->wellsDataVal[i];

		// mark point has simulated
		order[in] = order[toSim];
		order[toSim] = in;
		toSim--;
	}

	printf_dbg2("wellsNPoints: %d\tgeneral->nd: %d\n",general->wellsNPoints,general->nd);

	/* codigo reescrito ..................................... FIM */

	printf_dbg("dssim(): grid Points: %d\twells Points: %d\t toSim Points: %d\t should be: %d\n",
			general->nxyz, general->nd, toSim, general->nxyz - general->nd);
	printf_dbg2("\tdssim(): Starting simulation now\n");
	/* !MAIN LOOP OVER ALL THE NODES: */
	ierr = 0;
	in = 0;
	zmean = 0.f;
	zvariance = 0.f;
	while( toSim > 0 ) {
		
		if(toSim == (general->nxyz * 0.75)) {
			printf_dbg("\tdsssim(): 1/4 completed.\n");
		} else if(toSim == (general->nxyz / 2)) {
			printf_dbg("\tdsssim(): 1/2 completed.\n");
		} else if(toSim == (general->nxyz / 4)) {
			printf_dbg("\tdsssim(): 3/4 completed.\n");
 		}

		/* generate point to simulate */
		in = ((int) tsi_random() % toSim) +1; /* +1 because of fortran offsets */
		index = order[in];  /* point to be simulated */

		/* mark has simulated */
		order[in] = order[toSim];
		order[toSim] = index;
		toSim--;
		

		if(index < 1)
			printf("dssim(): ERROR, INDEX(%d) < 1\n", index);
		/* if value has a value allready (like in the case of a hard data guiven point), skip simulation */
		/* if( sim[index] >= general->tmin &&
		    sim[index] <= general->tmax) { */
		/*
		 if (sim[index] > general->nosim_value + 1e-20f || 
				 sim[index] < general->nosim_value * 2.f) {
			++ierr;
			printf("dssim(): ERROR: index(%d) allready has valid data: %f\n",index,sim[index]);
			continue;
		}
		*/

		/* get relative x,y,z from index */
		get3Dcoords(index, general->nx, general->nxy, &ix, &iy, &iz);

		/* getting absolute coords from relative coords */
		xx = general->xmn + (float) (ix - 1) * general->xsiz;
		yy = general->ymn + (float) (iy - 1) * general->ysiz;
		zz = general->zmn + (float) (iz - 1) * general->zsiz;

		/* !Now, we'll simulate the point ix,iy,iz.  First, get the close data */
		/* !and make sure that there are enough to actually simulate a value, */
		/* !we'll only keep the closest "ndmax" data, and look for previously */
		/* !simulated grid nodes: */
        /*
         * TSI NODE: SSTRAT IS ALLWAYS 1
		if (search->sstrat == 0) {
			srchsupr(xx, yy, zz, search->radsqd, covariance->isrot,
					krige_vars->rotmat, nsbtosr, ixsbtosr, iysbtosr, 
					izsbtosr, search->noct, general->x, 
					general->y, general->z, general->wt, nisb,
					nxsup, xmnsup, xsizsup,
					nysup, ymnsup, ysizsup,
					nzsup, zmnsup, zsizsup,
					&search->nclose, general->close);
			if (search->nclose < search->ndmin) {
				ierr++;
				printf("dssim(): SKIP2  %d\n",index);
				continue;
			}
			if (search->nclose > search->ndmax) {
				search->nclose = search->ndmax;
			}
		}
        */
		srchnod(ix, iy, iz, &sim[1], general, search, covtable_lookup);
		/* !WARNING:Para ter em atencao; bai c/ NOSIMVALUE, do simple kriging*/
		kinicial = ktype;
		if (ktype == 5 && bestAICube[index] == general->nosim_value) {
			ktype = 0;
		}
		/* !Calculate the conditional mean and standard deviation.  This will be */
		/* !done with kriging if there are data, otherwise, the global mean and */
		/* !standard deviation will be used: */
		if (ktype == 5) {
			global_mean = bestAICube[index];
		} else {
			global_mean = simulation->vmedexp;
		}
		if (search->nclose + covtable_lookup->ncnode < 1) {
			cmean = simulation->vmedexp;
			std_deviation = sqrt(simulation->vvarexp);
		} else {
			/* !Perform the kriging.  Note that if there are fewer than four data */
			/* !then simple kriging is prefered so that the variance of the */
			/* !realization does not become artificially inflated: */
			lktype = ktype;
			if (	ktype == 1 && 
					search->nclose + covtable_lookup->ncnode < 4) {
				lktype = 0;
			}
			/* !Estimacao em xo (SDSIM) */
			/* aceder ao bestCorrCube, apenas se ktype >= 4 */
			if(ktype == 5) {
				clcorr = bestCorrCube[index];
			}
				
			krige(ix, iy, iz, xx, yy, zz, lktype, global_mean, 
					&cmean, &std_deviation, // these are the output vars of krige
					&bestAICube[1], clcorr,
					general, search, simulation,
					covariance, covtable_lookup, krige_vars);

			/* this was used when dss did more than one simulation */
			/*
			if (simulation->nsim > 0) {
				if (covtable_lookup->icmean == 1) {
					cmean -= zmean / simulation->nsim - simulation->vmedexp;
				}
				// Computing 2nd power
				d__1 = zmean / simulation->nsim;
				cpdev = zvariance / simulation->nsim - d__1 * d__1;
				if (cpdev > 0.f) {
					cpdev = sqrt(cpdev);
					if (covtable_lookup->icvar == 1) {
						std_deviation = std_deviation * sqrt(simulation->vvarexp) / cpdev;
					}
				}
			}
			*/
		}
		ktype = kinicial;
		vmy = compute_gaussian_equiv( cmean, general->nd, general->vrtr, general->vrgtr);

		/* !Gera um valor aleatorio com distribuicao Gaussiana */
		vms = 0;
		for (i = 1; i <= covtable_lookup->ntry; ++i) {
			p = tsi_random_real();
			gauinv(&p, &xp); // XXX: not checking return code
			xp = xp * std_deviation + vmy;
			/* !Transformada inversa final (iv)               (SDSIM) */
			simval = backtr(xp, general->ntr, general->vrtr, general->vrgtr,
					general->min_value, general->max_value);
			vms += simval;
		}
		vms /= covtable_lookup->ntry;
		/* !Reter o valor simulado  (SDSIM) */
		if (covtable_lookup->ntry > 1) {
			sim[index] = simval + (cmean - vms);
		} else {
			sim[index] = simval;
		}
		if (sim[index] < general->min_value) {
			sim[index] = general->min_value;
		} else if (sim[index] > general->max_value) {
			sim[index] = general->max_value;
		}
		/* !Condicionamento as medias locais */
		/* used only when doing more than 1 simulation */
		/*
		zmean += sim[index]; 
		r__1 = sim[index];
		zvariance += r__1 * r__1;
		*/
	
	}	/* !END MAIN LOOP OVER NODES: */

	printf_dbg("dssim(): DEBUG: SKIP points: %d\n",ierr);

	ierr = 0;
	for(i = 0; i < general->wellsNPoints; i++) {
		in = general->wellsDataPos[i];
		simval = sim[in] - general->wellsDataVal[i];
		if(simval != 0) {
			printf_dbg("sim[%d] - wellsData = %f\n",in,simval);
			ierr++;
		}
	}	
	printf_dbg("dssim(): sim grid disrespects %d wells data points\n",ierr);

	/* !Return to the main program: */
	return 0;
} /* sdsim_ */

