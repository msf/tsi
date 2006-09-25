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

static int five = 5;

int dssim(float *sim, float *bestAICube, float *bestCorrCube, int *order, int *mask_data,
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
	float xsizsup, ysizsup, zsizsup;
	int kinicial;
	int ixsbtosr[1000], iysbtosr[1000], izsbtosr[1000];
	double p, zvariance;
	int in;
	int ix, iz, iy;
	float xp, xx, yy, zz;
	float vmy = 0, vms, sec2, sec3;
	int nsec, nisb[125], ierr;
	float cmean, gmean;

	double zmean;
	int index, nxsup, nysup, nzsup;
	float cstdev;
	int lktype;
	float xmnsup, ymnsup, zmnsup, clcorr;
	float simval = 0;
	int nsbtosr;


	/* Parameter adjustments */
	--mask_data;
	--order;
	--bestCorrCube;
	--bestAICube;
	--sim;

	printf_dbg2("\tdssim() called\n");

	/* Function Body */
	for (i = 1; i <= covariance->nst[0]; ++i) {
		setrot(covariance->ang1[i - 1], covariance->ang2[i - 1], covariance->ang3[i - 1],
				covariance->anis1[i - 1], covariance->anis2[i - 1], i, krige_vars->rotmat);
	}
	covariance->isrot = 5;
	setrot(search->sang1, search->sang2, search->sang3,
			search->sanis1, search->sanis2, covariance->isrot, 
			krige_vars->rotmat);
	/* !Set up the super block search: */
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

	printf_dbg2("dssim() calling covtable\n");
	/* !Set up the covariance table and the spiral search: */
	covtable(&order[1], &sim[1], general, search, covariance, covtable_lookup, krige_vars);

	/* codigo reescrito ..................................... */

	/* prepare the random path */
	/* init the sim table with NOSVALUE (no simulated value) */
	for(i = 1; i <= general->nxyz; ++i) {
		order[i] = i;
		sim[i] = general->nosvalue;
	}

	toSim = general->nxyz;
	index = in = 0;

	/* copy wells data to simulation grid */
	for(i = 0; i < (int) general->wellsNPoints; i++) {
		in = general->wellsDataPos[i];
		sim[in] = (float) general->wellsDataVal[i];

		// mark point has simulated
		order[in] = order[toSim];
		order[toSim] = in;
		toSim--;
	}

	printf_dbg2("wellsNPoints: %d\tgeneral->nd: %d\n",general->wellsNPoints,general->nd);

	/* codigo reescrito ..................................... FIM */

/* !Assign a flag so that the node out of the mask boundaries does not get simulated: */
	if (general->imask == 1) {
		for (i = 1; i <= general->nxyz; ++i) {
			if (mask_data[i] == 0) {
				sim[i] = general->nosvalue;
			}
		}
	}

	printf_dbg("dssim(): grid Points: %d\twells Points: %d\t toSim Points: %d\t should be: %d\n",
			general->nxyz, general->nd, toSim, general->nxyz - general->nd);
	printf_dbg2("\tdssim(): Starting simulation now\n");
	/* !MAIN LOOP OVER ALL THE NODES: */
	ierr = 0;
	in = 0;
	zmean = 0.f;
	zvariance = 0.f;
	while( toSim > 0 ) {
		
		if(toSim == (general->nxyz * 0.75))
			printf_dbg("\tdsssim(): 1/4 completed.\n");
		else if(toSim == (general->nxyz / 2))
			printf_dbg("\tdsssim(): 1/2 completed.\n");
		else if(toSim == (general->nxyz / 4))
			printf_dbg("\tdsssim(): 3/4 completed.\n");

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
		 if (sim[index] > general->nosvalue + 1e-20f || 
				 sim[index] < general->nosvalue * 2.f) {
			++ierr;
			printf("dssim(): ERROR: index(%d) allready has valid data: %f\n",index,sim[index]);
			continue;
		}
		*/

		/* get relative x,y,z from index */
		iz = (index - 1) / general->nxy + 1;
		iy = (index - (iz - 1) * general->nxy - 1) / general->nx + 1;
		ix = index - (iz - 1) * general->nxy - (iy - 1) * general->nx;

		/* getting absolute coords from relative coords */
		xx = general->xmn + (float) (ix - 1) * general->xsiz;
		yy = general->ymn + (float) (iy - 1) * general->ysiz;
		zz = general->zmn + (float) (iz - 1) * general->zsiz;

		/* !Now, we'll simulate the point ix,iy,iz.  First, get the close data */
		/* !and make sure that there are enough to actually simulate a value, */
		/* !we'll only keep the closest "ndmax" data, and look for previously */
		/* !simulated grid nodes: */
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
		srchnod(ix, iy, iz, &sim[1], general, search, covtable_lookup);
		/* !WARNING:Para ter em atencao; bai c/ NOSIMVALUE, do simple kriging*/
		kinicial = general->ktype;
		if (general->ktype == 5 && bestAICube[index] == general->nosvalue) {
			general->ktype = 0;
		}
		/* !Calculate the conditional mean and standard deviation.  This will be */
		/* !done with kriging if there are data, otherwise, the global mean and */
		/* !standard deviation will be used: */
		if (general->ktype == 2 || general->ktype >= 4) {
			gmean = bestAICube[index];
		} else {
			gmean = simulation->vmedexp;
		}
		if (search->nclose + covtable_lookup->ncnode < 1) {
			cmean = simulation->vmedexp;
			cstdev = sqrt(simulation->vvarexp);
		} else {
			/* !Perform the kriging.  Note that if there are fewer than four data */
			/* !then simple kriging is prefered so that the variance of the */
			/* !realization does not become artificially inflated: */
			lktype = general->ktype;
			if (	general->ktype == 1 && 
					search->nclose + covtable_lookup->ncnode < 4) {
				lktype = 0;
			}
			/* !Estimacao em xo (SDSIM) */
			/* aceder ao bestCorrCube, apenas se ktype >= 4 */
			if(general->ktype == 5) {
				clcorr = bestCorrCube[ix + (iy - 1) * general->nx + (iz - 1) * general->nxy];
			}
				
			krige(ix, iy, iz, xx, yy, zz, lktype, gmean, 
					&cmean, &cstdev, // these are the output vars of krige
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
						cstdev = cstdev * sqrt(simulation->vvarexp) / cpdev;
					}
				}
			}
			*/
		}
		general->ktype = kinicial;
		/* !Calculo do equivalente valor gaussiano        (SDSIM) */
		if (cmean <= general->vrtr[0]) {
			vmy = general->vrgtr[0];
			goto L9999;
		}
		if (cmean >= general->vrtr[general->nd - 1]) {
			vmy = general->vrgtr[general->nd - 1];
			goto L9999;
		}
		for (i = 1; i < general->nd; ++i) {
			if (cmean >= general->vrtr[i - 1] && cmean < general->vrtr[i]) {
				vmy = general->vrgtr[i - 1] + (cmean - general->vrtr[i - 1]) * 
					(general->vrgtr[i] - general->vrgtr[i - 1]) / (general->vrtr[i] - general->vrtr[i - 1]);
				break;
			}
		}
L9999:
		/* !Gera um valor aleatorio com distribuicao Gaussiana */
		vms = 0;
		for (i = 1; i <= covtable_lookup->ntry; ++i) {
			p = tsi_random_real();
			gauinv(&p, &xp, &ierr);
			xp = xp * cstdev + vmy;
			/* !Transformada inversa final (iv)               (SDSIM) */
			simval = backtr(xp, general->ntr, general->vrtr, general->vrgtr,
					general->zmin, general->zmax, general->ltail, general->ltpar,
					general->utail, general->utpar);
			vms += simval;
		}
		vms /= covtable_lookup->ntry;
		/* !Reter o valor simulado  (SDSIM) */
		if (covtable_lookup->ntry > 1) {
			sim[index] = simval + (cmean - vms);
		} else {
			sim[index] = simval;
		}
		if (sim[index] < general->tmin) {
			sim[index] = general->tmin;
		} else if (sim[index] > general->tmax) {
			sim[index] = general->tmax;
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

