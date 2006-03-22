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
/*   cova3            Calculates the covariance given a variogram model */
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

#include "dss.h"
#include "acorni.h"
#include "profile.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)

#define TRUE_ (1)
#define FALSE_ (0)

/* Table of constant values */

static int five = 5;

int sdsim(float *sim, float *bestAICube, float *bestCorrCube, int *order, int *mask_data,
		general_vars_t			*	general,
		search_vars_t			*	search,
		simulation_vars_t		*	simulation,
		covariance_vars_t		*	covariance,
		covtable_lookup_vars_t	*	covtable_lookup,
		krige_vars_t			*	krige_vars)
{
	/* System generated locals */
	int i__3;
	float r__1;
	double d__1;


	/* Local variables */
	float xsizsup, ysizsup, zsizsup;
	int kinicial;
	int ixsbtosr[1000], iysbtosr[1000], izsbtosr[1000];
	int j;
	double p, zvariance;
	int in;
	int is, ix, iz, iy;
	float xp, xx, yy, zz;
	int ind, kkk;
	float vmy = 0, vms, sec2, sec3;
	int nsec, nisb[125], ierr;
	float cmean, gmean;

	float cpdev;
	double zmean;
	int index, nxsup, nysup, nzsup;
	int infoct;
	float cstdev;
	int lktype;
	float xmnsup, ymnsup, zmnsup, clcorr;
	float simval = 0;;
	int nsbtosr;

#ifdef PROFILE
	profile.sdsim++;
	profBegin("sdsim");
#endif

	/* Parameter adjustments */
	--mask_data;
	--order;
	--bestCorrCube;
	--bestAICube;
	--sim;

	/* Function Body */
	for (is = 1; is <= covariance->nst[0]; ++is) {
		setrot(&covariance->ang1[is - 1], &covariance->ang2[is - 1],
				&covariance->ang3[is - 1], &covariance->anis1[is - 1],
				&covariance->anis2[is - 1], &is, &five, krige_vars->rotmat);
	}
	covariance->isrot = 5;
	setrot(&search->sang1, &search->sang2, &search->sang3,
			&search->sanis1, &search->sanis2, &covariance->isrot, &five, 
			krige_vars->rotmat);
	/* !Set up the super block search: */
	if (search->sstrat == 0) {
		nsec = 1;
		setsupr(&general->nx, &general->xmn, &general->xsiz, &general->ny,
				&general->ymn, &general->ysiz, &general->nz, &general->zmn,
				&general->zsiz, &general->nd, general->x, general->y, 
				general->z__, general->vr, general->wt, &nsec, general->sec,
				&sec2, &sec3, &five, &five, &five, nisb, &nxsup, &xmnsup,
				&xsizsup, &nysup, &ymnsup, &ysizsup, &nzsup, &zmnsup, &zsizsup)
			;
		picksup(&nxsup, &xsizsup, &nysup, &ysizsup, &nzsup, &zsizsup,
				&covariance->isrot, &five, krige_vars->rotmat, &search->radsqd,
				&nsbtosr, ixsbtosr, iysbtosr, izsbtosr);
	}
	/* !Set up the covariance table and the spiral search: */
	/*   printf("Step3 - calling covtable\n"); */
	covtable(&order[1], &sim[1], general, search, covariance, covtable_lookup, krige_vars);
	/* !nxyz2 = nxyz / 2 */

	/* codigo reescrito ..................................... */

	/* create a random path */
	for(ind = 1; ind <= general->nxyz; ++ind) {
		sim[ind] = random() % general->nxyz;
		order[ind] = ind;
	}
	/* permuta o array order pela ordem do caminho indicado em sim[] */
	sort_permute_int(0, general->nxyz -1, &sim[1], &order[1]);
	/* agora temos em order[] as posicoes a serem simuladas */

	/* iniciar a grid de simulação a NOSVALUES (no simulated value) */
	for (ind = 1; ind <= general->nxyz; ++ind) {
		sim[ind] = general->nosvalue;
	}

	/* introduzir os dados dos poc,os (wells/harddata) na grid */
	for(ind = 0; ind < general->wellsNPoints; ind++) {
		sim[general->wellsDataPos[ind]] = (float) general->wellsDataVal[ind];
	}

	/* codigo reescrito ..................................... FIM */
/* !Assign a flag so that the node out of the mask boundaries does not get simulated: */
	if (general->imask == 1) {
		for (ind = 1; ind <= general->nxyz; ++ind) {
			if (mask_data[ind] == 0) {
				sim[ind] = general->nosvalue;
			}
		}
	}
	/* !MAIN LOOP OVER ALL THE NODES: */
	simulation->nsim = index = 0;
	zmean = 0.f;
	zvariance = 0.f;
	for (in = 1; in <= general->nxyz; ++in) {
		/* verificar se a ultima simulacao foi bem sucedida */
		if( in > 1 && sim[index] < 0)
				printf("dssim(): ERRO: ponto %d nao foi simulado!!\n",index);

		index =  order[in];
		/* if value has a value allready (like in the case of a hard data guiven point), skip simulation */
		/* if( sim[index] >= general->tmin &&
		    sim[index] <= general->tmax) { */
		if( sim[index] != general->nosvalue) {
			/* value was allready  simulated! */
			++simulation->nsim;
			continue;
		}

		iz = (index - 1) / general->nxy + 1;
		iy = (index - (iz - 1) * general->nxy - 1) / general->nx + 1;
		ix = index - (iz - 1) * general->nxy - (iy - 1) * general->nx;

		xx = general->xmn + (float) (ix - 1) * general->xsiz;
		yy = general->ymn + (float) (iy - 1) * general->ysiz;
		zz = general->zmn + (float) (iz - 1) * general->zsiz;

		/* !Now, we'll simulate the point ix,iy,iz.  First, get the close data */
		/* !and make sure that there are enough to actually simulate a value, */
		/* !we'll only keep the closest "ndmax" data, and look for previously */
		/* !simulated grid nodes: */
		if (search->sstrat == 0) {
			srchsupr(&xx, &yy, &zz, &search->radsqd, &covariance->isrot,
					&five, krige_vars->rotmat, &nsbtosr, ixsbtosr, iysbtosr, 
					izsbtosr, &search->noct, &general->nd, general->x, 
					general->y, general->z__, general->wt, nisb, &nxsup, 
					&xmnsup, &xsizsup, &nysup, &ymnsup, &ysizsup, &nzsup, 
					&zmnsup, &zsizsup, &search->nclose, general->close, 
					&infoct);
			if (search->nclose < search->ndmin) {
				printf("dssim(): SKIP2  %d\n",index);
				continue;
			}
			if (search->nclose > search->ndmax) {
				search->nclose = search->ndmax;
			}
		}
		srchnod(&ix, &iy, &iz, &sim[1], general, search, covtable_lookup);
		/* !WARNING:Para ter em atencaoi; valores -999.25 a krigagem passa a simples */
		kinicial = general->ktype;
		if (general->ktype == 5 && bestAICube[index] == -999.25f) {
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
		if (search->nclose + covtable_lookup->ncnode < 1.f) {
			cmean = simulation->vmedexp;
			cstdev = sqrt(simulation->vvarexp);
		} else {
			/* !Perform the kriging.  Note that if there are fewer than four data */
			/* !then simple kriging is prefered so that the variance of the */
			/* !realization does not become artificially inflated: */
			lktype = general->ktype;
			if (	general->ktype == 1 && 
					search->nclose + covtable_lookup->ncnode < 4.f) {
				lktype = 0;
			}
			/* !Estimacao em xo (SDSIM) */
			/* aceder ao bestCorrCube, apenas se ktype >= 4 */
			if(general->ktype == 5)
				clcorr = bestCorrCube[ix + (iy - 1) * general->nx + (iz - 1) * general->nxy];

			krige(&ix, &iy, &iz, &xx, &yy, &zz, &lktype, &gmean, &cmean, 
					&cstdev, &bestAICube[1], &clcorr,
					general, search, simulation,
					covariance, covtable_lookup, krige_vars);

			if (simulation->nsim > 0) {
				if (covtable_lookup->icmean == 1) {
					cmean -= zmean / simulation->nsim - simulation->vmedexp;
				}
				/* Computing 2nd power */
				d__1 = zmean / simulation->nsim;
				cpdev = zvariance / simulation->nsim - d__1 * d__1;
				if (cpdev > 0.f) {
					cpdev = sqrt(cpdev);
					if (covtable_lookup->icvar == 1) {
						cstdev = cstdev * sqrt(simulation->vvarexp) / cpdev;
					}
				}
			}
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
		i__3 = general->nd - 1;
		for (j = 1; j <= i__3; ++j) {
			if (cmean >= general->vrtr[j - 1] && cmean < general->vrtr[j])
			{
				vmy = general->vrgtr[j - 1] + (cmean - general->vrtr[j - 
						1]) * (general->vrgtr[j] - general->vrgtr[j - 1]) 
					/ (general->vrtr[j] - general->vrtr[j - 1]);
				break;
			}
		}
L9999:
		/* !Gera um valor aleatorio com distribuicao Gaussiana */
		vms = 0;
		for (kkk = 1; kkk <= covtable_lookup->ntry; ++kkk) {
			p = ((double)random()) / RAND_MAX;
			gauinv(&p, &xp, &ierr);
			xp = xp * cstdev + vmy;
			/* !Transformada inversa final (iv)               (SDSIM) */
			simval = backtr(&xp, &general->ntr, general->vrtr, 
					general->vrgtr, &general->zmin, &general->zmax, 
					&general->ltail, &general->ltpar, &general->utail, 
					&general->utpar);
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
		zmean += sim[index];
		/* Computing 2nd power */
		r__1 = sim[index];
		zvariance += r__1 * r__1;
	
	}	/* !END MAIN LOOP OVER NODES: */

	printf("dssim(): DEBUG: SKIP points: %d\n",simulation->nsim);

#ifdef PROFILE
	profEnd("sdsim");
#endif

	/* !Return to the main program: */
	return 0;
} /* sdsim_ */

