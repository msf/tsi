/* ----------------------------------------------------------------------- */
/*           Conditional Simulation of a 3-D Rectangular Grid */

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
#include <math.h>

#include "dss.h"
#include "dss_legacy.h"
//#include "acorni.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)

#define TRUE_ (1)
#define FALSE_ (0)

/* Table of constant values */


int sdsim(float *sim, float *bestAICube, float *bestCorrCube, int *order, int *mask_data,
		general_vars_t			*	general,
		search_vars_t			*	search,
		simulation_vars_t		*	simulation,
		covariance_vars_t		*	covariance,
		covtable_lookup_vars_t	*	covtable_lookup,
		krige_vars_t			*	krige_vars)
{
	/* System generated locals */
	int i__1, i__2, i__3, i__4, i__5;
	float r__1, r__2, r__3;
	double d__1;


	/* Local variables */
	float xsizsup, ysizsup, zsizsup;
	int kinicial;
	int ixsbtosr[1000], iysbtosr[1000], izsbtosr[1000];
	float c__, d__, e, f, g, h__;
	int j;
	double p, zvariance;
	int testecomp;
	int id, in;
	int is, ix, jx, jy, jz, iz, iy;
	float xp, xx, yy, zz;
	int id2, ind, kkk, nnz, nny, nnx;
	float vmy, vms, sec2, sec3;
	int nsec, nisb[125], isim, ierr;
	float test, tiny, test2, cmean, gmean;

	float cpdev;
	double zmean;
	int index, irepo, imult, nvsim, nxsup, nysup, nzsup;
	int infoct;
	float cstdev;
	int lktype;
	float xmnsup, ymnsup, zmnsup, clcorr;
	float simval;
	int testind;
	int nsbtosr;

        int five = 5;
        int one = 1;

	/* Parameter adjustments */
	--mask_data;
	--order;
	--bestCorrCube;
	--bestAICube;
	--sim;

	/* Function Body */
	i__1 = covariance->nst[0];
	for (is = 1; is <= i__1; ++is) {
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
	/* !MAIN LOOP OVER ALL THE SIMULATIONS: */
	/*   printf("Step4 - out of covtable\n"); */
	i__1 = simulation->nsim;
	for (isim = 1; isim <= i__1; ++isim) {
		/* !Work out a random path for this realization: */
		i__2 = general->nxyz;
		for (ind = 1; ind <= i__2; ++ind) {
			sim[ind] = (float) acorni(general->ixv);
			order[ind] = ind;
		}
		/* ! 33 + 1 = SEED_POS */
		/* !ixv(1) = PARAMS(34) */
		/* !The multiple grid search works with multiples of 4 (yes, that is */
		/* !somewhat arbitrary): */

		/* TODO: 
		 * documentar: o que e' que isto esta' a fazer ???
		 */
		if (search->mults == 1) {
			i__2 = search->nmult;
			for (imult = 1; imult <= i__2; ++imult) {
				i__3 = 1;
			   	i__4 = general->nz / (imult << 2);
				nnz = MAX(i__3,i__4);
			   	i__4 = general->ny / (imult << 2);
				nny = MAX(i__3,i__4);
			   	i__4 = general->nx / (imult << 2);
				nnx = MAX(i__3,i__4);
				jz = 1;
				jy = 1;
				jx = 1;
				i__3 = nnz;
				for (iz = 1; iz <= i__3; ++iz) {
					if (nnz > 1) {
						jz = iz * imult << 2;
					}
					i__4 = nny;
					for (iy = 1; iy <= i__4; ++iy) {
						if (nny > 1) {
							jy = iy * imult << 2;
						}
						i__5 = nnx;
						for (ix = 1; ix <= i__5; ++ix) {
							if (nnx > 1) {
								jx = ix * imult << 2;
							}
							index = jx + (jy - 1) * general->nx + (jz - 1) * 
								general->nxy;
							sim[index] += imult;
						}
					}
				}
			}
		}

		/*  printf("Step6 - calling sortem\n"); */
		sortemi(&one, &general->nxyz, &sim[1], &one, &order[1], &c__, &d__, &e, &f, &g, &h__);
		/* !Initialize the simulation: */
		/*   printf("Step6 - out of sortem\n"); */
		nvsim = 0;
		i__2 = general->nxyz;
		for (ind = 1; ind <= i__2; ++ind) {
			sim[ind] = general->nosvalue;
		}
		/* !Assign the data to the closest grid node: */
		tiny = 1e-4f;
		testecomp = FALSE_;
		i__2 = general->nd;
		/*   printf("Step7\n"); */
		for (id = 1; id <= i__2; ++id) {
			getindx(&general->nx, &general->xmn, &general->xsiz, &general->x[id - 1], &ix, &testind);
			getindx(&general->ny, &general->ymn, &general->ysiz, &general->y[id - 1], &iy, &testind);
			getindx(&general->nz, &general->zmn, &general->zsiz, &general->z__[id - 1], &iz, &testind);
/*			GETINDX(&general->nx, &general->xmn, &general->xsiz, &general->x[id - 1], &ix);
			GETINDX(&general->ny, &general->ymn, &general->ysiz, &general->y[id - 1], &iy);
			GETINDX(&general->nz, &general->zmn, &general->zsiz, &general->z__[id - 1], &iz);
*/			ind = ix + (iy - 1) * general->nx + (iz - 1) * general->nxy;
			xx = general->xmn + (float) (ix - 1) * general->xsiz;
			yy = general->ymn + (float) (iy - 1) * general->ysiz;
			zz = general->zmn + (float) (iz - 1) * general->zsiz;
			test = (r__1 = xx - general->x[id - 1], (double) fabs(r__1)) + (r__2 = yy 
					- general->y[id - 1], (double) fabs(r__2)) + (r__3 = zz - 
					general->z__[id - 1], (double) fabs(r__3));
			/* !Assign this data to the node (unless there is a closer data): */
			if (search->sstrat == 1) {
				if (sim[ind] >= 0.f) {
					id2 = (int) (sim[ind] + .5f);
					test2 = (r__1 = xx - general->x[id2 - 1], (double) fabs(r__1)) + (
							r__2 = yy - general->y[id2 - 1], (double) fabs(r__2)) + (
							r__3 = zz - general->z__[id2 - 1], (double) fabs(r__3));
					if (test <= test2) {
						sim[ind] = (float) id;
					}
				} else {
					sim[ind] = (float) id;
				}
			}
			/* !Assign a flag so that this node does not get simulated: */
			if (search->sstrat == 0 && test <= tiny) {
				sim[ind] = general->nosvalue * 10.f;
			}
		}
		/* !Now, enter data values into the simulated grid: */
		/*   printf("Step8\n"); */
		i__2 = general->nxyz;
		for (ind = 1; ind <= i__2; ++ind) {
			id = (int) (sim[ind] + .5f);
			if (id > 0) {
				sim[ind] = general->vr[id - 1];
			}
		}
		i__4 = general->nxyz / 10;
		i__2 = 1;
		i__3 = MIN(i__4,1000); 
		irepo = MAX(i__2,i__3);
		/* !Assign a flag so that the node out of the mask boundaries does not get simulated: */
		if (general->imask == 1) {
			i__2 = general->nxyz;
			for (ind = 1; ind <= i__2; ++ind) {
				if (mask_data[ind] == 0) {
					sim[ind] = general->nosvalue;
				}
			}
		}
		/* !MAIN LOOP OVER ALL THE NODES: */
		simulation->nsim = 0;
		zmean = 0.f;
		zvariance = 0.f;
		i__2 = general->nxyz;
		for (in = 1; in <= i__2; ++in) {
			/* !Figure out the location of this point and make sure it has */
			/* !not been assigned a value already: */
			index = (int) (order[in] + .5f);
			if (sim[index] > general->nosvalue + 1e-20f || sim[index] < 
					general->nosvalue * 2.f) {
				goto L5;
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
					goto L5;
				}
				if (search->nclose > search->ndmax) {
					search->nclose = search->ndmax;
				}
			}
			srchnod(&ix, &iy, &iz, &sim[1],
					general,
					search,
					covtable_lookup);
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
				if (general->ktype == 1 && search->nclose + covtable_lookup->ncnode <
						4.f) {
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
					goto L9999;
				}
			}
L9999:
			/* !Gera um valor aleatorio com distribuicao Gaussiana */
			vms = 0.f;
			i__3 = covtable_lookup->ntry;
			for (kkk = 1; kkk <= i__3; ++kkk) {
				p = acorni(general->ixv);
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
			if (sim[index] <= general->tmin) {
				sim[index] = general->tmin;
			}
			if (sim[index] >= general->tmax) {
				sim[index] = general->tmax;
			}
			/* !Condicionamento as medias locais */
			++simulation->nsim;
			zmean += sim[index];
			/* Computing 2nd power */
			r__1 = sim[index];
			zvariance += r__1 * r__1;
			/* !END MAIN LOOP OVER NODES: */
L5:
			;
		}
		/* !Back transform each value and write results: */
		/*
		ne = 0;
		av = 0.f;
		ss = 0.f;
		r__1 = (float) ne;
		av /= (double) MAX(r__1,1.f);
		ss = ss / ((double) MAX(r__1,1.f)) - av * av;
		*/
		/* !END MAIN LOOP OVER SIMULATIONS: */
	}

	/* !Return to the main program: */
	return 0;
} /* sdsim_ */

