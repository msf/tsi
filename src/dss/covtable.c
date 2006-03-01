#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


/* ----------------------------------------------------------------------- */
/*               Establish the Covariance Look up Table */
/* The idea is to establish a 3-D network that contains the covariance */
/* value for a range of grid node offsets that should be at as large */
/* as twice the search radius in each direction.  The reason it has to */
/* be twice as large as the search radius is because we want to use it */
/* to compute the data covariance matrix as well as the data-point */
/* covariance matrix. */

/* Secondly, we want to establish a search for nearby nodes that */
/* in order of closeness as defined by the variogram. */

/* INPUT VARIABLES: */
/*   xsiz,ysiz,zsiz  Definition of the grid being considered */
/*   MAXCTX,Y,Z      Number of blocks in covariance table */
/*   covariance table parameters */

/* OUTPUT VARIABLES:  covtab()         Covariance table */

/* EXTERNAL REFERENCES: */
/*   sqdist          Computes 3-D anisotropic squared distance */
/*   sortem          Sorts multiple arrays in ascending order */
/*   cova3           Computes the covariance according to a 3-D model */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 * sortem 
 * cova3
 * sqdist
 */

/** CUBOS utilizados
 * tmp - allocado e dealocado na funcao 
 * order
 */ 

/** CUBOS _nao_ utilizados
 * sim
 * lvm
 * restantes
 * -
 */

/** structs globais utilizadas:
 * general 
 * search
 * covariance (passado ao cova3 e sqdist)
 * covtable_lookup
 * krige_vars (passado ao cova3: krige_vars->rotmat & krige_vars->ccb)
 */




int covtable(int *order, float * tmp,
		general_vars_t * general, 
		search_vars_t * search, 
		covariance_vars_t * covariance, 
		covtable_lookup_vars_t * covtable_lookup, 
		krige_vars_t * krige_vars)
{
        /* Table of constant values */
        float c_b2 = 0.f;
        int one = 1;
        int four = 4;
        int five = 5;

	/* System generated locals */
	int i__1, i__2, i__3;

	/* Local variables */
	/*	float * tmp; */
	float c__, d__, e, f, g, h__;
	int i__, j, k, ic, jc, kc, il, ix, iy, iz;
	float xx, yy, zz;
	int loc;
	double hsqd;

#ifdef PROFILE
	profile.ctable++;
	profBegin("ctable");
#endif

	/* Parameter adjustments */
	--order;

	/* Function Body */
	i__1 = (general->nx-1)/2;
	covtable_lookup->nctx = i__1; /* min(i__1,i__2);*/
	i__1 = (general->ny-1)/2;
	covtable_lookup->ncty = i__1; /* min(i__1,i__2); */
	i__1 = (general->nz-1)/2;
	covtable_lookup->nctz = i__1; /* min(i__1,i__2); */

	/* NOTE: If dynamically allocating memory, and if there is no shortage */
	/*       it would a good idea to go at least as far as the radius and */
	/*       twice that far if you wanted to be sure that all covariances */
	/*       in the left hand covariance matrix are within the table look-up. */

	/* Initialize the covariance subroutine and cbb at the same time: */
	/*    printf("Calling cova3\n"); */
	cova3(&c_b2, &c_b2, &c_b2, &c_b2, &c_b2, &c_b2, &one, covariance->nst, &four,
			covariance->c0, covariance->it, covariance->cc, covariance->aa, &one, 
			&five, krige_vars->rotmat, &covariance->cmax, &krige_vars->cbb);
	/* 		Now, set up the table and keep track of the node offsets that are */
	/* 		within the search radius: */
	/*    printf("loop 1/3\n"); */
 
	/* covtable needs a tmp array the size of aicube 
	tmp = (float *) malloc(general->nxyz * sizeof(float));
	*/

	covtable_lookup->nlooku = 0;
	i__1 = covtable_lookup->nctx;
	for (i__ = -covtable_lookup->nctx; i__ <= i__1; ++i__) {
		xx = i__ * general->xsiz;
		ic = covtable_lookup->nctx + 1 + i__;
		i__2 = covtable_lookup->ncty;
		for (j = -covtable_lookup->ncty; j <= i__2; ++j) {
			yy = j * general->ysiz;
			jc = covtable_lookup->ncty + 1 + j;
			i__3 = covtable_lookup->nctz;
			for (k = -covtable_lookup->nctz; k <= i__3; ++k) {
				zz = k * general->zsiz;
				kc = covtable_lookup->nctz + 1 + k;

				cova3(&c_b2, &c_b2, &c_b2, &xx, &yy, &zz, &one, 
						covariance->nst, &four, covariance->c0, covariance->it, 
						covariance->cc, covariance->aa, &one, &five, 
						krige_vars->rotmat, &covariance->cmax, 
						&covtable_lookup->covtab[getPos(ic,jc,kc, general->nx, general->nxy)]);

				hsqd = sqdist(&c_b2, &c_b2, &c_b2, &xx, &yy, &zz, &
						covariance->isrot, &five, krige_vars->rotmat);

				if ((float) hsqd <= search->radsqd) {
					++covtable_lookup->nlooku;
					/*		    printf("nlooku: %i\n", covtable_lookup->nlooku); 
					 						We want to search by closest variogram distance (and use the 
					 						anisotropic Euclidean distance to break ties: */
					tmp[covtable_lookup->nlooku] = -(covtable_lookup->covtab[getPos(ic,jc,kc, general->nx, general->nxy)] - (float) hsqd * 1e-10f);
					order[covtable_lookup->nlooku] = getPos(ic,jc,kc, general->nx, general->nxy);
				}

			}
		}
	}

	/* 		Finished setting up the look-up table, now order the nodes such */
	/* 		that the closest ones, according to variogram distance, are searched */
	/* 		first. Note: the "loc" array is used because I didn't want to make */
	/* 		special allowance for 2 byte integers in the sorting subroutine: */
	/*    printf("calling sortem\n"); */

	sortemi(&one, &covtable_lookup->nlooku, &tmp[1], &one, &order[1],
               &c__, &d__, &e, &f, &g, &h__);
	/*    printf("loop 2/3\n"); */
 	
 	/* tmp array is no longer needed. 
 	free(tmp);
 	*/

	i__1 = covtable_lookup->nlooku;  
	for (il = 1; il <= i__1; ++il) {
		loc = order[il];

		iz = (loc - 1) / (general->nxy) + 1;
		iy = (loc - (iz - 1) * (general->nxy) - 1) / general->nx + 1;
		ix = loc - (iz - 1) * (general->nxy) - (iy - 1) * general->nx;

		covtable_lookup->iznode[il - 1] = iz;
		covtable_lookup->iynode[il - 1] = iy;
		covtable_lookup->ixnode[il - 1] = ix;
		/*      printf("loop 2/3: %i %i %i\n", il, (int)covtable_lookup->nlooku, general->nxyz); */

	}

	if (covtable_lookup->nodmax > 64) {
		fprintf(stderr, "covtable(): covtable_lookup->nodmax above limit of 64, setting to 64!");
		/* apagar */
		/*            write(ldbg,*) */
		/*            write(ldbg,*) 'The maximum number of close nodes = ',nodmax */
		/*            write(ldbg,*) 'this was reset from your specification due ' */
		/*            write(ldbg,*) 'to storage limitations.' */
		covtable_lookup->nodmax = 64;
	}

#ifdef PROFILE
	profEnd("ctable");
#endif
	return 0;
} /* ctable */


