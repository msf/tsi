/* ----------------------------------------------------------------------- */
/*            Builds and Solves the SK or OK Kriging System */
/* INPUT VARIABLES: */
/*   ix,iy,iz        index of the point currently being simulated */
/*   xx,yy,zz        location of the point currently being simulated */
/* OUTPUT VARIABLES: */
/*   cmean           kriged estimate */
/*   cstdev          kriged standard deviation */
/* EXTERNAL REFERENCES: ksol   Gaussian elimination system solution */
/* ----------------------------------------------------------------------- */


/** Funcoes utilizadas
 * ksol
 * cova3
 */

/** CUBOS utilizados
 * bestAICube quando ktype == 5
 */ 

/** CUBOS _nao_ utilizados
 * sim
 * tmp
 * order
 * clc
 * e os cov*
 */

/** structs globais utilizadas:
 * generl_1
 * clooku_1
 * krigev_1
 * search_1
 * simula_1
 */

#include <math.h>
#include "dss.h"
#include "dss_legacy.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define FALSE (0)

/* Table of constant values */

static int one = 1;

int krige(int *ix, int *iy, int *iz, float *xx, float *yy, float *zz,
		int *lktype, float *gmean, float *cmean, float * cstdev, 
		float *bestAICube, float *clcorr,
		general_vars_t * general,
		search_vars_t * search,
		simulation_vars_t * simulation,
		covariance_vars_t * covariance,
		covtable_lookup_vars_t * covtable_lookup,
		krige_vars_t * krige_vars)
{
	/* System generated locals */
	int i1, i2;

	/* Local variables */
	int i, j;
	float x1, y1, z1, x2, y2, z2;
	int na, ii, jj, kk, in, ix1, iy1, iz1, ix2, iy2, iz2, ind, neq;
	float cov, vrea[129];
	float edmin, edmax;
	int index;
	float sfmin, sfmax;
	int ising;
	int first;
	float sumwts;


	/* Parameter adjustments */
	--bestAICube;

	/* init vars */
	ix1 = iy1 = iz1 = ix2 = iy2 = iz2 = 0;
	
	/* Function Body */
	first = FALSE;
	na = search->nclose + covtable_lookup->ncnode;
	neq = na;
	if (*lktype == 1) {
		neq += 1;
	}
	else if (*lktype == 3) {
		neq += 2;
	}
	else if (*lktype >= 4) {
		neq += 1;
		if (*lktype == 5) {
			general->colocorr = *clcorr;
		}
	}
	/* Set up kriging matrices: */
	in = 0;
	i1 = na;
	for (j = 1; j <= na; ++j) {
		/* Sort out the actual location of point "j" */
		if ((float) j <= search->nclose) {
			index = (int) general->close[j - 1];
			x1 = general->x[index - 1];
			y1 = general->y[index - 1];
			z1 = general->z[index - 1];
			krige_vars->vra[j - 1] = general->vr[index - 1];
			vrea[j - 1] = general->sec[index - 1];
		} else {
			/* It is a previously simulated node (keep index for table look-up): */
			index = j - search->nclose;
			x1 = covtable_lookup->cnodex[index - 1];
			y1 = covtable_lookup->cnodey[index - 1];
			z1 = covtable_lookup->cnodez[index - 1];
			krige_vars->vra[j - 1] = covtable_lookup->cnodev[index - 1];
			ind = covtable_lookup->icnode[index - 1];
			ix1 = *ix + (covtable_lookup->ixnode[ind - 1] - covtable_lookup->nctx - 1);
			iy1 = *iy + (covtable_lookup->iynode[ind - 1] - covtable_lookup->ncty - 1);
			iz1 = *iz + (covtable_lookup->iznode[ind - 1] - covtable_lookup->nctz - 1);
			index = ix1 + (iy1 - 1) * general->nx + (iz1 - 1) * general->nxy;
			if(general->ktype == 5)
				vrea[j - 1] = bestAICube[index];
		}
		if (*lktype == 0) {
			krige_vars->vra[j - 1] -= *gmean;
		}
		if (*lktype == 2) {
			krige_vars->vra[j - 1] -= vrea[j - 1];
		}
		if (*lktype >= 4) {
			krige_vars->vra[j - 1] -= simulation->vmedexp;
		}
		i2 = j;
		for (i = 1; i <= i2; ++i) {
			/* Sort out the actual location of point "i" */
			if ((float) i <= search->nclose) {
				index = (int) general->close[i - 1];
				x2 = general->x[index - 1];
				y2 = general->y[index - 1];
				z2 = general->z[index - 1];
			} else {
				/* It is a previously simulated node (keep index for table
				   look-up): */
				index = i - search->nclose;
				x2 = covtable_lookup->cnodex[index - 1];
				y2 = covtable_lookup->cnodey[index - 1];
				z2 = covtable_lookup->cnodez[index - 1];
				ind = covtable_lookup->icnode[index - 1];
				ix2 = *ix + (covtable_lookup->ixnode[ind - 1] - covtable_lookup->nctx - 1);
				iy2 = *iy + (covtable_lookup->iynode[ind - 1] - covtable_lookup->ncty - 1);
				iz2 = *iz + (covtable_lookup->iznode[ind - 1] - covtable_lookup->nctz - 1);
			}
			/* Now, get the covariance value: */
			++in;
			/* Decide whether or not to use the covariance look-up table: */
			if ((float) j <= search->nclose || (float) i <= search->nclose) 
			{
				cov = cova3(x1, y1, z1, x2, y2, z2, covariance->nst,
						covariance->c0, covariance->it, covariance->cc, 
						covariance->aa, krige_vars->rotmat, &
						covariance->cmax);
				krige_vars->a[in - 1] = (double) cov;
			} else {
				/* Try to use the covariance look-up (if the distance is in range): */
				ii = covtable_lookup->nctx + 1 + (ix1 - ix2);
				jj = covtable_lookup->ncty + 1 + (iy1 - iy2);
				kk = covtable_lookup->nctz + 1 + (iz1 - iz2);
				if (ii < 1 || ii > general->nx || jj < 1 || jj > general->ny || kk < 1 || kk 
						> general->nz) {
					cov = cova3(x1, y1, z1, x2, y2, z2, covariance->nst,
							covariance->c0, covariance->it, covariance->cc, 
							covariance->aa, krige_vars->rotmat, &covariance->cmax);
				} else {
					cov = covtable_lookup->covtab[getPos(ii,jj,kk, general->nx, general->nxy)];
				}
				krige_vars->a[in - 1] = (double) cov;
			}
		}
		/* Get the RHS value (possibly with covariance look-up table): */
		if ((float) j <= search->nclose) {
			cov = cova3(*xx, *yy, *zz, x1, y1, z1, covariance->nst,
					covariance->c0, covariance->it, covariance->cc, covariance->aa,
					krige_vars->rotmat, &covariance->cmax);
			krige_vars->r__[j - 1] = (double) cov;
		} else {
			/* Try to use the covariance look-up (if the distance is in range): */
			ii = covtable_lookup->nctx + 1 + (*ix - ix1);
			jj = covtable_lookup->ncty + 1 + (*iy - iy1);
			kk = covtable_lookup->nctz + 1 + (*iz - iz1);
			if (ii < 1 || ii > general->nx || jj < 1 || jj > general->ny || kk < 1 || kk > 
					general->nz) {
				cov = cova3(*xx, *yy, *zz, x1, y1, z1, covariance->nst,
						covariance->c0, covariance->it, covariance->cc, covariance->aa, 
						krige_vars->rotmat, &covariance->cmax);
			} else {
				cov = covtable_lookup->covtab[getPos(ii,jj,kk, general->nx, general->nxy)];
			}
			krige_vars->r__[j - 1] = (double) cov;
		}
		krige_vars->rr[j - 1] = krige_vars->r__[j - 1];
	}
	/* Addition of OK constraint: */
	if (*lktype == 1 || *lktype == 3) {
		i1 = na;
		for (i = 1; i <= i1; ++i) {
			++in;
			krige_vars->a[in - 1] = 1.f;
		}
		++in;
		krige_vars->a[in - 1] = 0.f;
		krige_vars->r__[na] = 1.f;
		krige_vars->rr[na] = 1.f;
	}
	/* Addition of the External Drift Constraint: */
	if (*lktype == 3) {
		edmin = 999999.f;
		edmax = -999999.f;
		i1 = na;
		for (i = 1; i <= i1; ++i) {
			++in;
			krige_vars->a[in - 1] = vrea[i - 1];
			if (krige_vars->a[in - 1] < edmin) {
				edmin = krige_vars->a[in - 1];
			}
			else if (krige_vars->a[in - 1] > edmax) {
				edmax = krige_vars->a[in - 1];
			}
		}
		++in;
		krige_vars->a[in - 1] = 0.f;
		++in;
		krige_vars->a[in - 1] = 0.f;
		ind = *ix + (*iy - 1) * general->nx + (*iz - 1) * general->nxy;
		krige_vars->r__[na + 1] = (double) bestAICube[ind];
		krige_vars->rr[na + 1] = krige_vars->r__[na + 1];
		if (edmax - edmin < 1e-20f) {
			--neq;
		}
	}
	/* Addition of Collocated Cosimulation Constraint: */
	else if (*lktype >= 4) {
		sfmin = 1e21f;
		sfmax = -1e21f;
		i1 = na;
		for (i = 1; i <= i1; ++i) {
			++in;
			krige_vars->a[in - 1] = (double) general->colocorr * 
				krige_vars->r__[i - 1];
			if (krige_vars->a[in - 1] < sfmin) {
				sfmin = krige_vars->a[in - 1];
			}
			else if (krige_vars->a[in - 1] > sfmax) {
				sfmax = krige_vars->a[in - 1];
			}
		}
		++in;
		krige_vars->a[in - 1] = 1.f;
		ii = na + 1;
		krige_vars->r__[ii - 1] = (double) general->colocorr;
		krige_vars->rr[ii - 1] = krige_vars->r__[ii - 1];
		/* apagar */
		/*           if((sfmax-sfmin).lt.EPSLON) neq = neq - 1 */
	}
	/* 		Write out the kriging Matrix if Seriously Debugging: */
	/* 		Solve the Kriging System: */
	if (neq == 1 && *lktype != 3) {
		krige_vars->s[0] = krige_vars->r__[0] / krige_vars->a[0];
		ising = 0;
	} else {
		ksol(&one, &neq, &one, krige_vars->a, krige_vars->r__, krige_vars->s, & ising);
	}
	/* 		Write a warning if the matrix is singular: */
	if (ising != 0) {
		if (general->idbg >= 1) {
			/* apagar */
			/*                  write(ldbg,*) 'WARNING SGSIM: singular matrix' */
			/*                  write(ldbg,*) '               for node',ix,iy,iz */
		}
		*cmean = *gmean;
		*cstdev = 1.f;
		return 0;
	}
	/* 		Compute the estimate and kriging variance.  Recall that kriging type */
	/* 			0 = Simple Kriging: */
	/* 			1 = Ordinary Kriging: */
	/* 			2 = Locally Varying Mean: */
	/* 			3 = External Drift: */
	/* 			4 = Collocated Cosimulation: */

	*cmean = 0.f;
	*cstdev = krige_vars->cbb;
	sumwts = 0.f;
	i1 = na;
	for (i = 1; i <= i1; ++i) {
		*cmean += (float) krige_vars->s[i - 1] * krige_vars->vra[i - 1];
		*cstdev -= (float) (krige_vars->s[i - 1] * krige_vars->rr[i - 1]);
		sumwts += (float) krige_vars->s[i - 1];
	}
	if (*lktype == 0) {
		*cmean += *gmean;
	}
	else if (*lktype == 1) {
		*cstdev -= (float) krige_vars->s[na];
	}
	else if (*lktype == 2) {
		*cmean += *gmean;
	}
	else if (*lktype >= 4) {
		ind = *ix + (*iy - 1) * general->nx + (*iz - 1) * general->nxy;
		*cmean += (float) krige_vars->s[na] * (bestAICube[ind] - simulation->vmedexp);
		*cstdev -= (float) (krige_vars->s[na] * krige_vars->rr[na]);
		*cmean += simulation->vmedexp;
	}
	/* Error message if negative variance: */
	if (*cstdev < 0.f) {
		/* apagar */
		/*            write(ldbg,*) 'ERROR: Negative Variance: ',cstdev */
		*cstdev = 0.f;
	}
	*cstdev = sqrt((double) MAX(*cstdev, 0.f)); /* sqrt((dmax(*cstdev,0.f))); */
	/*      Write out the kriging Weights if Seriously Debugging: */

	/*      if(idbg.ge.3) then */
	/*            do i=1,na */
	/*                  write(ldbg,140) i,vra(i),s(i) */
	/*            end do */
	/* 140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4) */
	/*            if(lktype.ge.4) write(ldbg,141) bestAICube(ind),s(na+1) */
	/* 141        fomat(' Sec Data  value ',f8.4,' weight ',f8.4) */
	/*            wrirte(ldbg,142) gmean,cmean,cstdev */
	/* 142        format(' Global mean ',f8.4,' conditional ',f8.4, */
	/*     +             ' std dev ',f8.4) */
	/*      end if */
	/* Finished Here: */
	return 0;
} /* krige_ */

