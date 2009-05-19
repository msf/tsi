#include <math.h>
#include "dss.h"
#include "dss_legacy.h"
#include "stdlib.h"
#include "memdebug.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define FALSE (0)



/* ----------------------------------------------------------------------- */
/*            Builds and Solves the SK or OK Kriging System */
/* INPUT VARIABLES: */
/*   ix,iy,iz        index of the point currently being simulated */
/*   xx,yy,zz        location of the point currently being simulated */
/* OUTPUT VARIABLES: */
/*   cmean           kriged estimate */
/*   std_deviation          kriged standard deviation */
/* EXTERNAL REFERENCES: ksol   Gaussian elimination system vra_weight */
/* ----------------------------------------------------------------------- */
int krige(int ix, int iy, int iz, float xx, float yy, float zz,
		int lktype, float global_mean, float *cmean, float * std_deviation, 
		float *bestAICube, float clcorr,
		general_vars_t 	* general,
		harddata_t		* harddata,
		covariance_vars_t * covariance,
        covtable_lookup_vars_t * covtable_lookup,
        krige_vars_t 	* krige_vars,
		search_node_t	*search_node)
{
	/* System generated locals */
	int i1;

	/* Local variables */
	int i, j;
	float x1, y1, z1, x2, y2, z2;
	int na, ii, jj, kk, in, ix1, iy1, iz1, ix2, iy2, iz2, ind, neq;
	float cov;
	int index;
	float sfmin, sfmax;
	float sumwts;

    /* matrices and temp data */
    double *cov_table  = (double *) tsi_malloc(16 + na * na * sizeof(double));
    double *cov_vector = (double *) tsi_malloc(8 + na * sizeof(double));
    double *vra_weight = (double *) tsi_malloc(8 + na * sizeof(double));
    double *rr =   (double *) tsi_malloc(8 + na * sizeof(double));
    double *vra =  (float *)  tsi_malloc(na * sizeof(float));



	/* Parameter adjustments */
	--bestAICube;

	/* init vars */
	ix1 = iy1 = iz1 = ix2 = iy2 = iz2 = 0;
	
	/* Function Body */
	na = covtable_lookup->ncnode;

	neq = na;
	if (lktype == ORDINARY_KRIG || lktype == CO_KRIG) {
		neq += 1;
	}

	/* Set up kriging matrices: */
	in = 0;
	for (j = 1; j <= na; ++j) {
		/* It is a previously simulated node (keep index for table look-up): */
		index = j - 1;
		ind = search_node[index].index;
		x1 = search_node[index].x;
		y1 = search_node[index].y;
		z1 = search_node[index].z;
		vra[index] = search_node[index].value;

		ix1 = covtable_lookup->ixnode[ind];
		iy1 = covtable_lookup->iynode[ind];
		iz1 = covtable_lookup->iznode[ind];

		ix1 += ix - covtable_lookup->nctx -1;
		iy1 += iy - covtable_lookup->ncty -1;
		iz1 += iz - covtable_lookup->nctz -1;

		index = getPos(ix1, iy1, iz1, general->nx, general->nxy);
		if(lktype == CO_KRIG) {
			vra[index] -= harddata->average;
		}
		if (lktype == SIMPLE_KRIG) {
			vra[index] -= global_mean;
		}
		for (i = 1; i <= j; ++i) {
			/* Sort out the actual location of point "i" */
			/* It is a previously simulated node (keep index for table
			   look-up */
			index = i - 1;
			ind = search_node[index].index;
			x2 = search_node[index].x;
			y2 = search_node[index].y;
			z2 = search_node[index].z;

			ix2 = covtable_lookup->ixnode[ind];
			iy2 = covtable_lookup->iynode[ind];
			iz2 = covtable_lookup->iznode[ind];

			ix2 += ix - covtable_lookup->nctx -1;
			iy2 += iy - covtable_lookup->ncty -1;
			iz2 += iz - covtable_lookup->nctz -1;

		
			/* Now, get the covariance value: */
			++in;

			/* Try to use the covariance look-up (if the distance is in range): */
			ii = covtable_lookup->nctx + 1 + (ix1 - ix2);
			jj = covtable_lookup->ncty + 1 + (iy1 - iy2);
			kk = covtable_lookup->nctz + 1 + (iz1 - iz2);
			if (ii < 1 || ii > general->nx ||
				jj < 1 || jj > general->ny ||
				kk < 1 || kk > general->nz) {
				double cmax, c;
				c = cova3(x1, y1, z1, x2, y2, z2, covariance->varnum,
						covariance->nugget, covariance->variogram,
						krige_vars->rotmat, &cmax);
				cov = (float) c;
				covariance->cmax = (float) cmax;
			} else {
				cov = covtable_lookup->covtab[getPos(ii,jj,kk, general->nx, general->nxy)];
			}
			cov_table[in - 1] = (double) cov;
		}
		/* Get the RHS value (possibly with covariance look-up table): */
		/* Try to use the covariance look-up (if the distance is in range): */
		ii = covtable_lookup->nctx + 1 + (ix - ix1);
		jj = covtable_lookup->ncty + 1 + (iy - iy1);
		kk = covtable_lookup->nctz + 1 + (iz - iz1);
		if (ii < 1 || ii > general->nx ||
			jj < 1 || jj > general->ny ||
			kk < 1 || kk > general->nz) {
			double cmax;
			cov_vector[j - 1] = cova3(xx, yy, zz, x1, y1, z1, covariance->varnum,
					covariance->nugget, covariance->variogram, 
					krige_vars->rotmat, &cmax);
			covariance->cmax = (float) cmax;
		} else {
			cov_vector[j - 1] = (double) covtable_lookup->covtab[getPos(ii,jj,kk, general->nx, general->nxy)];
		}
		
		rr[j - 1] = cov_vector[j - 1];
	}
	/* Addition of OK constraint: */
	if (lktype == ORDINARY_KRIG ) {
		i1 = na;
		for (i = 1; i <= i1; ++i) {
			++in;
			cov_table[in - 1] = 1.f;
		}
		++in;
		cov_table[in - 1] = 0.f;
		cov_vector[na] = 1.f;
		rr[na] = 1.f;
	}
	/* Addition of Collocated Cosimulation Constraint: */
	else if (lktype == CO_KRIG ) {
		sfmin = 1e21f;
		sfmax = -1e21f;
		for (i = 1; i <= na; ++i) {
			++in;
			cov_table[in - 1] = (double) clcorr * cov_vector[i - 1];
			if (cov_table[in - 1] < sfmin) {
				sfmin = (float) cov_table[in - 1];
			}
			else if (cov_table[in - 1] > sfmax) {
				sfmax = (float) cov_table[in - 1];
			}
		}
		++in;
		cov_table[in - 1] = 1.f;
		ii = na + 1;
		cov_vector[ii - 1] = (double) clcorr;
		rr[ii - 1] = cov_vector[ii - 1];
	}
	/* 		Write out the kriging Matrix if Seriously Debugging: */
	/* 		Solve the Kriging System: */
	if (neq == 1) {
		vra_weight[0] = cov_vector[0] / cov_table[0];
    } else {
		int ising;
		double *rp, *ap, *sp;   //LPL test...
        ap = cov_table;
        rp = cov_vector;
        sp = vra_weight;
        ap--;
        rp--;
        sp--;
        ising = ksol(neq, ap, rp, sp);
        /* 		Write a warning if the matrix is singular: */
        if (ising != 0) {
            printf("krige() ERROR: singular matrix for node (%d,%d,%d)\n",
                    ix, iy, iz);

            *cmean = global_mean;
            *std_deviation = 1.f;
            tsi_free(rr);
            tsi_free(cov_vector);
            tsi_free(vra_weight);
            tsi_free(cov_table);
            tsi_free(vra);
            return 0;
        }
    }
	/* 		Compute the estimate and kriging variance.  Recall that kriging type */
	/* 			0 = Simple Kriging: */
	/* 			1 = Ordinary Kriging: */
	/* 			2 = Locally Varying Mean: */
	/* 			3 = External Drift: */
	/* 			4 = Collocated Cosimulation: */

	*cmean = 0.f;
	*std_deviation = krige_vars->cbb;
	sumwts = 0.f;
	for (i = 0; i < na; ++i) {
		*cmean += (float) vra_weight[i] * vra[i];
		*std_deviation -= (float) (vra_weight[i] * rr[i]);
		sumwts += (float) vra_weight[i];
	}
	if (lktype == SIMPLE_KRIG) {
		*cmean += global_mean;
	}
	else if (lktype == ORDINARY_KRIG) {
		*std_deviation -= (float) vra_weight[na];
	}
	else if (lktype == CO_KRIG) {
		ind = getPos(ix, iy, iz, general->nx, general->nxy);
		*cmean += (float) vra_weight[na] * (bestAICube[ind] - harddata->average);
		*std_deviation -= (float) (vra_weight[na] * rr[na]);
		*cmean += harddata->average;
	}
	/* Error message if negative variance: */
	if (*std_deviation < -ZERO_THRESHOLD) {
		/* apagar */
		/*            write(ldbg,*) 'ERROR: Negative Variance: ',std_deviation */
		fprintf(stderr,"krige(): ERROR: Negative Variance: %f\n", *std_deviation);
	}
	*std_deviation = sqrt((double) MAX(*std_deviation, 0.f));
	/*      Write out the kriging Weights if Seriously Debugging: */

    /* NOTE: this can give information regarding purpose of cryptic vars like:
     *  vra[]
     *  s[]
     */
	/*      if(idbg.ge.3) then */
	/*            do i=1,na */
	/*                  write(ldbg,140) i,vra(i),s(i) */
	/*            end do */
	/* 140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4) */
	/*            if(lktype.ge.4) write(ldbg,141) bestAICube(ind),s(na+1) */
	/* 141        fomat(' Sec Data  value ',f8.4,' weight ',f8.4) */
	/*            wrirte(ldbg,142) global_mean,cmean,std_deviation */
	/* 142        format(' Global mean ',f8.4,' conditional ',f8.4, */
	/*     +             ' std dev ',f8.4) */
	/*      end if */
	/* Finished Here: */
    tsi_free(cov_table);
    tsi_free(cov_vector);
    tsi_free(vra_weight);
    tsi_free(rr);
    tsi_free(vra);
	return 0;
} /* krige_ */
