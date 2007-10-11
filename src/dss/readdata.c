#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"
#include "dss_legacy.h"
#include "debug.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)




/* ----------------------------------------------------------------------- */
/*     Read primary data and secondary data */
/* ----------------------------------------------------------------------- */

/** Funcoes utilizadas
 * sortem()
 * gauinv()
 * getIndex()
 */

/** CUBOS utilizados
 * lvm
 * hard_data
 */ 

/** CUBOS _nao_ utilizados
 * sim
 * clc
 * tmp
 * order
 * mask
 * bcm_data
 * bai_data
 */

/** structs globais utilizadas:
 * general
 * simulation
 * search
 */


int readdata(float *hard_data,
             unsigned int  hard_data_size,
	     general_vars_t * general, 
	     search_vars_t * search,
	     simulation_vars_t * simulation)
{
	/* Table of constant values */


	/* System generated locals */
	float r1;


	/* Local variables */
	int i, j;
	double av; /* weighted average */
	double ss; /* weighted variance */
	double w, cp;
	int nt;
	float vrg, twt;
	int iend, ierr;
	unsigned int nelem;
	double oldcp;
	float vmedg;
	float vvarg;
	int ivrl, iwt, istart;

	float *var;
	
	/* Parameter adjustments */
	--hard_data;

    printf_dbg2("readdata(): begin\n");
	/* Function Body */
	if (hard_data_size == 0) {
		/* !it means that there is no Hard Data file! */
		fprintf(stderr,"WARNING data file does not exist!");
		fprintf(stderr,"\t- creating an *unconditional simulation*\n");
		fprintf(stderr,"\t- Resetting ndmin, ndmax and itrans to 0\n");
		fprintf(stderr,"\t- Resetting sstrar to 1\n");
		search->ndmin = 0;
		search->ndmax = 0;
		//general->itrans = 0;
		//search->sstrat = 1;
	}
	
	var = (float *) tsi_malloc(sizeof(float)*general->nvari);
	general->nd = 0;
	av = 0;
	ss = 0;
	
	/* !Establish the reference histogram for the simulation (provided that */
	/* !we have data, and we are transforming the data): */
	if (general->itrans == 1) {
        /* ISMOOTH IS ALLWAYS 0 */
        /*
		if (general->ismooth == 1) {
			ivrl = general->isvr;
			iwt = general->iswt;
		} 
        else 
        */
        {
			ivrl = general->ivrl;
			iwt = general->iwt;
		}
		/* !Now, read in the actual data: */

		nt = 0;
		general->ntr = 0;
		twt = 0;
		nelem = 0;
        printf_dbg2("readdata(): read wells array %d\n", hard_data_size);
		while(nelem < hard_data_size) {

			for (j = 1; j <= general->nvari; ++j) {
				var[j-1] = hard_data[nelem + j];
			}
			nelem += general->nvari;
			/* !Trim this data? */
			if (var[ivrl - 1] < general->tmin || var[ivrl - 1] > general->tmax) {
				++nt;
				continue;
			}

			if (ivrl > general->nvari || iwt > general->nvari) {
				fprintf(stderr,"ERROR: too few columns in ref data\n");
				fprintf(stderr,"\taborting.\n");
				return -1; /* ERROR */
			}
			/* !Keep this data: Assign the data value and coordinate location: */
			general->vrtr[general->ntr] = var[ivrl - 1];
			if (iwt <= 0) {
				general->vrgtr[general->ntr] = 1;
			} else {
				general->vrgtr[general->ntr] = var[iwt - 1];
			}
			if (general->vrgtr[general->ntr] <= 0) {
				++nt;
				continue;
			}
			twt += general->vrgtr[general->ntr];
			++general->ntr;
			/* !Go back for another datum: */
		}

        printf_dbg2("readdata(): wells data transformation\n");
		if (general->ntr <= 1) {
			fprintf(stderr,"EROOR: too few data for transformation\n");
			fprintf(stderr,"\taborting.\n");
			return -1; /* ERROR */
		}
		/* !Sort data by value: */
		istart = 0;
		iend = general->ntr;
		sort_permute_float(istart, iend, general->vrtr, general->vrgtr);
		/* !Compute the cumulative probabilities and write transformation table */
		twt = MAX(twt,1e-20f); 
		oldcp = 0;
		cp = 0;
		vmedg = 0;// TODO: vmedg should be double, possible overflow
		vvarg = 0;
		for (j = 0; j < iend; ++j) {
			cp += (double) (general->vrgtr[j] / twt);
			w = (cp + oldcp) * .5f;
			ierr = gauinv(w, &vrg );
			if (ierr == 1) {
				vrg = general->nosim_value;
			}
			oldcp = cp;
			/* !Now, reset the weight to the normal scores value: */
			general->vrgtr[j] = vrg;
		}
		/* !Basic statistics calculation - mean and variance (DSSIM) */
		for (j = 0; j < iend; ++j) {
			vmedg += general->vrgtr[j]; 
		}
		vmedg /= iend - istart;
		for (j = 0; j < iend; ++j) {
			/* Computing 2nd power */
			r1 = general->vrgtr[j] - vmedg;
			vvarg += r1 * r1;
		}
		vvarg /= iend - istart;
	}

	if (hard_data_size != 0) {
		if ((general->ixl > general->nvari) ||
            (general->iyl > general->nvari) || 
		    (general->izl > general->nvari) ||
            (general->ivrl > general->nvari) ||
            (general->isecvr > general->nvari) ||
		    (general->iwt > general->nvari)) {
			fprintf(stderr,"ERROR: you have asked for a column number\n");
			fprintf(stderr,"\t\tgreater than available in file\n");
			return -1; /* ERROR */
		}
		/* !Read all the data until the end of the file: */
		twt = 0;
		general->nd = 0;
		nt = 0;
		nelem = 0;
		printf_dbg2("readdata(): ivrl: %d, ixl: %d, iyl: %d, izl: %d\n",
				general->ivrl,general->ixl,general->iyl,general->izl);
		while(nelem < hard_data_size) {

			for (j = 1; j <= general->nvari; ++j) {
				var[j - 1] = hard_data[nelem + j];
			}
			nelem += general->nvari;
			if (var[general->ivrl - 1] < general->tmin || 
			    var[general->ivrl - 1] > general->tmax) {
				++nt;
				continue;
			}
			if (general->nd >= general->maxdat) {
				fprintf(stderr,"ERROR execeeded MAXDAT - check config file\n");
				return -1; /* ERROR */
			}
			/* !Acceptable data, assign the value, X, Y, Z coordinates, and weight: */
			general->vr[general->nd] = var[general->ivrl - 1];

			if (general->ixl <= 0) {
				general->x[general->nd] = general->xmn;
			} else {
				general->x[general->nd] = var[general->ixl - 1];
			}
			if (general->iyl <= 0) {
				general->y[general->nd] = general->ymn;
			} else {
				general->y[general->nd] = var[general->iyl - 1];
			}
			if (general->izl <= 0) {
				general->z[general->nd] = general->zmn;
			} else {
				general->z[general->nd] = var[general->izl - 1];
			}
            /* NOT USED
			if (general->iwt <= 0) {
				general->wt[general->nd] = 1.f;
			} else {
				general->wt[general->nd] = var[general->iwt - 1];
			}
			if (general->isecvr <= 0) {
				general->sec[general->nd] = general->nosim_value;
			} else {
				general->sec[general->nd] = var[general->isecvr - 1];
			}
            */
			general->wt[general->nd] = 1.f;
			general->sec[general->nd] = general->nosim_value;

			twt += general->wt[general->nd];
			av += var[general->ivrl - 1] * general->wt[general->nd];
			ss += var[general->ivrl - 1] * var[general->ivrl - 1] * general->wt[general->nd];
			printf_dbg2("readdata: Wells Point: (%f, %f, %f) = %f\n",
					general->x[general->nd], 
                    general->y[general->nd], 
                    general->z[general->nd], 
                    general->vr[general->nd]);
			++general->nd;
		}
		
		tsi_free(var);

		simulation->vmedexp = 0;
		simulation->vvarexp = 0;
		for (i = 0; i < general->nd; ++i) {
			simulation->vmedexp += general->vr[i];
		}
		simulation->vmedexp /= general->nd;
		for (i = 0; i < general->nd; ++i) {
			/* Computing 2nd power */
			r1 = general->vr[i] - simulation->vmedexp;
			simulation->vvarexp += r1 * r1;
		}
		simulation->vvarexp /= general->nd;
		/* !Compute the averages and variances as an error check for the user: */
		av /= (double) MAX(twt,1e-20f); 
		ss = ss / ((double) MAX(twt,1e-20f)) - av * av;

		printf_dbg2("readdata(): vmedexp: %f, vvarexp: %f\n",
			   simulation->vmedexp, simulation->vvarexp);
		/* !write (*,*) ' Data for SGSIM: ' */
		/* !write (*,*) ' Number of acceptable data  = ',nd */
		/* !write (*,*) ' Number trimmed             = ',nt */
		/* !write (*,*) ' Weighted Average           = ',av */
		/* !write (*,*) ' Weighted Variance          = ',ss */
		printf_dbg2(" Number of acceptable data: %d\nNumber trimmed: %d\n Weighted Average: %f\n Weighted Variance: %f\n",
			general->nd, nt, av, ss);
	}
	

	return 0;
} /* readdata_ */


