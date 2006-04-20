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
 * getindx()
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


int readdata(float *lvm,
             double *hard_data,
             unsigned int  hard_data_size,
	     general_vars_t * general, 
	     search_vars_t * search,
	     simulation_vars_t * simulation)
{
	/* Table of constant values */


	/* System generated locals */
	int i1;
	float r1;


	/* Local variables */
	int i, j;
	double w, cp, av;
	int ix, nt, iy, iz;
	double ss;
	float var[4], vrg, twt;
	int iend, ierr;
	unsigned int nelem;
	double oldcp;
	float vmedg;
	int index, inull;
	float vvarg;
	int icolvr, icolwt, istart;
	/* float testind; */
	int testind;    /* LPL */


	/* Parameter adjustments */
	--hard_data;
	--lvm;

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
		general->itrans = 0;
		search->sstrat = 1;
	}
	general->nd = 0;
	av = 0;
	ss = 0;
	/* !Establish the reference histogram for the simulation (provided that */
	/* !we have data, and we are transforming the data): */
	if (general->itrans == 1) {
		/* !write(*,*) 'Setting up transformation table' */
		/* !Decide which file to use for establishing the transformation table: */
		if (general->ismooth == 1) {
			/* !MARTELADO! Alterar depois, passar do C para aqui */
			/* !tmpfl='smoothFile.out' */
			icolvr = general->isvr;
			icolwt = general->iswt;
		} else {
			/* !tmpfl = DATA_FILENAME */
			icolvr = general->ivrl;
			icolwt = general->iwt;
		}
		/* !Now, read in the actual data: */
		nt = 0;
		general->ntr = 0;
		twt = 0;
		nelem = 0;
        printf_dbg2("readdata(): L3 %d\n", hard_data_size);
L3:
		if (nelem == hard_data_size) {
			goto L4;
		}
		i1 = general->nvari;
		for (j = 1; j <= i1; ++j) {
			var[j - 1] = hard_data[nelem + j];
		}
		nelem += general->nvari;
		/* !Trim this data? */
		if (var[icolvr - 1] < general->tmin || var[icolvr - 1] > general->tmax) {
			++nt;
			goto L3;
		}
		++general->ntr;
		/* !Removed after dynamic memory allocation */
		/* ! */
		/* !Exceeded available storage? */
		/* !if(ntr.gt.MAXDAT) then */
		/* !	write(*,*) ' ERROR: not enough room for data ',MAXDAT */
		/* !	stop */
		/* !endif */
		if (icolvr > general->nvari || icolwt > general->nvari) {
			fprintf(stderr,"ERROR: too few columns in ref data\n");
			fprintf(stderr,"\taborting.\n");
			return -1; /* ERROR */
		}
		/* !Keep this data: Assign the data value and coordinate location: */
		general->vrtr[general->ntr - 1] = var[icolvr - 1];
		if (icolwt <= 0) {
			general->vrgtr[general->ntr - 1] = 1;
		} else {
			general->vrgtr[general->ntr - 1] = var[icolwt - 1];
		}
		if (general->vrgtr[general->ntr - 1] <= 0) {
			--general->ntr;
			++nt;
			goto L3;
		}
		twt += general->vrgtr[general->ntr - 1];
		/* !Go back for another datum: */
		goto L3;
L4:
        printf_dbg2("readdata(): L4\n");
		if (general->ntr <= 1) {
			fprintf(stderr,"EROOR: too few data for transformation\n");
			fprintf(stderr,"\taborting.\n");
			return -1; /* ERROR */
		}
		/* !Sort data by value: */
		istart = 1;
		iend = general->ntr;
		/* sortem(&istart, &iend, general->vrtr, &one, general->vrgtr, &c__, &d__, &e, &f, &g, &h__); */
		sort_permute_float(istart, iend, general->vrtr, general->vrgtr);
		/* !Compute the cumulative probabilities and write transformation table */
		twt = (double) MAX(twt,1e-20f); /* dmax(twt,1e-20f); */
		oldcp = 0;
		cp = 0;
		vmedg = 0;
		vvarg = 0;
		for (j = 0; j < iend; ++j) {
			cp += (double) (general->vrgtr[j] / twt);
			w = (cp + oldcp) * .5f;
			gauinv(&w, &vrg, &ierr);
			if (ierr == 1) {
				vrg = general->nosvalue;
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
		/* !do j=istart,iend */
		/* !	write(lout,201) vrtr(j),(vrgtr(j)-vmedg)/sqrt(vvarg) */
		/* !201				format(2F12.5) */
		/* !end do */
		/* !close(lout) */
	}
	if (hard_data_size != 0) {
		/* !it means that there is no Hard Data file! */
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
        printf_dbg2("readdata(): L5\n");
L5:
		if (nelem == hard_data_size) {
			goto L6;
		}
		for (j = 1; j <= general->nvari; ++j) {
			var[j - 1] = hard_data[nelem + j];
		}
		nelem += general->nvari;
		if (var[general->ivrl - 1] < general->tmin || var[general->ivrl - 1] 
				> general->tmax) {
			++nt;
			goto L5;
		}
		++general->nd;
		if (general->nd > general->maxdat) {
			fprintf(stderr,"ERROR execeeded MAXDAT - check config file\n");
			return -1; /* ERROR */
		}
		/* !Acceptable data, assign the value, X, Y, Z coordinates, and weight: */
		general->vr[general->nd - 1] = var[general->ivrl - 1];
		if (general->ixl <= 0) {
			general->x[general->nd - 1] = general->xmn;
		} else {
			general->x[general->nd - 1] = var[general->ixl - 1];
		}
		if (general->iyl <= 0) {
			general->y[general->nd - 1] = general->ymn;
		} else {
			general->y[general->nd - 1] = var[general->iyl - 1];
		}
		if (general->izl <= 0) {
			general->z__[general->nd - 1] = general->zmn;
		} else {
			general->z__[general->nd - 1] = var[general->izl - 1];
		}
		if (general->iwt <= 0) {
			general->wt[general->nd - 1] = 1.f;
		} else {
			general->wt[general->nd - 1] = var[general->iwt - 1];
		}
		if (general->isecvr <= 0) {
			general->sec[general->nd - 1] = general->nosvalue;
		} else {
			general->sec[general->nd - 1] = var[general->isecvr - 1];
		}
		twt += general->wt[general->nd - 1];
		av += var[general->ivrl - 1] * general->wt[general->nd - 1];
		ss += var[general->ivrl - 1] * var[general->ivrl - 1] * general->wt[
			general->nd - 1];
		goto L5;
L6:
        printf_dbg2("readdata(): L6\n");
		if (general->imask == 1) {
			fprintf(stderr,"imaks = 1\n");
		}
			/* !do index=1, MASK_DATA_SIZE */
			/* !if(MASK_DATA(index).eq.0) then */
			/* !	write(*,*) 'MASK_DATA=', MASK_DATA(index) */
			/* !end if */
			/* !	mask(index) = MASK_DATA(index) */
			/* !end do */
			/* !index=0 */
			/* !do index=1, MASK_SIZE_DATA */
			/* !	write(*,*) 'mask=', mask(index) */
			/* !	write(*,*) ' MASK_DATA=', MASK_DATA(index) */
			/* !end do */
		simulation->vmedexp = 0;
		simulation->vvarexp = 0;
		for (i = 0; i < general->nd; ++i) {
			simulation->vmedexp += general->vr[i];
		}
		simulation->vmedexp /= general->nd;
		i1 = general->nd;
		for (i = 0; i < general->nd; ++i) {
			/* Computing 2nd power */
			r1 = general->vr[i] - simulation->vmedexp;
			simulation->vvarexp += r1 * r1;
		}
		simulation->vvarexp /= general->nd;
		/* !Compute the averages and variances as an error check for the user: */
		av /= (double) MAX(twt,1e-20f); /* dmax(twt,1e-20f); */
		ss = ss / ((double) MAX(twt,1e-20f)) - av * av;
		/* !write (*,*) ' Data for SGSIM: ' */
		/* !write (*,*) ' Number of acceptable data  = ',nd */
		/* !write (*,*) ' Number trimmed             = ',nt */
		/* !write (*,*) ' Weighted Average           = ',av */
		/* !write (*,*) ' Weighted Variance          = ',ss */
	}
	/* !Read secondary attribute model if necessary: */
	if (general->ktype >= 2) {
		/* !write(*,*) 'Reading secondary attribute file' */
		/*
		if (*bai_data_size == 0) {
			fprintf(stderr,"WARNING secondary attribute file does not exist!\n");
			return -1;
		}
		*/
		index = 0;
		av = 0;
		ss = 0;
		nelem = 0;
		/* lvm is allready set to a pointer of BAI on dssdll */
		/* !BAI_DATA_SIZE = nx*ny*nz */
		/*
		i1 = *bai_data_size__;
		for (index = 1; index <= i1; ++index) {
			i2 = general->nvaril;
			for (j = 1; j <= i2; ++j) {
				var[j - 1] = bai_data[nelem + j];
			}
			nelem += general->nvaril;
			sim[index] = (float) index;
			lvm[index] = var[general->icollvm - 1];
			av += var[general->icollvm - 1];
			ss += var[general->icollvm - 1] * var[general->icollvm - 1];
		}
		*/
		r1 = (float) general->nxyz;
		av /= ((double) MAX(r1,1)); /*dmax(r1,1.f); */
		r1 = (float) general->nxyz;
		ss = ss / ((double) MAX(r1,1)) - av * av;
		if (general->ktype == 5) {
			/* do nothing, dssdll allready sets clc acordingly to BCM.
			 *  older method was to make a linear iteration
			 *  over the BCM array and copy each position to
			 *  the clc array.
			 */
		}
		/* !Do we need to work with data residuals? (Locally Varying Mean) */
		else if (general->ktype == 2 ||
				 general->ktype == 3) {
			for (i = 0; i < general->nd; ++i) {
				getindx(&general->nx, &general->xmn, &general->xsiz, &general->x[i], &ix, &testind);
				getindx(&general->ny, &general->ymn, &general->ysiz, &general->y[i], &iy, &testind);
				getindx(&general->nz, &general->zmn, &general->zsiz, &general->z__[i], &iz, &testind);
				index = ix + (iy - 1) * general->nx + (iz - 1) * general->nxy;
				general->sec[i] = lvm[index];
				/* !Calculation of residual moved to krige subroutine: vr(i)=vr(i)-sec(i) */
			}
		}
	}
	/* !Re-scale secondary variable to mean and variance */
	/* !of the primary variable */
	if (general->ktype >= 4) {
		simulation->vmedsec = 0;
		simulation->vvarsec = 0;
		inull = 0;
		for (i = 1; i <= general->nxyz; ++i) {
			if (lvm[i] != -999.25f) {   /* > -999.25? */ 
				simulation->vmedsec += lvm[i];
			} else {
				++inull;
			}
		}
		simulation->vmedsec /= general->nxyz - inull;
		for (i = 1; i <= general->nxyz; ++i) {
			if (lvm[i] != -999.25f) {
				/* Computing 2nd power */
				r1 = lvm[i] - simulation->vmedsec;
				simulation->vvarsec += r1 * r1;
			}
		}
		simulation->vvarsec /= general->nxyz - inull;
		for (i = 1; i <= general->nxyz; ++i) {
			if (lvm[i] != -999.25f) {
				lvm[i] = (lvm[i] - simulation->vmedsec) / sqrt(   /* fuck-up */
						simulation->vvarsec) * sqrt(simulation->vvarexp) + 
					simulation->vmedexp;
			}
		}
	}

	return 0;
} /* readdata_ */


