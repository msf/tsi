#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"
#include "dss_legacy.h"
#include "debug.h"
#include "log.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


#define IXL		1
#define IYL		2
#define IZL		3
#define IVRL	4
#define IWT		0
#define ISECVR	0
#define ITRANS	1
#define ISMOOTH 0
#define ISVR    1
#define ISWT    2


/* ----------------------------------------------------------------------- */
/*     Read primary data and secondary data */
/* ----------------------------------------------------------------------- */

/** Funcoes utilizadas
 * sortem()
 * gauinv()
 * getIndex()
 */


int readdata(log_t *l, float *hard_data,
             unsigned int  hard_data_size,
	     general_vars_t * general, 
	     search_vars_t * search,
	     simulation_vars_t * simulation)
{
	/* System generated locals */
	float r1;


	/* Local variables */
	int i, j;
	double avg; /* weighted average */
	double variance; /* weighted variance */
	double w;
	int nt;
	float vrg, twt;
	int iend, ierr;
	unsigned int nelem;
	int istart;

	float var[4];

    printf_dbg2("readdata(): begin\n");
	/* Function Body */
	if (hard_data_size == 0) {
		/* !it means that there is no Hard Data file! */
		fprintf(stderr,"WARNING data file does not exist!");
		fprintf(stderr,"\t- creating an *unconditional simulation*\n");
		fprintf(stderr,"\t- Resetting ndmin, ndmax to 0\n");
		fprintf(stderr,"\t- Resetting sstrar to 1\n");
		search->ndmin = 0;
		search->ndmax = 0;
	}
	
	/* !Establish the reference histogram for the simulation (provided that */
	/* !we have data, and we are transforming the data): */
	/* !Now, read in the actual data: */

    printf_dbg2("readdata(): read wells array %d\n", hard_data_size);
	i = 0;
	nt = 0;
	nelem = 0;
	simulation->vmedexp = 0;
	while(nelem < hard_data_size) {

		for (j = 0; j < NVARI; ++j) {
			var[j] = hard_data[nelem + j];
		}
		nelem += NVARI;
		/* !Trim this data? */
		if (var[3] < general->min_value || var[3] > general->max_value) {
			++nt;
			continue;
		}

		/* !Keep this data: Assign the data value and coordinate location: */
		/* !Acceptable data, assign the value, X, Y, Z coordinates, and weight: */
		general->x[i] = var[0];
		general->y[i] = var[1];
		general->z[i] = var[2];
		general->vr[i] = var[3];
		general->vrtr[i] = var[3];


		simulation->vmedexp += var[3];
		/*
		printf_dbg2("readdata: Wells Point: (%f, %f, %f) = %f\n",
				general->x[i], 
				general->y[i], 
				general->z[i], 
				general->vr[i]);
				*/
		i++;
	}
	general->maxdat = i;
	simulation->vmedexp /= general->maxdat;

    printf_dbg2("readdata(): wells data transformation\n");
	if (general->maxdat <= 1) {
		fprintf(stderr,"EROOR: too few data for transformation\n");
		fprintf(stderr,"\taborting.\n");
		return -1; /* ERROR */
	}
	/* !Sort data by value: */
	istart = 0;
	iend = general->maxdat;
	sort_permute_float( istart, iend, general->vrtr, general->vrgtr);
//	qsort(general->vrtr, general->maxdat, sizeof(float), cmpfloat);
	/* !Compute the cumulative probabilities and write transformation table */
	double cp = 0;
	double oldcp = 0;
	for (j = 0; j < general->maxdat; ++j) {
		cp = (double) (j+1) / (double) general->maxdat;
		w = (cp + oldcp) / 2;
		ierr = gauinv( w, &vrg);
		if (ierr == 1) {
			vrg = general->nosim_value;
		}
		oldcp = cp;
		/* !Now, reset the weight to the normal scores value: */
		general->vrgtr[j] = vrg;
		printf_dbg2("readdata(): vrtr[%u] = %f,\tvrgtr[%u] = %f\n",
				j, general->vrtr[j], j, general->vrgtr[j]);
	}

	/* !Read all the data until the end of the file: */

	simulation->vvarexp = 0;
	for (i = 0; i < general->maxdat; ++i) {
		/* Computing 2nd power */
		r1 = general->vr[i] - simulation->vmedexp;
		simulation->vvarexp += r1 * r1;
	}
	simulation->vvarexp /= general->maxdat;

	printf_dbg2("readdata(): vmedexp: %f, vvarexp: %f\n",
		   simulation->vmedexp, simulation->vvarexp);
	printf_dbg2(" Number of acceptable data: %d\nNumber trimmed: %d\n", general->maxdat, nt);

	return 0;
} /* readdata_ */


