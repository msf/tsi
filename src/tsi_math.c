#include <stdio.h>
#include <math.h>
#include "debug.h"
#include "tsi_math.h"

void update_best_grids(unsigned int size, float *BCM, float *BAI, float *new_ai, float *new_corr)
{
	int x, y, z;
	int i;

	for(i = 0; i < size; i++) {
		if(BCM[i] < new_corr[i])
		{
  		// new best correlation 
			BCM[i] = new_corr[i];
			BAI[i] = new_ai[i];
		}
	}
	
	/*
	for(x=0; x<BCM->x_max; x++)
		for(y=0; y<BCM->y_max; y++)
			for(z=0; z<BCM->z_max; z++) {
				if (BCM->values[BCM->Coords(x,y,z)] < new_corr->values[BCM->Coords(x,y,z)]) {
					// new best correlation 
					BCM->values[BCM->Coords(x,y,z)] = new_corr->values[BCM->Coords(x,y,z)];
					BAI->values[BCM->Coords(x,y,z)] = new_ai->values[BCM->Coords(x,y,z)];
				}
			}
	*/
}

/**
 * calculate the global correlation between two grids
 */ 
double grid_correlation(float *A, float *B, unsigned int size) {
    unsigned int i;
    double sum_AB = 0;
    double sum_A =  0;
    double sum_A2 = 0;
    double sum_A_2 = 0;
    double sum_B = 0;
    double sum_B2 = 0;
    double sum_B_2 = 0;
    double denom = 0;

    for (i = 0; i < size; i++) {
        sum_B  += B[i];
        sum_B2 += B[i] * B[i];
        sum_A  += A[i];
        sum_A2 += A[i] * A[i];
        sum_AB += A[i] * B[i];
    }
    sum_A_2 = sum_A * sum_A;
    sum_B_2 = sum_B * sum_B;

    //printf("sum_SEISMIC = %f\t sum_SY = %f\nsum_SEISMIC2 = %f\t sum_SY2 = %f\nsum_AB = %f\n", sum_A,sum_B,sum_A2,sum_B2,sum_AB);
    denom = sqrt(((size * sum_A2) - sum_A_2) *
                 ((size * sum_B2) - sum_B_2));

    if (denom > 0)
        return (((size * sum_AB) - (sum_A * sum_B)) / denom);
    else
        printf_dbg("calculateGlobalCorr: ERROR: division by zero, returning 0\n");
                                                                                                                                                                        
    return 0;
} /* grid_correlation */




/* Calc the nth root of x, x**(1/n) */
double nth_root(double xx, int nn) {
    int vi, i, doitagain;
    double vdelta,diff,testval,midpnt,lower,upper;
    double xu,epsilon,sign;

    epsilon = 0.00000001;
    xu = xx;
    sign = 1.0;
    if (xu < 0.0) {
        sign = -1.0;
        xu = -xu;
    }
    lower = 0.0;
    if (xu >= 1.0) {
        upper = 1.5*sqrt(xu);   /* use tagged sqrt() */
    } else {
        upper = 1.0;
    }
    doitagain = 1;

    i = 0;          /* binary search loop */
    while ((i < 1000) && (doitagain == 1)) {
        midpnt = lower + ((upper - lower)/2.0);
	testval = 1.0;
	for (vi = 0; vi < nn; vi ++) testval = midpnt * testval;
	diff = testval - xu;
	vdelta = diff;
	if (vdelta < 0.0) vdelta = 0.0 - vdelta;
        if (vdelta > epsilon) doitagain = 1;
        if (vdelta <= epsilon) doitagain = 0;
        if (diff < 0.0) lower = midpnt;
        if (diff > 0.0) upper = midpnt;
        if (diff == 0.0) doitagain = 0;
        i = i + 1;
    }
    return(sign*midpnt);
} /* nth_root */

/* end of file tsi_math.c */
