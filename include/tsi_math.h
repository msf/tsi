#ifndef _TSI_MATH_H
#define _TSI_MATH_H

void update_best_grids(unsigned int size, float *bcm, float *bai, float *new_ai, float *new_corr);

double grid_correlation(float *A, float *B, unsigned int size);

double nth_root(double xx, int nn);

#endif /* _TSI_MATH_H */
