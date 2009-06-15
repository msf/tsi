#ifndef _DSSUTILS_H_
#define _DSSUTILS_H_

#include "dss.h"

int getPos(int nx, int ny, int x, int y, int z);
void get3Dcoords(int nx , int ny, int index, int *x, int *y, int *z);

int getIndex(float min, float siz, float loc);
float getAbsolutePos(float base, float siz, int index);

float compute_gaussian_equiv(float cmean, unsigned size, harddata_point_t *point);

int cmpfloat(const void *a1, const void *b1);

int cmpharddata_point_val(const void *a, const void *b);
int cmpharddata_point_gauss_cprob(const void *a, const void *b);

int cmpvalue_index(const void *a, const void *b);


#endif
