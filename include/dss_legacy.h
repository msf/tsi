#ifndef _DSS_LEGACY_H
#define _DSS_LEGACY_H

#include "dss.h"
#include "log.h"

/* Float functions */
#define ZERO_THRESHOLD 1e-6
#define EQUAL(f1,f2) (f1 > f2 ? 0 : (f2 - f1) < ZERO_THRESHOLD)
#define ZERO(f) (f < ZERO_THRESHOLD ? f > (-ZERO_THRESHOLD) : 0)


/* readdata */
extern int load_harddata_file(log_t *l, char *filename, harddata_t *);

extern int readdata(log_t *l, harddata_t *, general_vars_t *);
/* dss */


/* dss_sim */
extern int dssim(float *, float *, float *, int *, int ktype,
	   general_vars_t 	*,
	   harddata_t		*, 
	   search_vars_t 	*,
	   simulation_vars_t *,
	   covariance_vars_t *, 
	   covtable_lookup_vars_t *, 
	   krige_vars_t 	*);

extern int setrot(float, float, float , float , float,  int, double rotmat[5][3][3]);

extern int srchnod(int, int, int, float *, general_vars_t *, search_vars_t *, 
                   covtable_lookup_vars_t *);

extern int gauinv(double, float *result);


/* this is in dss-utils.c */
extern int getPos(int, int, int, int, int);

extern void get3Dcoords(int, int, int, int*, int*, int*);

extern int getIndex(float min, float siz, float loc);

extern float getAbsolutePos(float base, float siz, int index);

extern float compute_gaussian_equiv(float cmean, unsigned size, harddata_point_t *point);

extern int cmpfloat(const void *a1, const void *b1);

extern int cmpharddata_point_val(const void *a, const void *b);

extern int cmpharddata_point_gauss_cprob(const void *a, const void *b);

int cmpvalue_index(const void *a, const void *b);
/* dss_supr */
extern int setsupr(int *, float *, float *, int *, float *, float *, int *, float *,
                   float *, int *, float *, float *, float *, float *, float *, int *,
                   float *, float *, float *, int *, int *, int *, int *, int *, float *, 
                   float *, int *, float *, float *, int *, float *, float *);

extern int srchsupr(float xloc, float yloc, float zloc, float radsqd, int irot,
                    double rotmat[5][3][3],int nsbtosr, int *ixsbtosr, int *iysbtosr,
                    int *izsbtosr, int noct, float *x, float *y, float *z, float *tmp,
                    int *nisb, int nxsup, float xmnsup, float xsizsup, int nysup,
                    float ymnsup, float ysizsup, int nzsup, float zmnsup, float zsizsup,
                    int *nclose, float *close);

extern int picksup(int nxsup, float xsizsup, int nysup, float ysizsup, int nzsup,
                   float zsizsup, int irot, double rotmat[5][3][3], float radsqd,
                   int *nsbtosr, int * ixsbtosr, int *iysbtosr, int *izsbtosr);

extern int sortem(int *, int *, float *, int *, float *, float *, float *, float *,
                  float *, float *, float *);

extern int sort_permute_float(int , int , float *, float *);



/* dss_krige */
extern int krige(int , int , int , float , float , float ,
		int , float , float *, float *, float *, float ,
		general_vars_t *,
		harddata_t	*,
		search_vars_t *,
		simulation_vars_t *,
		covariance_vars_t *,
		covtable_lookup_vars_t *,
		krige_vars_t *);

extern int krige1(int *, int *, int *, float * , float *, float *,
		float *, float *, float *, float *, float *,
		general_vars_t *,
		search_vars_t *,
		simulation_vars_t *,
		covariance_vars_t *,
		covtable_lookup_vars_t *,
		krige_vars_t *);

extern int krige5(int *, int *, int *, float * , float *, float *,
		float *, float *, float *, float *, float *,
		general_vars_t *,
		search_vars_t *,
		simulation_vars_t *,
		covariance_vars_t *,
		covtable_lookup_vars_t *,
		krige_vars_t *);

extern int ksol(int nright, int neq, int nsb, double *a, double *r, double *s);
extern int ksol_opt(int neq, double *a, double *r, double *s);

extern double cova3(float x1, float y1, float z1, float x2, float y2, float z2,
					int varnum, float nugget, variogram_t *variogram,
					double rotmat[5][3][3], double *cmax);

extern int covtable(int *, float *, general_vars_t *, search_vars_t *, covariance_vars_t *,
                    covtable_lookup_vars_t *, krige_vars_t *);

extern double sqdist(float, float, float, float, float, float, int, double rotmat[5][3][3]);

extern int sort_permute_int(int , int , float *, int *);



/* dss_backtr */
extern float backtr(float vrgs, int nt, harddata_point_t *point, float min_value, float max_value);

/*extern int sortemi(int *, int *,
           float *,
           int *, int *,
           float *, float *, float *, float *, float *, float *);
*/


#endif /* _DSS_LEGACY_H */
