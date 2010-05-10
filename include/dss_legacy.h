#ifndef _DSS_LEGACY_H
#define _DSS_LEGACY_H

#include "dss.h"
#include "dssutils.h"
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
           mask_t               *,
	   general_vars_t 	*,
	   harddata_t		*,
	   search_vars_t 	*,
	   covariance_vars_t *,
	   covtable_lookup_vars_t *,
	   krige_vars_t 	*);

extern int setrot(float, float, float , float , float,  int, double rotmat[5][3][3]);

extern int srchnod(int, int, int, float *,
		general_vars_t *,
		covtable_lookup_vars_t *,
		search_node_t *);

extern int gauinv(double, float *result);



/* dss_krige */
extern int krige(int , int , int , float , float , float ,
		int , float , float *, float *, float *, float ,
		general_vars_t *,
		harddata_t *,
		covariance_vars_t *,
		covtable_lookup_vars_t *,
		krige_vars_t *,
		search_node_t *);

extern int ksol(int neq, double *a, double *r, double *s);

extern double cova3(float x1, float y1, float z1, float x2, float y2, float z2,
					int varnum, float nugget, variogram_t *variogram,
					double rotmat[5][3][3], double *cmax);

extern int covtable(general_vars_t *, search_vars_t *, covariance_vars_t *,
                    covtable_lookup_vars_t *, krige_vars_t *);

extern double sqdist(float, float, float, float, float, float, int, double rotmat[5][3][3]);

extern int sort_permute_int(int , int , float *, int *);



/* dss_backtr */
extern float backtr(float vrgs, int nt, harddata_point_t *point, float min_value, float max_value);


#endif /* _DSS_LEGACY_H */
