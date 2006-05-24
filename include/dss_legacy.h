#ifndef _DSS_LEGACY_H
#define _DSS_LEGACY_H

/* Float functions */

#define ZERO_THRESHOLD 1e-6
#define EQUAL(f1,f2) (f1 > f2 ? 0 : (f2 - f1) < ZERO_THRESHOLD)
#define ZERO(f) (f < ZERO_THRESHOLD ? f > (-ZERO_THRESHOLD) : 0)

/* DSS legacy functions */

extern float backtr(float vrgs, int nt, float *vr, float *vrg, float zmin, float zmax, int ltail, float ltpar, int utail, float utpar);

float cova3(float x1, float y1, float z1, float x2, float y2, float z2,
		int *nst, float *c0, int *it, float *cc,
		float *aa, double rotmat[5][3][3], float *cmax);
/*
extern int cova3(float *, float *, float *,
          float *, float *, float *,
           int *,
          int *,
          int *,
           float *,
           int *,
           float *,
           float *,
          int *,
          int *,
           double *,
          float *,
		  float *);
*/
 
extern int covtable(int *, float *,
	     general_vars_t *,
	     search_vars_t *,
	     covariance_vars_t *,
	     covtable_lookup_vars_t *,
	     krige_vars_t *);


extern int dsslib(float *params, float *models, double *hard_data,
	          int hard_data_size, float *bcm_data, float *bai_data,
		  float *output_data);

extern inline int getPos(int, int, int, int, int);

extern int gauinv(double *, float *, int *);

extern int getIndex(float min, float siz, float loc);

extern int krige(int , int , int , float , float , float ,
		int , float , float *, float *, float *, float ,
		general_vars_t *,
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


extern int ksol(int *, int *, int *, double *, double *, double *, int *);

int picksup(int nxsup, float xsizsup, int nysup, float ysizsup,
		int nzsup, float zsizsup, int irot,
		double rotmat[5][3][3], float radsqd, int *nsbtosr, int * ixsbtosr,
		int *iysbtosr, int *izsbtosr);

extern int readdata(float *,
             unsigned int ,
             general_vars_t *,
             search_vars_t *,
             simulation_vars_t *);	

extern int readWellsData(general_vars_t *, float *, unsigned int);

extern int readparam(float *,
              float *, 
		      general_vars_t *,
              search_vars_t *,
              simulation_vars_t *,
		      covariance_vars_t *,
              covtable_lookup_vars_t *);



extern int dssim(float *,
          float *,
          float *,
          int *,
          int *,
		  general_vars_t *,
		  search_vars_t *,
		  simulation_vars_t *,
		  covariance_vars_t *,
		  covtable_lookup_vars_t *,
		  krige_vars_t *);


extern int setrot(float , float , float , float , float, 
		int, double rotmat[5][3][3]);

extern int setsupr(int *, float *, float *, int *,
			float *, float *, int *, float *, float *, int *, float *, 
			float *, float *, float *, float *, int *, float *, float *, float *,
			int *, int *, int *, int *, int *, float *, 
			float *, int *, float *, float *, int *, float *, float *);

extern int sortem(int *, int *,
           float *,
           int *, float *,
           float *, float *, float *, float *, float *, float *);

extern int sortemi(int *, int *,
           float *,
           int *, int *,
           float *, float *, float *, float *, float *, float *);

extern int sort_permute_float(int , int , float *, float *);
extern int sort_permute_int(int , int , float *, int *);

extern double sqdist(float , float , float ,
              float , float , float ,
              int , double rotmat[5][3][3]); 

extern int srchnod(int , int , int , float *,
		general_vars_t *, 
		search_vars_t *, 
		covtable_lookup_vars_t *);

extern int srchsupr(float xloc, float yloc, float zloc,
             float radsqd,
		     int irot, double rotmat[5][3][3],
             int nsbtosr, int *ixsbtosr, int *iysbtosr, int *izsbtosr,
             int noct,
             float *x, float *y, float *z,
             float *tmp,
             int *nisb,
             int nxsup, float xmnsup, float xsizsup,
             int nysup, float ymnsup, float ysizsup,
             int nzsup, float zmnsup, float zsizsup,
             int *nclose, float *close);


#endif /* _DSS_LEGACY_H */
