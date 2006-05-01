#ifndef _DSS_LEGACY_H
#define _DSS_LEGACY_H

/* Float functions */

#define ZERO_THRESHOLD 1e-6
#define EQUAL(f1,f2) (f1 > f2 ? 0 : (f2 - f1) < ZERO_THRESHOLD)
#define ZERO(f) (f < ZERO_THRESHOLD ? f > (-ZERO_THRESHOLD) : 0)

/* Original DSS functions */
//#define GETINDX(n, mn, sz, loc, ix) ((*(ix) = (int)((*(loc) - *(mn))/(*(sz)) + 1.5)) < 1 ? !(*(ix) = 1) : (*(ix) > *(n) ? !(*(ix) = *(n)) : 1))


/* DSS legacy functions */

extern double acorni(int *);

extern int *new_acorni(int );


extern double gcum(float );

extern int locate(float *, int *, int *, int *, float *, int *);

extern double powint(float *, float *, float *, float *, float *, float *);

extern double backtr(float *, int *, float *, float *, float *, float *, int *, float *, int *, float *);

extern int cova3(float *,
          float *,
          float *,
          float *,
          float *,
		  float *,
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

extern int getindx(int *, float *, float *, float *, int *, int *);

extern int krige(int *, int *, int *, float * , float *, float *,
		int *, float *, float *, float *, float *, float *,
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

extern int picksup(int *, float *, int *, float *,
			int *, float *, int *, int *, double *, float *, 
			int *, int *, int *, int *);


extern int readdata(float *,
             double *,
             unsigned int ,
             general_vars_t *,
             search_vars_t *,
             simulation_vars_t *);	

extern int readWellsData(general_vars_t *, double *, unsigned int);

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


extern int setrot(float *, float *, float *, float *, float *
			, int *, int *, double *);

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

extern double sqdist(float *,
              float *,
              float *,
              float *,
              float *,
              float *,
              int *,
              int *,
              double *);


extern int srchnod(int *, int *, int *, float *,
		general_vars_t *, 
		search_vars_t *, 
		covtable_lookup_vars_t *);

extern int srchsupr(float *xloc, float *yloc, float *zloc,
             float *radsqd,
		     int *irot, int *maxrot, double *rotmat,
             int *nsbtosr, int *ixsbtosr, int *iysbtosr, int *izsbtosr,
             int *noct,
             float *x, float *y, float *z,
             float *tmp,
             int *nisb,
             int *nxsup, float *xmnsup, float *xsizsup,
             int *nysup, float *ymnsup, float *ysizsup,
             int *nzsup, float *zmnsup, float *zsizsup,
             int *nclose, float *close,
             int *infoct);


#endif /* _DSS_LEGACY_H */
