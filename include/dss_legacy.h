#ifndef _DSS_LEGACY_H
#define _DSS_LEGACY_H


/* Float functions */
#define ZERO_THRESHOLD 1e-6
#define EQUAL(f1,f2) (f1 > f2 ? 0 : (f2 - f1) < ZERO_THRESHOLD)
#define ZERO(f) (f < ZERO_THRESHOLD ? f > (-ZERO_THRESHOLD) : 0)


/* DSS legacy functions */

/* dss */
extern int readWellsData(general_vars_t *, float *, unsigned int);



/* dss_sim */
extern int dssim(float *, float *, float *, int *, int *, general_vars_t *,
                 search_vars_t *, simulation_vars_t *, covariance_vars_t *,
                 covtable_lookup_vars_t *, krige_vars_t *);

extern int setrot(float, float, float , float , float,  int, double rotmat[5][3][3]);

extern int srchnod(int, int, int, float *, general_vars_t *, search_vars_t *, 
                   covtable_lookup_vars_t *);

extern int gauinv(double *, float *, int *);

extern int readdata(float *, unsigned int, general_vars_t *, search_vars_t *,
                    simulation_vars_t *);

extern int getPos(int, int, int, int, int);
extern void get3Dcoords(int, int, int, int*, int*, int*);

extern int getIndex(float min, float siz, float loc);

/* this is in dss-utils.c */
extern float compute_gaussian_equiv(float cmean, unsigned size, float *vrtr, float *vrgtr);



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

extern double cova3(float x1, float y1, float z1, float x2, float y2, float z2, int nst,
                   float c0, int *it, float *cc, float *aa, double rotmat[5][3][3],
                   double *cmax);

extern int covtable(int *, float *, general_vars_t *, search_vars_t *, covariance_vars_t *,
                    covtable_lookup_vars_t *, krige_vars_t *);

extern double sqdist(float, float, float, float, float, float, int, double rotmat[5][3][3]);

extern int sort_permute_int(int , int , float *, int *);



/* dss_backtr */
extern float backtr(float vrgs, int nt, float *vr, float *vrg, float zmin, float zmax,
                    int ltail, float ltpar, int utail, float utpar);

/*extern int sortemi(int *, int *,
           float *,
           int *, int *,
           float *, float *, float *, float *, float *, float *);
*/


#endif /* _DSS_LEGACY_H */
