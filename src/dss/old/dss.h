/*
 *  DSS.H
 *
 *
 */

#ifndef TSI_DSS
#define TSI_DSS

#define NVARI 0
#define IXL 1
#define IYL 2
#define IZL 3
#define IVRL 4
#define IWT 5
#define ISECVR 6
#define TMIN 7
#define TMAX 8
#define ITRANS 9 
#define ISMOOTH 10
#define ISVR 11
#define ISWT 12
#define ZMIN 13
#define ZMAX 14
#define LTAIL 15
#define LTPAR 16
#define UTAIL 17
#define UTPAR 18
#define NTRY 19
#define ICMEAN 20
#define ICVAR 21
#define NX 22
#define XMN 23
#define XSIZ 24
#define NY 25
#define YMN 26
#define YSIZ 27
#define NZ 28
#define ZMN 29
#define ZSIZ 30
#define NOSVALUE 31
#define IMASK 32
#define SEED 33
#define NDMIN 34
#define NDMAX 35
#define NODMAX 36
#define SSTRAT 37
#define MULTS 38
#define NMULTS 39
#define NOCT 40
#define RADIUS 41
#define RADIUS1 42
#define RADIUS2 43
#define SANG 44
#define SANG1 45
#define SANG2 46
#define KTYPE 47
#define COLOCORR 48
#define VARRED 49
#define NVARIL 50
#define ICOLLVM 51
#define VARNUM 52 /* covariance->nst[0] */
#define NUGGET 53 /* covariance->c0[0] */
#define DSS_TOTAL_PARAMS_NUM 54

/* Common Block Declarations */

struct general_vars {
    int nx, ny, nz;
    float xsiz, ysiz, zsiz, xmn, ymn, zmn;
	int nxy, nxyz;
    int nd, itrans, ntr, idbg;
    float x[7000], y[7000], z__[7000], vr[7000], wt[7000], vrtr[7000];
    float vrgtr[7000];
    float close[7000], sec[7000];
    int lin, lout, ldbg, llvm, icollvm, nvari, nvaril, ktype;
    float colocorr;
    int ltail;
    float ltpar;
    int utail;
    float utpar, zmin, zmax, varred;
    int ixl, iyl, izl, ivrl, iwt, isecvr;
    float tmin, tmax;
    int ismooth, isvr, iswt;
    float nosvalue;
    int imask;
	unsigned int wellsNPoints, *wellsDataPos; /* new wells hard data containers */
	double *wellsDataVal;
} generl_;

typedef struct general_vars general_vars_t;


struct search_vars {
    float radius,
          radsqd,
          sang1,
          sang2,
          sang3,
          sanis1,
          sanis2,
          nmult;
    int sstrat,
          noct,
          mults;
    int nclose;     /* LPL SHOULD BE INT!!!!! */
    int ndmin;
    int ndmax;
} search_;

typedef struct search_vars  search_vars_t;


struct simulation_vars {
    int nsim;
    float vmedexp, vvarexp, vmedsec, vvarsec;
} simula_;

typedef struct simulation_vars  simulation_vars_t;


struct covariance_vars {
    int nst[1], it[4];
    float cmax, c0[1], cc[4], aa[4], ang1[4], ang2[4], ang3[4], anis1[4], 
	    anis2[4];
    int isrot;
} cova3d_;

typedef struct covariance_vars covariance_vars_t;


struct covtable_lookup_vars {
    int nctx, ncty, nctz;
    float * covtab	/* was [300][300][560] */;
    int nlooku, ncnode, icnode[64];
    float cnodex[64], cnodey[64], cnodez[64], cnodev[64];
    int nodmax;
    int * ixnode, * iynode, * iznode;
    int ntry, icmean, icvar;
} clooku_;

typedef struct covtable_lookup_vars covtable_lookup_vars_t;


struct krige_vars {
    double rotmat[45]	/* was [5][3][3] */, r__[129], rr[129], s[129]
	    , a[16641];
    float vra[129], cbb;
} krigev_;

typedef struct krige_vars krige_vars_t;


struct file_vars {
    char transfl[256],
         smthfl[256],
         tmpfl[256],
         datafl[256],
         outfl[256], 
	 dbgfl[256],
         lvmfl[256],
         corrfl[256],
         str[256];
} files_;

typedef struct file_vars file_vars_t;


extern double acorni();

extern void newAcorni(int );


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


extern int coDss(float *, float *, double *, int *, float *, float *, float *);
extern int dss(float *, float *, double *, int *, float *);

extern int getPos(int, int, int, int, int);

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

extern int readparam(float *,
              float *, 
		      general_vars_t *,
              search_vars_t *,
              simulation_vars_t *,
		      covariance_vars_t *,
              covtable_lookup_vars_t *);

extern int readWellsData(general_vars_t *, double *, unsigned int);
		


extern int sdsim(float *,
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


#endif /* TSI_DSS */

