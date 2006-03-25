#ifndef _DSS_H
#define _DSS_H

#include "debug.h"
#include "grid_heap.h"
#include "registry.h"


typedef struct general_vars_type {
    int    nd,     /* init to 0, inc to harddata */
           ntr,     /* init to 0, inc to harddata */
           idbg,
           lin,
           lout,
           ldbg,
           llvm;
    float  close[7000]; /* same size as harddata??? */

    /* used by "readdata" for harddata processing */
    int    maxdat;  /* harddata size... */
    float  *x,      /* same size as harddata */
           *y,      /* same size as harddata */
           *z__,    /* same size as harddata */
           *vr,     /* same size as harddata */
           *wt,     /* same size as harddata */
           *vrtr,   /* same size as harddata */
           *vrgtr,  /* same size as harddata */
           *sec;    /* same size as harddata */
	/* new wells hard data containers */
	unsigned int wellsNPoints, *wellsDataPos; 
	double *wellsDataVal;

    /* acorni data */
    int    *ixv;

    /* parameters from registry */
    int    nvari,     /* HARDDATA:NVARI */
           ixl,       /* HARDDATA:IXL */
           iyl,       /* HARDDATA:IYL */
           izl,       /* HARDDATA:IZL */
           ivrl,      /* HARDDATA:IVRL */
           iwt,       /* HARDDATA:IWT */
           isecvr;    /* HARDDATA:ISECVR */
    float  tmin,      /* HARDDATA:TMIN */
           tmax;      /* HARDDATA:TMAX */

    int    itrans,    /* HDTRANS:ITRANS */
           ismooth,   /* HDTRANS:ISMOOTH */
           isvr,      /* HDTRANS:ISVR */
           iswt;      /* HDTRANS:ISWT */
    int    ltail;     /* HDTRANS:LTAIL */
    float  ltpar;     /* HDTRANS:LTPAR */
    int    utail;     /* HDTRANS:UTAIL */
    float  utpar,     /* HDTRANS:UTPAR */
           zmin,      /* HDTRANS:ZMIN */
           zmax;      /* HDTRANS:ZMAX */

    int    nx,        /* GRID:XNUMBER */
           ny,        /* GRID:YNUMBER */
           nz;        /* GRID:ZNUMBER */
    float  xsiz,      /* GRID:XSIZE */
           ysiz,      /* GRID:YSIZE */
           zsiz,      /* GRID:ZSIZE */
           xmn,       /* GRID:XMN */
           ymn,       /* GRID:YMN */
           zmn;       /* GRID:ZMN */
    int    nxy,       /* = nx * ny */
           nxyz;      /* = nx * ny * nz */

    int    icollvm,   /* SOFT:ICOLLVM */
           nvaril;    /* SOFT:NVARIL */

    float  nosvalue;  /* MASK:NULL_VALUE */
    int    imask;     /* MASK:USE_MASK*/

    int    ktype;     /* KRIG:KTYPE */
    float  colocorr,  /* KRIG:COLOCORR */
           varred;    /* KRIG:VARRED */
} general_vars_t;


typedef struct search_vars_type {
    int  nclose;     /* LPL SHOULD BE INT!!!!! */

    /* parameters from registry */
    int    ndmin,    /* SEARCH:NDMIN */
           ndmax;    /* SEARCH:NDMIN */
    int    sstrat,   /* SEARCH_SSTRAT */
           noct,     /* SEARCH:NOCT */
           mults;    /* SEARCH:MULTS */
    float  radius,   /* SEARCH:RADIUS */
           radsqd,   /* radius^2 */
           sang1,    /* SEARCH:SANG */
           sang2,    /* SEARCH:SANG1 */
           sang3,    /* SEARCH:SANG2 */
           sanis1,   /* (SEARCH:RADIUS1) / radius */
           sanis2,   /* (SEARCH:RADIUS2) / radius */
           nmult;    /* SEARCH:NMULTS */
} search_vars_t;


typedef struct simulation_vars_type {
    float  vmedexp,
           vvarexp,
           vmedsec,
           vvarsec;

    /* parameters from registry */
    int    nsim;         /* const = 1 */
} simulation_vars_t;


typedef struct covariance_vars_type {
    int    isrot;
    float  cmax;

    /* parameters from registry */
    int    nst[1];   /* VARIOGRAM:NUMBER */
    float  c0[1];    /* VARIOGRAM:NUGGET */

    int    *it;      /* VARIOGRAMn:TYPE */
    float  *cc,      /* VARIOGRAMn:COV */
           *ang1,    /* VARIOGRAMn:ANG1 */
           *ang2,    /* VARIOGRAMn:ANG2 */
           *ang3,    /* VARIOGRAMn:ANG3 */
           *aa,      /* VARIOGRAMn:AA */
           *anis1,   /* VARIOGRAMn:AA1 / (aa>1e-20 ? aa : 1e-20) */ 
           *anis2;   /* VARIOGRAMn:AA2 / (aa>1e-20 ? aa : 1e-20) */
} covariance_vars_t;


typedef struct covtable_lookup_vars_type {
    int    nctx,
           ncty,
           nctz;
    int    nlooku,
           ncnode,
           icnode[64];
    float  cnodex[64],
           cnodey[64],
           cnodez[64],
           cnodev[64];

    /* auxiliar grids allocated in run_dss */
    float  *covtab;
    int    *ixnode,
           *iynode,
           *iznode;
    
    /* parameters from registry */
    int    nodmax;  /* SEARCH:NODMAX */
    int    ntry,    /* QUALITY:NTRY */
           icmean,  /* QUALITY:ICMEAN */
           icvar;   /* QUALITY:ICVAR */
} covtable_lookup_vars_t;


typedef struct krige_vars_type {
    double rotmat[45], /* was [5][3][3] */
           r__[129],
           rr[129],
           s[129],
           a[16641];
    float  vra[129],
           cbb;
} krige_vars_t;


typedef struct file_vars_type {
    char   transfl[256],
           smthfl[256],
           tmpfl[256],
           datafl[256],
           outfl[256], 
	   dbgfl[256],
           lvmfl[256],
           corrfl[256],
           str[256];
} file_vars_t;


typedef struct dss_type {
    registry  *reg;
    grid_heap *heap;

    /* auxiliar grids */
    int covtab_idx,
        ixnode_idx,
        iynode_idx,
        iznode_idx,
        order_idx;

    /* DSS legacy types */
    general_vars_t          *general;
    search_vars_t           *search;
    simulation_vars_t       *simulation;
    covariance_vars_t       *covariance;
    covtable_lookup_vars_t  *clookup;
    krige_vars_t            *krige;
    file_vars_t             *files;

    /* harddata */
    double *harddata;
    int    harddata_size;
} dss;


dss *new_dss(registry *r, grid_heap *h);

int setup_dss(dss *d, float *AI);

int run_dss(dss *d, float *AI);

int run_codss(dss *d, float *currBAI, float *currBCM, float *AI);

void delete_dss(dss *d);


#endif /* _DSS_H */
