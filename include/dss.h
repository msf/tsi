#ifndef _DSS_H
#define _DSS_H

#include "debug.h"
#include "grid_heap.h"
#include "log.h"
#include "registry.h"


typedef struct general_vars_type {
    int    nd,      /* init to 0, inc to harddata */
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
           *z,      /* same size as harddata */
           *vr,     /* same size as harddata */
           *wt,     /* same size as harddata */
           *vrtr,   /* same size as harddata */
           *vrgtr,  /* same size as harddata */
           *sec;    /* same size as harddata */
	
    /* new wells hard data containers */
    unsigned int wellsNPoints, *wellsDataPos; 
    float *wellsDataVal;

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
    int    nclose;   /* LPL SHOULD BE INT!!!!! */

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
    float  vmedexp,    /* average experimental (wells) */
           vvarexp,    /* variance... */
           vmedsec,    /* secondary (AI data) */
           vvarsec;    /* secondary...*/

    /* parameters from registry */
    int    nsim;         /* const = 1, number of simulations */
    int    nsim_bk;
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

    /* backup values */
    int    nodmax_bk;  /* SEARCH:NODMAX */
    int    ntry_bk,    /* QUALITY:NTRY */
           icmean_bk,  /* QUALITY:ICMEAN */
           icvar_bk;   /* QUALITY:ICVAR */
} covtable_lookup_vars_t;


typedef struct krige_vars_type {
    double rotmat[5][3][3];

    double *rr,
           *r__,    /* used to set and solve a system of equations */
           *s,      /* search->nclose * covtable_lookp->ncnode */
           *a;

    float  *vra,
           *vrea,
           cbb;
           
    int    last_na;
} krige_vars_t;


typedef struct dss_type {
    /* auxiliar objects */
    registry *reg;            /* reference to the registry */
    grid_heap *heap;          /* reference to the grid heap */
    log_t *l;                   /* reference to the log */

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

    /* harddata */
    float *harddata;
    unsigned int    harddata_size;
} dss;

struct dss_check {
    void *general,
         *search,
         *simulation,
         *covariance,
         *clookup,
         *krige;
};

dss *new_dss(registry *r, grid_heap *h, log_t *l);

int run_dss(dss *d, float *AI);

int run_codss(dss *d, float *currBAI, float *currBCM, float *AI);

void delete_dss(dss *d);

int dss_parameters(dss *d, registry *r);

#endif /* _DSS_H */
