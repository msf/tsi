#ifndef _DSS_H
#define _DSS_H

#include "debug.h"
#include "grid_heap.h"
#include "log.h"
#include "registry.h"


typedef struct general_vars_type {
    int    nd,      // number of acceptance harddata 
           ntr;     // init to 0, inc to harddata 
    float  close[7000]; // same size as harddata??? 

    /* used by "readdata" for harddata processing */
    int    maxdat;  // harddata size... 
    float  *x,      // same size as harddata 
           *y,      // same size as harddata 
           *z,      // same size as harddata 
           *vr,     // hard_data, orig order
           *wt,     // weights
           *vrtr,   // hard_data, sorted by value 
           *vrgtr,  // weights of sorted, transformed to comulative probabilty. 
           *sec;    // NO SIM VALUE -> no secondary data.

    /* new wells hard data containers */
    unsigned int wellsNPoints, *wellsDataPos; 
    float *wellsDataVal;

    /* parameters from registry */
#define NVARI	4
#define IXL		1
#define IYL		2
#define IZL		3
#define IVRL	4
#define IWT		0
#define ISECVR	0
    float  min_value;      /* HARDDATA:TMIN */
    float  max_value;      /* HARDDATA:TMAX */

#define ITRANS	1
#define ISMOOTH 0
#define ISVR    1
#define ISWT    2
#define LTAIL	1
#define UTAIL	1

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

    float  nosim_value;  /* MASK:NULL_VALUE */

    int    ktype;     /* =1 (and 5)- KRIG:KTYPE */
} general_vars_t;


typedef struct search_vars_type {
    int    nclose;   /* LPL SHOULD BE INT!!!!! */

    /* parameters from registry */
    int    ndmin,    /* SEARCH:NDMIN */
           ndmax;    /* SEARCH:NDMAX */
    float  radius,   /* SEARCH:RADIUS */
           radsqd,   /* radius^2 */
           sang1,    /* SEARCH:SANG */
           sang2,    /* SEARCH:SANG1 */
           sang3,    /* SEARCH:SANG2 */
           sanis1,   /* (SEARCH:RADIUS1) / radius */
           sanis2;   /* (SEARCH:RADIUS2) / radius */
} search_vars_t;


typedef struct simulation_vars_type {
    float  vmedexp,    /* average experimental (wells) */
           vvarexp,    /* variance... */
           vmedsec,    /* secondary (AI data) */
           vvarsec;    /* secondary...*/

} simulation_vars_t;


typedef struct covariance_vars_type {
    int    isrot;
    float  cmax;

    /* parameters from registry */
    int    nst;   /* VARIOGRAM:NUMBER */
    float  c0;    /* VARIOGRAM:NUGGET */

    int    *it;      /* VARIOGRAMn:TYPE */
    float  *cc,      /* VARIOGRAMn:COV */
           *ang1,    /* VARIOGRAMn:ANG1 */
           *ang2,    /* VARIOGRAMn:ANG2 */
           *ang3,    /* VARIOGRAMn:ANG3 */
           *aa,      /* VARIOGRAMn:AA */
           *anis1,   /* VARIOGRAMn:AA1 / (aa>1e-20 ? aa : 1e-20) */ 
           *anis2;   /* VARIOGRAMn:AA2 / (aa>1e-20 ? aa : 1e-20) */

    struct variogram_type *variogram;
} covariance_vars_t;

typedef struct variogram_type {
    int     type; /* it */
    float   cov;  /* cc */
    float   ang1;
    float   ang2;
    float   ang3;
    float   aa;
    float   anis1;
    float   anis2;
} variogram_t;

typedef struct covtable_lookup_vars_type {
    int    nctx,
           ncty,
           nctz;
    int    nlooku,
           ncnode,
           icnode[64]; // max size = NODMAX
    float  cnodex[64], // max size = NODMAX
           cnodey[64], // max size = NODMAX
           cnodez[64], // max size = NODMAX
           cnodev[64]; // max size = NODMAX

    /* auxiliar grids allocated in run_dss */
    float  *covtab;
    int    *ixnode,
           *iynode,
           *iznode;
    
    /* parameters from registry */
    int    nodmax;  /* SEARCH:NODMAX */
    int    ntry;    /* QUALITY:NTRY */

} covtable_lookup_vars_t;


typedef struct krige_vars_type {
    double rotmat[5][3][3];

    double *rr,
           *r,    /* used to set and solve a system of equations */
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
    log_t *l;                 /* reference to the log */

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
