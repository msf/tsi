#ifndef _DSS_H
#define _DSS_H

#include "debug.h"
#include "grid_heap.h"
#include "log.h"
#include "registry.h"


#define SIMPLE_KRIG		0
#define ORDINARY_KRIG	1
#define CO_KRIG			2

typedef struct harddata_point_type {
	float x,y,z;
	float val;				// value of point
	float gauss_cprob; 		// gaussian comulative probability
} harddata_point_t;

typedef struct {
	short x;
	short y;
	short z;
	float value;
} value_point_t;

typedef struct {
	unsigned int index;
	float   x, y, z;
	float	value;
} search_node_t;

typedef struct {
	unsigned int index;
	float	value;
} value_index_t;

typedef struct harddata_type {
	int		point_count;
	harddata_point_t *point;
	int		in_grid_count;
	value_index_t	*in_grid;
	float	min_value, max_value;
	double	average;
	double	variance;
} harddata_t;

typedef struct general_vars_type {
#define NVARI	4
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

typedef struct covariance_vars_type {
    int    isrot;
    float  cmax;

    /* parameters from registry */
    int    varnum;   /* VARIOGRAM:NUMBER */
    float  nugget;    /* VARIOGRAM:NUGGET */

    struct variogram_type *variogram;
} covariance_vars_t;

#define VARIOGRAM_TYPE_SPHERICAL	1
#define VARIOGRAM_TYPE_EXPONENCIAL	2
#define VARIOGRAM_TYPE_GAUSSIAN		3
#define VARIOGRAM_TYPE_POWER		4
#define VARIOGRAM_TYPE_HOLE			5
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
           ncnode;
	search_node_t	node[64]; // max size = NODMAX

    /* auxiliar grids allocated in run_dss */
    float  *covtab;
    short  *ixnode,
           *iynode,
           *iznode;
    
    /* parameters from registry */
    int    nodmax;  /* SEARCH:NODMAX */
    int    ntry;    /* QUALITY:NTRY */

} covtable_lookup_vars_t;


typedef struct krige_vars_type {
    double rotmat[5][3][3];
    float  cbb;
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
    covariance_vars_t       *covariance;
    covtable_lookup_vars_t  *clookup;
    krige_vars_t            *krige;
	harddata_t				*harddata;

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

int dss_parameters(dss *d, log_t *l, registry *r);

#endif /* _DSS_H */
