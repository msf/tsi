#include <stdlib.h>
#include "debug.h"
#include "memdebug.h"
#include "grid_heap.h"
#include "registry.h"
#include "dss.h"
#include "dss_legacy.h"
#include "tsi_io.h"

int  dss_parameters(dss *d, registry *r )
{
    int i, varnum;
    char variogram[16];
    float aa1, aa2, aa, sill;
    reg_key *k;

    /* harddata parameters */
    if ((k = get_key(r, "HARDDATA", "NVARI")) == NULL) return 1;
    d->general->nvari = get_int(k);
    
    if ((k = get_key(r, "HARDDATA", "IXL")) == NULL) return 1;
    d->general->ixl = get_int(k);

    if ((k = get_key(r, "HARDDATA", "IYL")) == NULL) return 1;
    d->general->iyl = get_int(k);
    
    if ((k = get_key(r, "HARDDATA", "IZL")) == NULL) return 1;
    d->general->izl = get_int(k);
    
    if ((k = get_key(r, "HARDDATA", "IVRL")) == NULL) return 1;

    d->general->ivrl = get_int(k);
    if ((k = get_key(r, "HARDDATA", "IWT")) == NULL) return 1;
    d->general->iwt = get_int(k);

    if ((k = get_key(r, "HARDDATA", "ISECVR")) == NULL) return 1;
    d->general->isecvr = get_int(k);

    if ((k = get_key(r, "HARDDATA", "TMIN")) == NULL) return 1;
    d->general->tmin = get_float(k);

    if ((k = get_key(r, "HARDDATA", "TMAX")) == NULL) return 1;
    d->general->tmax = get_float(k);


    /* harddata transformations */
    if ((k = get_key(r, "HDTRANS", "ITRANS")) == NULL) return 1;
    d->general->itrans = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ISMOOTH")) == NULL) return 1;
    d->general->ismooth = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ISVR")) == NULL) return 1;
    d->general->isvr = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ISWT")) == NULL) return 1;
    d->general->iswt = get_int(k);

    if ((k = get_key(r, "HDTRANS", "ZMIN")) == NULL) return 1;
    d->general->zmin = get_float(k);

    if ((k = get_key(r, "HDTRANS", "ZMAX")) == NULL) return 1;
    d->general->zmax = get_float(k);

    if ((k = get_key(r, "HDTRANS", "LTAIL")) == NULL) return 1;
    d->general->ltail = get_int(k);

    if ((k = get_key(r, "HDTRANS", "LTPAR")) == NULL) return 1;
    d->general->ltpar = get_float(k);

    if ((k = get_key(r, "HDTRANS", "UTAIL")) == NULL) return 1;
    d->general->utail = get_int(k);

    if ((k = get_key(r, "HDTRANS", "UTPAR")) == NULL) return 1;
    d->general->utpar = get_float(k);

    if ((d->general->ltail != 1) && (d->general->ltail != 2)) {
        printf_dbg("ERROR invalid lower tail option: %d\n", (int) d->general->ltail);
        printf_dbg("\t\tonly allowed 1 or 2 - see manual.\n");
        return 1;
    }
    if ((d->general->utail != 1) && (d->general->utail != 2) && (d->general->utail != 4)) {
        printf_dbg("ERROR invalid upper tail option: %d\n", (int) d->general->utail);
        printf_dbg("\t\tonly allowed 1, 2 or 4 - see manual.\n");
        return 1;
    }
    if ((d->general->utail == 4) && (d->general->utpar < 1)) {
        printf_dbg("ERROR invalid power for hyperbolic tail: %f\n", d->general->utpar);
        printf_dbg("\t\tmust be greater than 1.0.\n");
        return 1;
    }
    if ((d->general->ltail == 2) && (d->general->ltpar < 0)) {
        printf_dbg("ERROR invalid power for power model: %f\n",d->general->ltpar);
        printf_dbg("\t\tmust be greater than 0.0.\n");
        return 1;
    }
    if ((d->general->utail == 2) && (d->general->utpar < 0)) {
        printf_dbg("ERROR invalid power for power model: %f\n",d->general->utpar);
        printf_dbg("\t\tmust be greater than 0.0.\n");
        return 1;
    }

    
    /* simulations quality */
    if ((k = get_key(r, "QUALITY", "NTRY")) == NULL) return 1;
    d->clookup->ntry = get_int(k);

    if ((k = get_key(r, "QUALITY", "ICVAR")) == NULL) return 1;
    d->clookup->icvar = get_int(k);

    if ((k = get_key(r, "QUALITY", "ICMEAN")) == NULL) return 1;
    d->clookup->icmean = get_int(k);

    
    /* grid parameters */
    if ((k = get_key(r, "GRID", "XNUMBER")) == NULL) return 1;
    d->general->nx = get_int(k);

    if ((k = get_key(r, "GRID", "YNUMBER")) == NULL) return 1;
    d->general->ny = get_int(k);

    if ((k = get_key(r, "GRID", "ZNUMBER")) == NULL) return 1;
    d->general->nz = get_int(k);

    if ((k = get_key(r, "GRID", "XCOORD")) == NULL) return 1;
    d->general->xmn = get_float(k);

    if ((k = get_key(r, "GRID", "YCOORD")) == NULL) return 1;
    d->general->ymn = get_float(k);

    if ((k = get_key(r, "GRID", "ZCOORD")) == NULL) return 1;
    d->general->zmn = get_float(k);

    if ((k = get_key(r, "GRID", "XSIZE")) == NULL) return 1;
    d->general->xsiz = get_float(k);

    if ((k = get_key(r, "GRID", "YSIZE")) == NULL) return 1;
    d->general->ysiz = get_float(k);

    if ((k = get_key(r, "GRID", "ZSIZE")) == NULL) return 1;
    d->general->zsiz = get_float(k);

    d->general->nxy = d->general->nx * d->general->ny;
    d->general->nxyz = d->general->nxy * d->general->nz;

    
    /* mask */
    if ((k = get_key(r, "MASK", "NULL_VALUE")) == NULL) return 1;
    d->general->nosvalue = get_float(k);

    if ((k = get_key(r, "MASK", "USE_MASK")) == NULL) return 1;
    d->general->imask = get_int(k);


    /* search parameters */
    if ((k = get_key(r, "SEARCH", "NDMIN")) == NULL) return 1;
    d->search->ndmin = get_int(k);

    if ((k = get_key(r, "SEARCH", "NDMAX")) == NULL) return 1;
    d->search->ndmax = get_int(k);

    if ((k = get_key(r, "SEARCH", "NODMAX")) == NULL) return 1;
    d->clookup->nodmax = get_int(k);

    if ((k = get_key(r, "SEARCH", "SSTRAT")) == NULL) return 1;
    d->search->sstrat = get_int(k);
    if (d->search->sstrat == 1) 
		d->search->ndmax = 0;

    if ((k = get_key(r, "SEARCH", "MULTS")) == NULL) return 1;
    d->search->mults = get_int(k);

    if ((k = get_key(r, "SEARCH", "NMULTS")) == NULL) return 1;
    d->search->nmult = get_float(k);

    if ((k = get_key(r, "SEARCH", "NOCT")) == NULL) return 1;
    d->search->noct = get_int(k);

    if ((k = get_key(r, "SEARCH", "RADIUS")) == NULL) return 1;
    d->search->radius = get_float(k);
    if (d->search->radius < 1e-20) {
        printf("load_dss_configs(): ERROR - radius (%f) must be greater than EPSLON (1E-20)\n", d->search->radius);
        return 1;
    }
    d->search->radsqd = d->search->radius * d->search->radius;

    if ((k = get_key(r, "SEARCH", "RADIUS1")) == NULL) return 1;
    d->search->sanis1 = get_float(k) / d->search->radius;

    if ((k = get_key(r, "SEARCH", "RADIUS2")) == NULL) return 1;
    d->search->sanis2 = get_float(k) / d->search->radius;

    if ((k = get_key(r, "SEARCH", "SANG1")) == NULL) return 1;
    d->search->sang1 = get_float(k);

    if ((k = get_key(r, "SEARCH", "SANG2")) == NULL) return 1;
    d->search->sang2 = get_float(k);

    if ((k = get_key(r, "SEARCH", "SANG3")) == NULL) return 1;
    d->search->sang3 = get_float(k);


    /* kriging parameters */
    if ((k = get_key(r, "KRIG", "TYPE")) == NULL) return 1;
    d->general->ktype = get_int(k);

    if ((k = get_key(r, "KRIG", "COLOCORR")) == NULL) return 1;
    d->general->colocorr = get_float(k);

    if ((k = get_key(r, "KRIG", "VARRED")) == NULL) return 1;
    d->general->varred = get_float(k);


    /* softdata parameters */
    if ((k = get_key(r, "SOFT", "NVARIL")) == NULL) return 1;
    d->general->nvaril = get_int(k);

    if ((k = get_key(r, "SOFT", "ICOLLVM")) == NULL) return 1;
    d->general->icollvm = get_int(k);


    /* variogram models */
    if ((k = get_key(r, "VARIOGRAM", "NUMBER")) == NULL) return 1;
    varnum = d->covariance->nst[0] = get_int(k);
    printf_dbg("load_dss_configs(): varnum = %d\n", varnum);

    if ((k = get_key(r, "VARIOGRAM", "NUGGET")) == NULL) return 1;
    sill = d->covariance->c0[0] = get_float(k);

    d->covariance->it = (int *) tsi_malloc(varnum * sizeof(int));
    d->covariance->cc = (float *) tsi_malloc(varnum * sizeof(float));
    d->covariance->aa = (float *) tsi_malloc(varnum * sizeof(float));
    d->covariance->ang1 = (float *) tsi_malloc(varnum * sizeof(float));
    d->covariance->ang2 = (float *) tsi_malloc(varnum * sizeof(float));
    d->covariance->ang3 = (float *) tsi_malloc(varnum * sizeof(float));
    d->covariance->anis1 = (float *) tsi_malloc(varnum * sizeof(float));
    d->covariance->anis2 = (float *) tsi_malloc(varnum * sizeof(float));
    if (!d->covariance->it ||
        !d->covariance->cc ||
        !d->covariance->aa ||
        !d->covariance->ang1 ||
        !d->covariance->ang2 ||
        !d->covariance->ang3 ||
        !d->covariance->anis1 ||
        !d->covariance->anis2) {
        printf("load_dss_configs(): ERROR - Failed to allocate space for variograms!\n");
        return 1;
    }

    for (i = 0; i < varnum; i++) {
        sprintf(variogram, "VARIOGRAM%d", i+1);
        if ((k = get_key(r, variogram, "TYPE")) == NULL) return 1;
        d->covariance->it[i] = get_int(k);
		printf_dbg2("%s - type: %d\t %d\n",variogram,get_int(k),d->covariance->it[i]);
        if (d->covariance->it[i] == 4) {
            printf("load_dss_configs(): ERROR - A power model is not allowed! Choose a different model and re-start.\n");
            return 1;
        }

        if ((k = get_key(r, variogram, "COV")) == NULL) return 1;
        d->covariance->cc[i] = get_float(k);
        sill += d->covariance->cc[i];

        if ((k = get_key(r, variogram, "ANG1")) == NULL) return 1;
        d->covariance->ang1[i] = get_float(k);

        if ((k = get_key(r, variogram, "ANG2")) == NULL) return 1;
        d->covariance->ang2[i] = get_float(k);

        if ((k = get_key(r, variogram, "ANG3")) == NULL) return 1;
        d->covariance->ang3[i] = get_float(k);

        if ((k = get_key(r, variogram, "AA")) == NULL) return 1;
        aa = d->covariance->aa[i] = get_float(k);
        if (aa < 1e-20)
			aa = 1e-20;

        if ((k = get_key(r, variogram, "AA1")) == NULL) return 1;
        aa1 = get_float(k);
        d->covariance->anis1[i] = aa1 / aa; 

        if ((k = get_key(r, variogram, "AA2")) == NULL) return 1;
        aa2 = get_float(k);
        d->covariance->anis2[i] = aa2 / aa;
		printf_dbg("load_dss_configs(): variogram %d\n \ttype:\t%d\n \tcov:\t%.2f\n \tang1:\t%.2f\n \tang2:\t%.2f\n \tang3:\t%.2f\n \taa:\t%.2f\n \taa1:\t%.2f\n \taa2:\t%.2f\n\n",
				i,d->covariance->it[i],d->covariance->cc[i], 
				d->covariance->ang1[i],d->covariance->ang2[i],d->covariance->ang3[i],aa,aa1,aa2);
    } /* for */

    if ((sill < 1) || (sill > 1)) {
        printf("load_dss_configs(): WARNING: The sill of the variogram is not 1.0! sill = %f\n", sill);
    }

    d->simulation->nsim = 1;   /* number of simulations */

	return 0;
}

