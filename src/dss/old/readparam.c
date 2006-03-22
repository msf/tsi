#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"

#undef PROFILE

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)



/* !----------------------------------------------------------------------- */
/* ! Initialization and Read Parameters */
/* ! */
/* ! The input parameters and data are read in from their files. Some quick */
/* ! error checking is performed and the statistics of all the variables */
/* ! being considered are written to standard output. */
/* !----------------------------------------------------------------------- */

/** Funcoes utilizadas
 */

/** CUBOS utilizados
 * ??
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 * generl_1
 * simula_1
 * clooku_1
 * search_1
 * cova3d_1
 */


int readparam(float *params, float *models, 
		general_vars_t *general,
		search_vars_t *search,
		simulation_vars_t *simulation,
		covariance_vars_t *covariance,
		covtable_lookup_vars_t *covt_lookup)
{
	/* System generated locals */
	int i__1, i__;
	float r__1;


	/* Local variables */
	float aa1, aa2;
	float sill;
	int testfl;
	float radius1, radius2;


#ifdef PROFILE
	profile.readparm++;
	profBegin("readparm");
#endif

	/* Parameter adjustments */
	--models;
	--params;

	/* Function Body */
	general->nvari = params[NVARI+1];
	general->ixl = params[IXL+1];
	general->iyl = params[IYL+1];
	general->izl = params[IZL+1];
	general->ivrl = params[IVRL+1];
	general->iwt = params[IWT+1];
	general->isecvr = params[ISECVR+1];
	general->tmin = params[TMIN+1];
	general->tmax = params[TMAX+1];
	general->itrans = params[ITRANS+1];
	general->ismooth = params[ISMOOTH+1];
	general->isvr = params[ISVR+1];
	general->iswt = params[ISWT+1];
	general->zmin = params[ZMIN+1];
	general->zmax = params[ZMAX+1];
	general->ltail = params[LTAIL+1];
	general->ltpar = params[LTPAR+1];
	general->utail = params[UTAIL+1];
	general->utpar = params[UTPAR+1];
	general->nx = params[NX+1];
	general->xmn = params[XMN+1];
	general->xsiz = params[XSIZ+1];
	general->ny = params[NY+1];
	general->ymn = params[YMN+1];
	general->ysiz = params[YSIZ+1];
	general->nz = params[NZ+1];
	general->zmn = params[ZMN+1];
	general->zsiz = params[ZSIZ+1];
	general->nxy = general->nx * general->ny;
	general->nxyz = general->nx * general->ny * general->nz;
	general->nosvalue = params[NOSVALUE+1];
	general->imask = params[IMASK+1];
	general->ktype = params[KTYPE+1];
	general->colocorr = params[COLOCORR+1];
	general->varred = params[VARRED+1];
	general->nvaril = params[NVARIL+1];
	general->icollvm = params[ICOLLVM+1];

	/* !MARTELADO para só fazer 1 simulação */
	simulation->nsim = 1;

	covt_lookup->ntry = params[NTRY+1];
	covt_lookup->icmean = params[ICMEAN+1];
	covt_lookup->icvar = params[ICVAR+1];
	covt_lookup->nodmax = params[NODMAX+1];
	
	search->ndmin = (int) params[NDMIN+1];   /* LPL */
	search->ndmax = (int) params[NDMAX+1];
	/*search->ndmin = params[NDMIN+1];
	search->ndmax = params[NDMAX+1]; */
	search->sstrat = params[SSTRAT+1];
	if (search->sstrat == 1) {
		search->ndmax = 0;
	}
	search->mults = params[MULTS+1];
	search->nmult = params[NMULTS+1];
	/*search->noct =  params[NOCT+1]; */          /* LPL */
	search->noct = (int) params[NOCT+1];
	search->radius = params[RADIUS+1];
	radius1 = params[RADIUS1+1];
	radius2 = params[RADIUS2+1];
	if (search->radius < 1e-20f) {
		fprintf(stderr,"radius (%f), must be greater than EPSLON (%f)\n",search->radius, 1e-20f);
		return -1; /* ERROR */
	}
	search->radsqd = search->radius * search->radius;
	search->sanis1 = radius1 / search->radius;
	search->sanis2 = radius2 / search->radius;
	search->sang1 = params[SANG+1];
	search->sang2 = params[SANG1+1];
	search->sang3 = params[SANG2+1];

	covariance->nst[0] = params[VARNUM+1];
	covariance->c0[0] = params[NUGGET+1];
	sill = covariance->c0[0];
	i__1 = covariance->nst[0];
	for (i__ = 1; i__ <= i__1; ++i__) {
		covariance->it[i__ - 1] = models[(i__ << 3) - 7];
		covariance->cc[i__ - 1] = models[(i__ << 3) - 6];
		covariance->ang1[i__ - 1] = models[(i__ << 3) - 5];
		covariance->ang2[i__ - 1] = models[(i__ << 3) - 4];
		covariance->ang3[i__ - 1] = models[(i__ << 3) - 3];
		covariance->aa[i__ - 1] = models[(i__ << 3) - 2];
		aa1 = models[(i__ << 3) - 1];
		aa2 = models[i__ * 8];
		r__1 = covariance->aa[i__ - 1];
		covariance->anis1[i__ - 1] = aa1 / ((double) MAX(r__1,1e-20f));  /*dmax(r__1,1e-20f)*/;
		r__1 = covariance->aa[i__ - 1];
		covariance->anis2[i__ - 1] = aa2 / ((double) MAX(r__1,1e-20f)); /*dmax(r__1,1e-20f) */;
		sill += covariance->cc[i__ - 1];
		if (covariance->it[i__ - 1] == 4) {
			fprintf(stderr,"A power model is NOT allowedi\n");
			fprintf(stderr,"Choose a diferente model and re-start\n");
			return -1; /* ERROR */
		}
	}
	/* Warn the user if the sill is different than 1.0: */
	if (sill > 1.f || sill < 1.f) {
		fprintf(stderr,"WARNING the sill of your variogram is not 1.0\n");
		fprintf(stderr,"\t\tsill = %f\n",sill);
	}
	if (general->ltail != 1 && general->ltail != 2) {
		fprintf(stderr,"ERROR invalid lower tail option: %d\n",(int)general->ltail);
		fprintf(stderr,"\t\tonly allowed 1 or 2 - see manual.\n");
		testfl = TRUE;
	}
	if (general->utail != 1 && general->utail != 2 && general->utail != 4) {
		fprintf(stderr,"ERROR invalid upper tail option: %d\n",(int) general->utail);
		fprintf(stderr,"\t\tonly allowed 1, 2 or 4 - see manual.\n");
		testfl = TRUE;
	}
	if (general->utail == 4 && general->utpar < 1.f) {
		fprintf(stderr,"ERROR invalid power for hyperbolic tail: %f\n",general->utpar);
		fprintf(stderr,"\t\tmust be greater than 1.0.\n");
		testfl = TRUE;
	}
	if (general->ltail == 2 && general->ltpar < 0.f) {
		fprintf(stderr,"ERROR invalid power for power model: %f\n",general->ltpar);
		fprintf(stderr,"\t\tmust be greater than 0.0.\n");
		testfl = TRUE;
	}
	if (general->utail == 2 && general->utpar < 0.f) {
		fprintf(stderr,"ERROR invalid power for power model: %f\n",general->utpar);
		fprintf(stderr,"\t\tmust be greater than 0.0.\n");
		testfl = TRUE;
	}

#ifdef PROFILE
	profEnd("readparm");
#endif

	return 0;
} /* readparm_ */

