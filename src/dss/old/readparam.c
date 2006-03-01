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

#include <stdio.h>
#include <stdlib.h>

#include "dss.h"
#include "acorni.h"
#include "profile.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE_ (1)

int readparam(float *params, float *models, 
		general_vars_t *general,
		search_vars_t *search,
		simulation_vars_t *simulation,
		covariance_vars_t *covariance,
		covtable_lookup_vars_t *covt_lookup)
{
	/* System generated locals */
	int i__1;
	float r__1;


	/* Local variables */
	int nst1_pos__, seed_pos__, noct_pos__, tmin_pos__, tmax_pos__,
			   ivrl_pos__, zmin_pos__, zmax_pos__, isvr_pos__, iswt_pos__, 
			   ntry_pos__, xsiz_pos__, ysiz_pos__, zsiz_pos__, sang1_pos__, 
			   sang2_pos__, sang3_pos__, i__;
	double p;
	int icvar_pos__, ltail_pos__, imask_pos__, ndmin_pos__, 
			   ndmax_pos__, utail_pos__, nvari_pos__, ltpar_pos__, utpar_pos__, 
			   ktype_pos__, nmult_pos__, mults_pos__;
	float aa1, aa2;
	int icmean_pos__, varred_pos__, nodmax_pos__, radius_pos__, 
			   isecvr_pos__, nvaril_pos__, itrans_pos__, sstrat_pos__, 
			   radius1_pos__, radius2_pos__, icollvm_pos__, ismooth_pos__;
	float sill;
	int colocorr_pos__, nosvalue_pos__, c01_pos__;
	int testfl;
	int nx_pos__, ny_pos__, nz_pos__;
	float radius1, radius2;
	int ixl_pos__, iyl_pos__, izl_pos__, xmn_pos__, iwt_pos__, 
			   ymn_pos__, zmn_pos__;


#ifdef PROFILE
	profile.readparm++;
	profBegin("readparm");
#endif

	/* Parameter adjustments */
	--models;
	--params;

	/* Function Body */
	nvari_pos__ = 1;
	ixl_pos__ = 2;
	iyl_pos__ = 3;
	izl_pos__ = 4;
	ivrl_pos__ = 5;
	iwt_pos__ = 7;
	isecvr_pos__ = 7;
	tmin_pos__ = 8;
	tmax_pos__ = 9;
	itrans_pos__ = 10;
	ismooth_pos__ = 11;
	isvr_pos__ = 12;
	iswt_pos__ = 13;
	zmin_pos__ = 14;
	zmax_pos__ = 15;
	ltail_pos__ = 16;
	ltpar_pos__ = 17;
	utail_pos__ = 18;
	utpar_pos__ = 19;
	ntry_pos__ = 20;
	icmean_pos__ = 21;
	icvar_pos__ = 22;
	nx_pos__ = 23;
	xmn_pos__ = 24;
	xsiz_pos__ = 25;
	ny_pos__ = 26;
	ymn_pos__ = 27;
	ysiz_pos__ = 28;
	nz_pos__ = 29;
	zmn_pos__ = 30;
	zsiz_pos__ = 31;
	nosvalue_pos__ = 32;
	imask_pos__ = 33;
	seed_pos__ = 34;
	ndmin_pos__ = 35;
	ndmax_pos__ = 36;
	nodmax_pos__ = 37;
	sstrat_pos__ = 38;
	mults_pos__ = 39;
	nmult_pos__ = 40;
	noct_pos__ = 41;
	radius_pos__ = 42;
	radius1_pos__ = 43;
	radius2_pos__ = 44;
	sang1_pos__ = 45;
	sang2_pos__ = 46;
	sang3_pos__ = 47;
	ktype_pos__ = 48;
	colocorr_pos__ = 49;
	varred_pos__ = 50;
	nvaril_pos__ = 51;
	icollvm_pos__ = 52;
	nst1_pos__ = 53;
	c01_pos__ = 54;
	general->nvari = params[nvari_pos__];
	general->ixl = params[ixl_pos__];
	general->iyl = params[iyl_pos__];
	general->izl = params[izl_pos__];
	general->ivrl = params[ivrl_pos__];
	general->iwt = params[iwt_pos__];
	general->isecvr = params[isecvr_pos__];
	general->tmin = params[tmin_pos__];
	general->tmax = params[tmax_pos__];
	general->itrans = params[itrans_pos__];
	/* !	read(lin,'(a256)',err=98) transfl */
	general->ismooth = params[ismooth_pos__];
	/* !      read(lin,'(a256)',err=98) smthfl */
	general->isvr = params[isvr_pos__];
	general->iswt = params[iswt_pos__];
	general->zmin = params[zmin_pos__];
	general->zmax = params[zmax_pos__];
	general->ltail = params[ltail_pos__];
	general->ltpar = params[ltpar_pos__];
	general->utail = params[utail_pos__];
	general->utpar = params[utpar_pos__];
	/* !MARTELADO para só fazer 1 simulação */
	simulation->nsim = 1;
	covt_lookup->ntry = params[ntry_pos__];
	covt_lookup->icmean = params[icmean_pos__];
	covt_lookup->icvar = params[icvar_pos__];
	general->nx = params[nx_pos__];
	general->xmn = params[xmn_pos__];
	general->xsiz = params[xsiz_pos__];
	general->ny = params[ny_pos__];
	general->ymn = params[ymn_pos__];
	general->ysiz = params[ysiz_pos__];
	general->nz = params[nz_pos__];
	general->zmn = params[zmn_pos__];
	general->zsiz = params[zsiz_pos__];
	general->nxy = general->nx * general->ny;
	general->nxyz = general->nx * general->ny * general->nz;
	general->nosvalue = params[nosvalue_pos__];
	general->imask = params[imask_pos__];
	/*	iaco_1.ixv[0] = params[seed_pos__];*/

	int rand = random();
	newAcorni(rand);
	for (i__ = 1; i__ <= 1000; ++i__) {
		p = acorni();
	}
	search->ndmin = params[ndmin_pos__];
	search->ndmax = params[ndmax_pos__];
	covt_lookup->nodmax = params[nodmax_pos__];
	search->sstrat = params[sstrat_pos__];
	if (search->sstrat == 1) {
		search->ndmax = 0.f;
	}
	search->mults = params[mults_pos__];
	search->nmult = params[nmult_pos__];
	search->noct = params[noct_pos__];
	search->radius = params[radius_pos__];
	radius1 = params[radius1_pos__];
	radius2 = params[radius2_pos__];
	if (search->radius < 1e-20f) {
		fprintf(stderr,"radius (%f), must be greater than EPSLON (%f)\n",search->radius, 1e-20f);
		return -1; /* ERROR */
	}
	search->radsqd = search->radius * search->radius;
	search->sanis1 = radius1 / search->radius;
	search->sanis2 = radius2 / search->radius;
	search->sang1 = params[sang1_pos__];
	search->sang2 = params[sang2_pos__];
	search->sang3 = params[sang3_pos__];
	general->ktype = params[ktype_pos__];
	general->colocorr = params[colocorr_pos__];
	general->varred = params[varred_pos__];
	/* !      read(lin,'(a256)',err=98) corrfl */
	/* !      read(lin,'(a256)',err=98) lvmfl */
	general->nvaril = params[nvaril_pos__];
	general->icollvm = params[icollvm_pos__];
	covariance->nst[0] = params[nst1_pos__];
	covariance->c0[0] = params[c01_pos__];
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
		covariance->anis1[i__ - 1] = aa1 / ((double) MAX(r__1,1e-20f));  //dmax(r__1,1e-20f);
		r__1 = covariance->aa[i__ - 1];
		covariance->anis2[i__ - 1] = aa2 / ((double) MAX(r__1,1e-20f)); //dmax(r__1,1e-20f);
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
		testfl = TRUE_;
	}
	if (general->utail != 1 && general->utail != 2 && general->utail != 4) {
		fprintf(stderr,"ERROR invalid upper tail option: %d\n",(int) general->utail);
		fprintf(stderr,"\t\tonly allowed 1, 2 or 4 - see manual.\n");
		testfl = TRUE_;
	}
	if (general->utail == 4 && general->utpar < 1.f) {
		fprintf(stderr,"ERROR invalid power for hyperbolic tail: %f\n",general->utpar);
		fprintf(stderr,"\t\tmust be greater than 1.0.\n");
		testfl = TRUE_;
	}
	if (general->ltail == 2 && general->ltpar < 0.f) {
		fprintf(stderr,"ERROR invalid power for power model: %f\n",general->ltpar);
		fprintf(stderr,"\t\tmust be greater than 0.0.\n");
		testfl = TRUE_;
	}
	if (general->utail == 2 && general->utpar < 0.f) {
		fprintf(stderr,"ERROR invalid power for power model: %f\n",general->utpar);
		fprintf(stderr,"\t\tmust be greater than 0.0.\n");
		testfl = TRUE_;
	}

#ifdef PROFILE
	profEnd("readparm");
#endif

	return 0;
} /* readparm_ */

