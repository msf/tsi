#include <math.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)




/* ----------------------------------------------------------------------- */
/*                Solution of a System of Linear Equations */
/*                **************************************** */

/* INPUT VARIABLES: */

/*   nright,nsb       number of columns in right hand side matrix. */
/*                      for KB2D: nright=1, nsb=1 */
/*   neq              number of equations */
/*   a()              upper triangular left hand side matrix (stored */
/*                      columnwise) */
/*   r()              right hand side matrix (stored columnwise) */
/*                      for kb2d, one column per variable */
/* OUTPUT VARIABLES: */

/*   s()              solution array, same dimension as  r  above. */
/*   ising            singularity indicator */
/*                      0,  no singularity problem */
/*                     -1,  neq .le. 1 */
/*                      k,  a null pivot appeared at the kth iteration */
/* PROGRAM NOTES: */

/*   1. Requires the upper triangular left hand side matrix. */
/*   2. Pivots are on the diagonal. */
/*   3. Does not search for max. element for pivot. */
/*   4. Several right hand side matrices possible. */
/*   5. USE for ok and sk only, NOT for UK. */
/* ----------------------------------------------------------------------- */

/**
 * funcoes utilizadas
 * -
 */

/** CUBOS utilizados
 * -
 */ 

/** CUBOS _nao_ utilizados
 * sim
 * tmp
 * order
 * clc
 * lvm
 * e os cov*
 */

/** structs globais utilizadas:
 * krigev_1 (krigev_1.{a,r__,s})
 */




int ksol(int *nright, int *neq, int *nsb, double *a, double *r, double *s,
		int *ising)
{
	/* System generated locals */
	int i__1;
	double d__1;

	/* Local variables */
	int i__, j, k, m1;
	double ak, ap;
	int ii, ij = 0, kk, in, ll, nm, nn, lp, iv, km1, ll1, nm1, llb, 
			   ijm;
	double tol, piv;


	/* If there is only one equation then set ising and return: */

	/* Parameter adjustments */
	--s;
	--r;
	--a;

	/* Function Body */
	if (*neq <= 1) {
		*ising = -1;
		return 0;
	}

	/* Initialize: */

	tol = 1e-7f;
	*ising = 0;
	nn = *neq * (*neq + 1) / 2;
	nm = *nsb * *neq;
	m1 = *neq - 1;
	kk = 0;

	/* Start triangulation: */

	for (k = 1; k <= m1; ++k) {
		kk += k;
		ak = a[kk];
		if (fabs(ak) < tol) {
			*ising = k;
			return 0;
		}
		km1 = k - 1;
		for (iv = 1; iv <= *nright; ++iv) {
			nm1 = nm * (iv - 1);
			ii = kk + nn * (iv - 1);
			piv = 1.f / a[ii];
			lp = 0;
			for (i__ = k; i__ <= m1; ++i__) {
				ll = ii;
				ii += i__;
				ap = a[ii] * piv;
				++lp;
				ij = ii - km1;
				for (j = i__; j <= m1; ++j) {
					ij += j;
					ll += j;
					a[ij] -= ap * a[ll];
				}
				/* optimizacoes. */
				i__1 = lp + nm1;
				/* fim das opts */
				for (llb = k; llb <= nm; llb += *neq) {
					in = llb + i__1;
					ll1 = llb + nm1;
					r[in] -= ap * r[ll1];
				}
			}
		}
	}

	/* Error checking - singular matrix: */

	ijm = ij - nn * (*nright - 1);
	if ((d__1 = a[ijm], fabs(d__1)) < tol) {
		*ising = *neq;
		return 0;
	}

	/* Finished triangulation, start solving back: */

	for (iv = 1; iv <= *nright; ++iv) {
		nm1 = nm * (iv - 1);
		ij = ijm + nn * (iv - 1);
		piv = 1.f / a[ij];
		for (llb = *neq; llb <= nm; llb += *neq) {
			ll1 = llb + nm1;
			s[ll1] = r[ll1] * piv;
		}
		i__ = *neq;
		kk = ij;
		for (ii = 1; ii <= m1; ++ii) {
			kk -= i__;
			piv = 1.f / a[kk];
			--i__;
			for (llb = i__; llb <= nm; llb += *neq) {
				ll1 = llb + nm1;
				in = ll1;
				ap = r[in];
				ij = kk;
				for (j = i__; j <= m1; ++j) {
					ij += j;
					++in;
					ap -= a[ij] * s[in];
				}
				s[ll1] = ap * piv;
			}
		}
	}

	/* Finished solving back, return: */

	return 0;
} /* ksol_ */

