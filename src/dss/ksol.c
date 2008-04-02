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

int ksol(int neq, double *a, double *r, double *s)
{
	/* System generated locals */
	int i1;
	double d1;

	int nright, nsb;

	/* Local variables */
	int i, j, k, m1, ising;
	double ak, ap;
	int ii, ij, kk, in, ll, nm, nn, lp, iv, km1, ll1, nm1, llb, ijm;
	double tol, piv;

	nright = 1;
	nsb = 1;

	/* If there is only one equation then set ising and return: */

	/* Function Body */
	if (neq <= 1) {
		ising = -1;
		return ising;
	}

	/* Initialize: */

	ij = 0;
	tol = 1e-7f;
	ising = 0;
	nn = neq * (neq + 1) / 2;
	nm = nsb * neq;
	m1 = neq - 1;
	kk = 0;

	/* Start triangulation: */

	for (k = 1; k <= m1; ++k) {
		kk += k;
		ak = a[kk];
		if (fabs(ak) < tol) {
			ising = k;
			return ising;
		}
		km1 = k - 1;
		for (iv = 0; iv < nright; ++iv) {
			nm1 = nm * iv;
			ii = kk + nn * iv;
			piv = 1.f / a[ii];
			lp = 0;
			for (i = k; i <= m1; ++i) {
				ll = ii;
				ii += i;
				ap = a[ii] * piv;
				++lp;
				ij = ii - km1;
				for (j = i; j <= m1; ++j) {
					ij += j;
					ll += j;
					a[ij] -= ap * a[ll];
				}
				/* optimizacoes. */
				i1 = lp + nm1;
				/* fim das opts */
				for (llb = k; llb <= nm; llb += neq) {
					in = llb + i1;
					ll1 = llb + nm1;
					r[in] -= ap * r[ll1];
				}
			}
		}
	}

	/* Error checking - singular matrix: */

	ijm = ij - nn * (nright - 1);
	if ((d1 = a[ijm], fabs(d1)) < tol) {
		ising = neq;
		return ising;
	}

	/* Finished triangulation, start solving back: */

	for (iv = 0; iv < nright; ++iv) {
		nm1 = nm * iv;
		ij = ijm + nn * iv;
		piv = 1.f / a[ij];
		for (llb = neq; llb <= nm; llb += neq) {
			ll1 = llb + nm1;
			s[ll1] = r[ll1] * piv;
		}
		i = neq;
		kk = ij;
		for (ii = 1; ii <= m1; ++ii) {
			kk -= i;
			piv = 1.f / a[kk];
			--i;
			for (llb = i; llb <= nm; llb += neq) {
				ll1 = llb + nm1;
				in = ll1;
				ap = r[in];
				ij = kk;
				for (j = i; j <= m1; ++j) {
					ij += j;
					++in;
					ap -= a[ij] * s[in];
				}
				s[ll1] = ap * piv;
			}
		}
	}

	/* Finished solving back, return: */

	return ising;
} /* ksol_ */

