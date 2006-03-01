/* ----------------------------------------------------------------------- */
/*                      Quickersort Subroutine */
/*                      ********************** */

/* This is a subroutine for sorting a real array in ascending order. This */
/* is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen */
/* in collected algorithms of the ACM. */

/* The method used is that of continually splitting the array into parts */
/* such that all elements of one part are less than all elements of the */
/* other, with a third part in the middle consisting of one element.  An */
/* element with value t is chosen arbitrarily (here we choose the middle */
/* element). i and j give the lower and upper limits of the segment being */
/* split.  After the split a value q will have been found such that */
/* a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then */
/* performs operations on the two segments (i,q-1) and (q+1,j) as follows */
/* The smaller segment is split and the position of the larger segment is */
/* stored in the lt and ut arrays.  If the segment to be split contains */
/* two or fewer elements, it is sorted and another segment is obtained */
/* from the lt and ut arrays.  When no more segments remain, the array */
/* is completely sorted. */

/* INPUT PARAMETERS: */
/*   ib,ie        start and end index of the array to be sorted */
/*   a            array, a portion of which has to be sorted. */
/*   iperm        0 no other array is permuted. */
/*                1 array b is permuted according to array a */
/*                2 arrays b,c are permuted. */
/*                3 arrays b,c,d are permuted. */
/*                4 arrays b,c,d,e are permuted. */
/*                5 arrays b,c,d,e,f are permuted. */
/*                6 arrays b,c,d,e,f,g are permuted. */
/*                7 arrays b,c,d,e,f,g,h are permuted. */
/*               >7 no other array is permuted. */
/*   b,c,d,e,f,g,h  arrays to be permuted according to array a. */

/* OUTPUT PARAMETERS: */
/*    a      = the array, a portion of which has been sorted. */
/*    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm) */

/* NO EXTERNAL ROUTINES REQUIRED: */
/* ----------------------------------------------------------------------- */

/** funcoes utilizadas
 */

/** CUBOS utilizados
 *  nota: sao passados cubos por parametro para o sortem em outras funcoes
 */ 

/** CUBOS _nao_ utilizados
 */

/** structs globais utilizadas:
 */

#include <stdio.h>
#include "profile.h"

int sortem(int *ib,  int *ie, float *a, int * iperm, float *b, float *c__,
		float *d__, float *e, float *f, float *g, float *h__)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, j, k, m, p, q;
	float  ta, tb, tc, td, te, tf, tg, th, xa, xf, xg;
	int lt[64];
	float xh, xe, xd, xc, xb;
	int ut[64], iring;

#ifdef PROFILE
	profile.sortem++;
	profBegin("sortem");
#endif

	/* The dimensions for lt and ut have to be at least log (base 2) n */


	/* Initialize: */

	/* Parameter adjustments */
	--h__;
	--g;
	--f;
	--e;
	--d__;
	--c__;
	--b;
	--a;

	/* Function Body */
	j = *ie;
	m = 1;
	i__ = *ib;
	iring = *iperm + 1;
	if (*iperm > 7) {
		iring = 1;
	}

	/* If this segment has more than two elements  we split it */

L10:
	if ((i__1 = j - i__ - 1) < 0) {
		goto L100;
	} else if (i__1 == 0) {
		goto L90;
	} else {
		goto L15;
	}

	/* p is the position of an arbitrary element in the segment we choose the */
	/* middle element. Under certain circumstances it may be advantageous */
	/* to choose p at random. */

L15:
	p = (j + i__) / 2;
	ta = a[p];
	a[p] = a[i__];
	switch (iring) {
		case 1:  goto L21;
		case 2:  goto L19;
		case 3:  goto L18;
		case 4:  goto L17;
		case 5:  goto L16;
		case 6:  goto L161;
		case 7:  goto L162;
		case 8:  goto L163;
	}
L163:
	th = h__[p];
	h__[p] = h__[i__];
L162:
	tg = g[p];
	g[p] = g[i__];
L161:
	tf = f[p];
	f[p] = f[i__];
L16:
	te = e[p];
	e[p] = e[i__];
L17:
	td = d__[p];
	d__[p] = d__[i__];
L18:
	tc = c__[p];
	c__[p] = c__[i__];
L19:
	tb = b[p];
	b[p] = b[i__];
L21:

	/* Start at the beginning of the segment, search for k such that a(k)>t */

	q = j;
	k = i__;
L20:
	++k;
	if (k > q) {
		goto L60;
	}
	if (a[k] <= ta) {
		goto L20;
	}

	/* Such an element has now been found now search for a q such that a(q)<t */
	/* starting at the end of the segment. */

L30:
	if (a[q] < ta) {
		goto L40;
	}
	--q;
	if (q > k) {
		goto L30;
	}
	goto L50;

	/* a(q) has now been found. we interchange a(q) and a(k) */

L40:
	xa = a[k];
	a[k] = a[q];
	a[q] = xa;
	switch (iring) {
		case 1:  goto L45;
		case 2:  goto L44;
		case 3:  goto L43;
		case 4:  goto L42;
		case 5:  goto L41;
		case 6:  goto L411;
		case 7:  goto L412;
		case 8:  goto L413;
	}
L413:
	xh = h__[k];
	h__[k] = h__[q];
	h__[q] = xh;
L412:
	xg = g[k];
	g[k] = g[q];
	g[q] = xg;
L411:
	xf = f[k];
	f[k] = f[q];
	f[q] = xf;
L41:
	xe = e[k];
	e[k] = e[q];
	e[q] = xe;
L42:
	xd = d__[k];
	d__[k] = d__[q];
	d__[q] = xd;
L43:
	xc = c__[k];
	c__[k] = c__[q];
	c__[q] = xc;
L44:
	xb = b[k];
	b[k] = b[q];
	b[q] = xb;
L45:

	/* Update q and search for another pair to interchange: */

	--q;
	goto L20;
L50:
	q = k - 1;
L60:

	/* The upwards search has now met the downwards search: */

	a[i__] = a[q];
	a[q] = ta;
	switch (iring) {
		case 1:  goto L65;
		case 2:  goto L64;
		case 3:  goto L63;
		case 4:  goto L62;
		case 5:  goto L61;
		case 6:  goto L611;
		case 7:  goto L612;
		case 8:  goto L613;
	}
L613:
	h__[i__] = h__[q];
	h__[q] = th;
L612:
	g[i__] = g[q];
	g[q] = tg;
L611:
	f[i__] = f[q];
	f[q] = tf;
L61:
	e[i__] = e[q];
	e[q] = te;
L62:
	d__[i__] = d__[q];
	d__[q] = td;
L63:
	c__[i__] = c__[q];
	c__[q] = tc;
L64:
	b[i__] = b[q];
	b[q] = tb;
L65:

	/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
	/* store the position of the largest segment in lt and ut */

	if (q << 1 <= i__ + j) {
		goto L70;
	}
	lt[m - 1] = i__;
	ut[m - 1] = q - 1;
	i__ = q + 1;
	goto L80;
L70:
	lt[m - 1] = q + 1;
	ut[m - 1] = j;
	j = q - 1;

	/* Update m and split the new smaller segment */

L80:
	++m;
	goto L10;

	/* We arrive here if the segment has  two elements we test to see if */
	/* the segment is properly ordered if not, we perform an interchange */

L90:
	if (a[i__] <= a[j]) {
		goto L100;
	}
	xa = a[i__];
	a[i__] = a[j];
	a[j] = xa;
	switch (iring) {
		case 1:  goto L95;
		case 2:  goto L94;
		case 3:  goto L93;
		case 4:  goto L92;
		case 5:  goto L91;
		case 6:  goto L911;
		case 7:  goto L912;
		case 8:  goto L913;
	}
L913:
	xh = h__[i__];
	h__[i__] = h__[j];
	h__[j] = xh;
L912:
	xg = g[i__];
	g[i__] = g[j];
	g[j] = xg;
L911:
	xf = f[i__];
	f[i__] = f[j];
	f[j] = xf;
L91:
	xe = e[i__];
	e[i__] = e[j];
	e[j] = xe;
L92:
	xd = d__[i__];
	d__[i__] = d__[j];
	d__[j] = xd;
L93:
	xc = c__[i__];
	c__[i__] = c__[j];
	c__[j] = xc;
L94:
	xb = b[i__];
	b[i__] = b[j];
	b[j] = xb;
L95:

	/* If lt and ut contain more segments to be sorted repeat process: */

L100:
	--m;
	if (m <= 0) {
		goto L110;
	}
	i__ = lt[m - 1];
	j = ut[m - 1];
	goto L10;
L110:

#ifdef PROFILE
	profEnd("sortem");
#endif

	return 0;
} /* sortem_ */

